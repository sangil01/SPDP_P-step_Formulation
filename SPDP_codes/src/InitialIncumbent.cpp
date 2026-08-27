#include "InitialIncumbent.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <utility>

#include "CpIncumbentAdapter.h"
#include "TwoIndexModelCore.h"
#include "gurobi_c++.h"

namespace spdp {

using detail::TwoIndexCoreModel;
using detail::TwoIndexCoreObjective;
using detail::TwoIndexCoreOptions;
using detail::build_two_index_model_core;

namespace {

struct InitialIncumbentSearchState {
    int initial_k = 0;
    int certified_k = 0;
    int candidate_k = 0;
    int increments_used = 0;
    int max_k_increments = 0;
    int max_candidate_k = 0;
    InitialIncumbentTimeoutAction timeout_action =
        InitialIncumbentTimeoutAction::Stop;
    bool finished = false;
    bool solution_found = false;
    bool stopped_on_time_limit = false;
    bool stopped_on_early_unknown = false;
    int early_unknown_k = 0;
    bool exhausted_increment_limit = false;
    bool exhausted_candidate_limit = false;
};

InitialIncumbentSearchState make_initial_incumbent_search_state(
    int initial_k,
    int max_k_increments,
    int max_candidate_k,
    InitialIncumbentTimeoutAction timeout_action
) {
    if (initial_k <= 0) {
        throw std::runtime_error(
            "Initial incumbent search requires a positive initial vehicle count."
        );
    }
    if (max_k_increments < 0) {
        throw std::runtime_error(
            "Initial incumbent search requires nonnegative maximum increments."
        );
    }
    if (max_candidate_k < initial_k) {
        throw std::runtime_error(
            "Initial incumbent search candidate limit is below the initial vehicle count."
        );
    }

    InitialIncumbentSearchState state;
    state.initial_k = initial_k;
    state.certified_k = initial_k;
    state.candidate_k = initial_k;
    state.max_k_increments = max_k_increments;
    state.max_candidate_k = max_candidate_k;
    state.timeout_action = timeout_action;
    return state;
}

InitialIncumbentSearchState advance_initial_incumbent_search(
    const InitialIncumbentSearchState& state,
    const InitialIncumbentSolveResult& solve_result
) {
    if (state.finished) {
        throw std::runtime_error(
            "Cannot advance a finished initial incumbent search."
        );
    }

    InitialIncumbentSearchState next = state;
    const FixedKSolveOutcome outcome = solve_result.outcome;
    if (outcome == FixedKSolveOutcome::Feasible) {
        next.finished = true;
        next.solution_found = true;
        return next;
    }

    if (outcome == FixedKSolveOutcome::ModelInvalid ||
        outcome == FixedKSolveOutcome::AdapterError) {
        next.finished = true;
        return next;
    }

    if (outcome == FixedKSolveOutcome::EarlyUnknown) {
        next.finished = true;
        next.stopped_on_early_unknown = true;
        next.early_unknown_k = state.candidate_k;
        return next;
    }

    if (outcome == FixedKSolveOutcome::TimedOutUnknown &&
        state.timeout_action == InitialIncumbentTimeoutAction::Stop) {
        next.finished = true;
        next.stopped_on_time_limit = true;
        return next;
    }

    // VI44 can be strengthened only by a contiguous infeasibility proof.
    // Advancing after a time limit leaves the certified lower bound unchanged.
    if (outcome == FixedKSolveOutcome::ProvenInfeasible &&
        state.candidate_k == state.certified_k) {
        next.certified_k = state.candidate_k + 1;
    }

    if (state.increments_used >= state.max_k_increments) {
        next.finished = true;
        next.exhausted_increment_limit = true;
        return next;
    }
    if (state.candidate_k >= state.max_candidate_k) {
        next.finished = true;
        next.exhausted_candidate_limit = true;
        return next;
    }

    next.candidate_k = state.candidate_k + 1;
    next.increments_used = state.increments_used + 1;
    return next;
}

// Match the 1e-9 feasibility tolerance used by the two-index core and route
// validator. A larger scaled tolerance could accept a route that validation rejects.
constexpr double kMakespanThresholdTolerance = 1e-9;

class MakespanThresholdCallback final : public GRBCallback {
public:
    MakespanThresholdCallback(double threshold, double tolerance)
        : threshold_(threshold), tolerance_(tolerance) {}

    bool stopped_by_feasible_callback() const {
        return stopped_by_feasible_callback_;
    }

    bool stopped_by_bound_callback() const {
        return stopped_by_bound_callback_;
    }

    double proving_bound() const {
        return proving_bound_;
    }

protected:
    void callback() override {
        try {
            if (where == GRB_CB_MIPSOL) {
                const double incumbent = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
                if (std::isfinite(incumbent) &&
                    incumbent <= threshold_ + tolerance_) {
                    stopped_by_feasible_callback_ = true;
                    abort();
                }
                return;
            }
            if (where == GRB_CB_MIP) {
                const double bound = getDoubleInfo(GRB_CB_MIP_OBJBND);
                if (std::isfinite(bound) &&
                    std::abs(bound) < 0.5 * GRB_INFINITY &&
                    bound > threshold_ + tolerance_) {
                    stopped_by_bound_callback_ = true;
                    proving_bound_ = bound;
                    abort();
                }
            }
        } catch (const GRBException&) {
            // Let Gurobi continue. Post-solve attributes still provide a safe
            // classification if callback information is unavailable.
        }
    }

private:
    double threshold_ = 0.0;
    double tolerance_ = 0.0;
    bool stopped_by_feasible_callback_ = false;
    bool stopped_by_bound_callback_ = false;
    double proving_bound_ = 0.0;
};

void add_fixed_vehicle_count_constraint(
    TwoIndexCoreModel& core,
    const MultiDiGraph& graph,
    int vehicle_count
) {
    GRBLinExpr actual_departures = 0.0;
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        if (edge.u == 0 &&
            graph.node(edge.v).kind == NodeSpec::Kind::Pickup) {
            actual_departures += core.y_vars[edge_id];
        }
    }
    core.model->addConstr(
        actual_departures == static_cast<double>(vehicle_count),
        "initial_incumbent_fixed_vehicle_count"
    );
}

void record_milp_objective_metadata(
    GRBModel& model,
    InitialIncumbentSolveResult& result
) {
    const int solution_count = model.get(GRB_IntAttr_SolCount);
    if (solution_count > 0) {
        const double objective_value = model.get(GRB_DoubleAttr_ObjVal);
        if (std::isfinite(objective_value)) {
            result.milp_has_objective_value = true;
            result.milp_objective_value = objective_value;
        }
    }
    const double objective_bound = model.get(GRB_DoubleAttr_ObjBound);
    if (std::isfinite(objective_bound) &&
        std::abs(objective_bound) < 0.5 * GRB_INFINITY) {
        result.milp_has_objective_bound = true;
        result.milp_best_objective_bound = objective_bound;
    }
}

void extract_milp_incumbent(
    const MultiDiGraph& graph,
    const TwoIndexCoreModel& core,
    InitialIncumbentSolveResult& result
) {
    result.total_duration = 0.0;
    result.total_original_cost = 0.0;
    result.edge_values.resize(graph.number_of_edges(), 0.0);
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const double value = core.y_vars[edge_id].get(GRB_DoubleAttr_X);
        result.edge_values[edge_id] = value;
        result.total_duration += graph.edges()[edge_id].data.time * value;
        result.total_original_cost += graph.edges()[edge_id].data.cost * value;
        const EdgeRecord& edge = graph.edges()[edge_id];
        if (value > 0.5 && edge.u == 0 &&
            graph.node(edge.v).kind == NodeSpec::Kind::Pickup) {
            ++result.vehicle_count;
        }
        if (value > 0.5) {
            result.active_edge_ids.push_back(static_cast<int>(edge_id));
        }
    }
}

void validate_initial_incumbent_solve_options(
    const InitialIncumbentSolveOptions& options
) {
    if (options.vehicle_count <= 0) {
        throw std::runtime_error(
            "Initial-incumbent generation requires a positive fixed vehicle count."
        );
    }
    if (!std::isfinite(options.solver_time_limit) ||
        options.solver_time_limit < 0.0) {
        throw std::runtime_error(
            "Initial-incumbent time limit must be finite and nonnegative."
        );
    }
    if (options.milp_mode == InitialIncumbentMilpMode::Makespan &&
        (!std::isfinite(options.milp_makespan_horizon_factor) ||
         options.milp_makespan_horizon_factor <= 1.0)) {
        throw std::runtime_error(
            "The MILP makespan horizon factor must be finite and greater than one."
        );
    }
}

InitialIncumbentSolveResult make_base_milp_result(
    const InitialIncumbentSolveOptions& options,
    GRBModel& model
) {
    InitialIncumbentSolveResult result;
    result.backend = InitialIncumbentBackend::TwoIndexMilp;
    result.milp_mode = options.milp_mode;
    result.milp_makespan_horizon_factor =
        options.milp_makespan_horizon_factor;
    result.status = model.get(GRB_IntAttr_Status);
    result.raw_status = result.status;
    result.status_name = std::to_string(result.status);
    result.hit_time_limit = result.status == GRB_TIME_LIMIT;
    result.hit_solution_limit = result.status == GRB_SOLUTION_LIMIT;
    result.runtime_seconds = model.get(GRB_DoubleAttr_Runtime);
    result.configured_time_limit_seconds = options.solver_time_limit;
    record_milp_objective_metadata(model, result);
    return result;
}

InitialIncumbentSolveResult solve_fixed_k_duration_impl(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const InitialIncumbentSolveOptions& options
) {
    TwoIndexCoreOptions core_options;
    core_options.objective = TwoIndexCoreObjective::Duration;
    core_options.binary_y = true;
    core_options.add_time_constraints = true;
    core_options.solver_time_limit = options.solver_time_limit;
    core_options.gurobi_threads = options.gurobi_threads;
    core_options.output_enabled = true;
    core_options.gurobi_log_path = options.gurobi_log_path;
    core_options.name_prefix = "initial_incumbent";
    TwoIndexCoreModel core = build_two_index_model_core(data, graph, core_options);

    add_fixed_vehicle_count_constraint(core, graph, options.vehicle_count);
    core.model->set(GRB_IntParam_SolutionLimit, 1);
    core.model->update();
    core.model->optimize();

    InitialIncumbentSolveResult result =
        make_base_milp_result(options, *core.model);
    result.milp_model_horizon = data.time_limit;
    result.has_feasible_solution = core.model->get(GRB_IntAttr_SolCount) > 0;
    if (result.has_feasible_solution) {
        result.outcome = FixedKSolveOutcome::Feasible;
        result.termination_name = "feasible";
    } else if (result.status == GRB_INFEASIBLE) {
        result.outcome = FixedKSolveOutcome::ProvenInfeasible;
        result.termination_name = "proven-infeasible";
        result.infeasible = true;
    } else if (result.hit_time_limit) {
        result.outcome = FixedKSolveOutcome::TimedOutUnknown;
        result.termination_name = "timed-out-unknown";
    } else {
        result.outcome = FixedKSolveOutcome::EarlyUnknown;
        result.termination_name = "early-unknown";
        result.early_unknown = true;
    }
    if (result.has_feasible_solution) {
        extract_milp_incumbent(graph, core, result);
    }
    return result;
}

InitialIncumbentSolveResult solve_fixed_k_makespan_impl(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const InitialIncumbentSolveOptions& options
) {
    const double threshold = data.time_limit;
    const double horizon = std::ceil(
        options.milp_makespan_horizon_factor * threshold
    );
    const double tolerance = kMakespanThresholdTolerance;

    TwoIndexCoreOptions core_options;
    core_options.objective = TwoIndexCoreObjective::None;
    core_options.binary_y = true;
    core_options.add_time_constraints = true;
    core_options.time_horizon = horizon;
    core_options.enforce_route_duration_limit = false;
    core_options.solver_time_limit = options.solver_time_limit;
    core_options.gurobi_threads = options.gurobi_threads;
    core_options.output_enabled = true;
    core_options.gurobi_log_path = options.gurobi_log_path;
    core_options.name_prefix = "initial_incumbent_makespan";
    TwoIndexCoreModel core = build_two_index_model_core(data, graph, core_options);

    add_fixed_vehicle_count_constraint(core, graph, options.vehicle_count);
    GRBVar makespan = core.model->addVar(
        0.0,
        horizon,
        1.0,
        GRB_CONTINUOUS,
        "initial_incumbent_makespan_M"
    );
    const NodeId end_node_id = graph.end_node_id();
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        if (edge.u == 0 || edge.v != end_node_id) {
            continue;
        }
        const GRBVar start_time = detail::two_index_time_var(
            core, edge.u, edge.data.start_state
        );
        core.model->addConstr(
            start_time + edge.data.time * core.y_vars[edge_id] <=
                makespan + horizon * (1.0 - core.y_vars[edge_id]),
            "initial_incumbent_makespan_return_" + std::to_string(edge_id)
        );
    }
    core.model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    core.model->update();

    MakespanThresholdCallback callback(threshold, tolerance);
    core.model->setCallback(&callback);
    core.model->optimize();

    InitialIncumbentSolveResult result =
        make_base_milp_result(options, *core.model);
    result.milp_model_horizon = horizon;
    result.milp_stopped_by_feasible_callback =
        callback.stopped_by_feasible_callback();
    result.milp_stopped_by_bound_callback =
        callback.stopped_by_bound_callback();
    if (callback.stopped_by_bound_callback() &&
        (!result.milp_has_objective_bound ||
         callback.proving_bound() > result.milp_best_objective_bound)) {
        result.milp_has_objective_bound = true;
        result.milp_best_objective_bound = callback.proving_bound();
    }

    const bool threshold_feasible =
        result.milp_has_objective_value &&
        result.milp_objective_value <= threshold + tolerance;
    const bool threshold_infeasible_by_bound =
        result.milp_has_objective_bound &&
        result.milp_best_objective_bound > threshold + tolerance;

    if (threshold_feasible) {
        result.outcome = FixedKSolveOutcome::Feasible;
        result.termination_name = "makespan-threshold-feasible";
        result.has_feasible_solution = true;
        extract_milp_incumbent(graph, core, result);
    } else if (threshold_infeasible_by_bound) {
        result.outcome = FixedKSolveOutcome::ProvenInfeasible;
        result.termination_name = "makespan-objective-bound-infeasible";
        result.infeasible = true;
    } else if (result.status == GRB_INFEASIBLE) {
        result.outcome = FixedKSolveOutcome::ProvenInfeasible;
        result.termination_name = "makespan-restricted-horizon-infeasible";
        result.infeasible = true;
    } else if (result.hit_time_limit) {
        result.outcome = FixedKSolveOutcome::TimedOutUnknown;
        result.termination_name = "timed-out-unknown";
    } else {
        result.outcome = FixedKSolveOutcome::EarlyUnknown;
        result.termination_name = "early-unknown";
        result.early_unknown = true;
    }
    return result;
}

}  // namespace

const char* initial_incumbent_milp_mode_name(InitialIncumbentMilpMode mode) {
    switch (mode) {
        case InitialIncumbentMilpMode::Duration:
            return "duration";
        case InitialIncumbentMilpMode::Makespan:
            return "makespan";
    }
    throw std::runtime_error("Unsupported initial-incumbent MILP mode.");
}

InitialIncumbentSolveResult solve_fixed_k_initial_incumbent(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const InitialIncumbentSolveOptions& options
) {
    validate_initial_incumbent_solve_options(options);
    switch (options.milp_mode) {
        case InitialIncumbentMilpMode::Duration:
            return solve_fixed_k_duration_impl(data, graph, options);
        case InitialIncumbentMilpMode::Makespan:
            return solve_fixed_k_makespan_impl(data, graph, options);
    }
    throw std::runtime_error("Unsupported initial-incumbent MILP mode.");
}

InitialIncumbentSolveResult solve_fixed_k_duration_initial_incumbent(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const InitialIncumbentSolveOptions& options
) {
    InitialIncumbentSolveOptions duration_options = options;
    duration_options.milp_mode = InitialIncumbentMilpMode::Duration;
    return solve_fixed_k_initial_incumbent(data, graph, duration_options);
}

InitialIncumbentSearchResult run_iterative_initial_incumbent_search(
    const InitialIncumbentSearchOptions& options,
    int max_candidate_k,
    const FixedKSolveCallback& solve_fixed_k
) {
    if (!solve_fixed_k) {
        throw std::runtime_error("Fixed-K incumbent callback is empty.");
    }
    InitialIncumbentSearchState state = make_initial_incumbent_search_state(
        options.initial_vehicle_count,
        options.max_k_increments,
        max_candidate_k,
        options.timeout_action
    );

    InitialIncumbentSearchResult search_result;
    search_result.initial_k = state.initial_k;

    while (!state.finished) {
        InitialIncumbentAttemptResult attempt;
        attempt.attempt_index = static_cast<int>(search_result.attempts.size());
        attempt.vehicle_count = state.candidate_k;
        attempt.solve_result = solve_fixed_k(state.candidate_k);
        search_result.total_runtime_seconds +=
            attempt.solve_result.runtime_seconds;
        search_result.last_candidate_k = state.candidate_k;

        if (attempt.solve_result.has_feasible_solution) {
            search_result.has_feasible_solution = true;
            search_result.incumbent_k = state.candidate_k;
            search_result.incumbent_result = attempt.solve_result;
        }

        search_result.attempts.push_back(std::move(attempt));
        state = advance_initial_incumbent_search(
            state,
            search_result.attempts.back().solve_result
        );
    }

    search_result.certified_k = state.certified_k;
    search_result.stopped_on_time_limit = state.stopped_on_time_limit;
    search_result.stopped_on_early_unknown = state.stopped_on_early_unknown;
    search_result.early_unknown_k = state.early_unknown_k;
    search_result.exhausted_increment_limit =
        state.exhausted_increment_limit;
    search_result.exhausted_candidate_limit =
        state.exhausted_candidate_limit;
    return search_result;
}

InitialIncumbentSearchResult solve_iterative_initial_incumbent(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const InitialIncumbentSearchOptions& options
) {
    if (!std::isfinite(options.per_attempt_time_limit) ||
        options.per_attempt_time_limit < 0.0) {
        throw std::runtime_error(
            "Initial-incumbent per-attempt time limit must be finite and nonnegative."
        );
    }

    auto log_path_for_k = [](const std::string& base_path_string,
                             const std::string& solver_suffix,
                             int vehicle_count) {
        if (base_path_string.empty()) {
            return std::string{};
        }
        const std::filesystem::path base_log_path(base_path_string);
        std::string stem = base_log_path.stem().string();
        const std::string suffix = "_" + solver_suffix;
        if (stem.size() >= suffix.size() &&
            stem.compare(stem.size() - suffix.size(), suffix.size(), suffix) == 0) {
            stem.erase(stem.size() - suffix.size());
            stem += "_k" + std::to_string(vehicle_count) + suffix;
        } else {
            stem += "_k" + std::to_string(vehicle_count);
        }
        return (base_log_path.parent_path() /
            (stem + base_log_path.extension().string())).string();
    };

    InitialIncumbentSearchResult result = run_iterative_initial_incumbent_search(
        options,
        std::max(options.initial_vehicle_count,
                 static_cast<int>(data.requests.size())),
        [&](int vehicle_count) {
            if (options.backend == InitialIncumbentBackend::CpSat) {
                CpSatSolveOptions cp_options = options.cp_sat;
                cp_options.vehicle_count = vehicle_count;
                cp_options.time_limit_seconds = options.per_attempt_time_limit;
                cp_options.log_file_path = log_path_for_k(
                    options.cp_sat_log_base_path, "cp_sat", vehicle_count
                );
                const CpFixedKSolveResult cp_result =
                    solve_fixed_k_cp_sat(data, cp_options);
                InitialIncumbentSolveResult converted;
                converted.backend = InitialIncumbentBackend::CpSat;
                converted.status = cp_result.raw_status;
                converted.raw_status = cp_result.raw_status;
                converted.status_name = cp_result.status_name;
                converted.termination_name = cp_result.termination_name;
                converted.runtime_seconds = cp_result.wall_time_seconds;
                converted.configured_time_limit_seconds =
                    cp_result.configured_time_limit_seconds;
                converted.cp_solve_mode = cp_result.solve_mode;
                converted.cp_threshold_horizon_factor =
                    cp_result.threshold_horizon_factor;
                converted.cp_model_horizon = cp_result.model_horizon;
                converted.cp_has_objective_value =
                    cp_result.has_objective_value;
                converted.cp_objective_value = cp_result.objective_value;
                converted.cp_has_objective_bound =
                    cp_result.has_objective_bound;
                converted.cp_best_objective_bound =
                    cp_result.best_objective_bound;
                converted.cp_stopped_by_feasible_observer =
                    cp_result.stopped_by_feasible_observer;
                converted.cp_stopped_by_bound_callback =
                    cp_result.stopped_by_bound_callback;
                converted.hit_time_limit = cp_result.hit_time_limit;
                converted.early_unknown = cp_result.early_unknown;
                converted.solution_info = cp_result.solution_info;
                converted.conflicts = cp_result.conflicts;
                converted.branches = cp_result.branches;
                converted.cp_build_stats = cp_result.build_stats;
                converted.vehicle_count = vehicle_count;
                converted.error_message = cp_result.error_message;
                switch (cp_result.outcome) {
                case CpSolveOutcome::ProvenInfeasible:
                    converted.outcome = FixedKSolveOutcome::ProvenInfeasible;
                    converted.infeasible = true;
                    return converted;
                case CpSolveOutcome::TimedOutUnknown:
                    converted.outcome = FixedKSolveOutcome::TimedOutUnknown;
                    return converted;
                case CpSolveOutcome::EarlyUnknown:
                    converted.outcome = FixedKSolveOutcome::EarlyUnknown;
                    return converted;
                case CpSolveOutcome::ModelInvalid:
                    converted.outcome = FixedKSolveOutcome::ModelInvalid;
                    return converted;
                case CpSolveOutcome::Feasible:
                    break;
                }

                const CpMappedIncumbent mapped = map_cp_incumbent_to_multigraph(
                    data, graph, cp_result.routes
                );
                if (!mapped.success) {
                    converted.outcome = FixedKSolveOutcome::AdapterError;
                    converted.termination_name = "adapter-error";
                    converted.adapter_status = "failed";
                    converted.error_message = mapped.error_message;
                    return converted;
                }
                converted.outcome = FixedKSolveOutcome::Feasible;
                converted.has_feasible_solution = true;
                converted.adapter_status = "passed";
                converted.active_edge_ids = mapped.active_edge_ids;
                converted.total_duration = mapped.total_duration;
                converted.total_original_cost = mapped.total_original_cost;
                converted.edge_values.assign(graph.number_of_edges(), 0.0);
                for (int edge_id : converted.active_edge_ids) {
                    converted.edge_values[static_cast<std::size_t>(edge_id)] = 1.0;
                }
                return converted;
            }
            InitialIncumbentSolveOptions attempt_options;
            attempt_options.vehicle_count = vehicle_count;
            attempt_options.solver_time_limit = options.per_attempt_time_limit;
            attempt_options.gurobi_threads = options.gurobi_threads;
            attempt_options.milp_mode = options.milp_mode;
            attempt_options.milp_makespan_horizon_factor =
                options.milp_makespan_horizon_factor;
            attempt_options.gurobi_log_path = log_path_for_k(
                options.gurobi_log_base_path, "gurobi", vehicle_count
            );
            return solve_fixed_k_initial_incumbent(
                data, graph, attempt_options
            );
        }
    );
    for (InitialIncumbentAttemptResult& attempt : result.attempts) {
        if (options.backend == InitialIncumbentBackend::TwoIndexMilp) {
            attempt.gurobi_log_path = log_path_for_k(
                options.gurobi_log_base_path,
                "gurobi",
                attempt.vehicle_count
            );
        }
    }
    return result;
}

InitialIncumbentSearchResult solve_iterative_duration_initial_incumbent(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const InitialIncumbentSearchOptions& options
) {
    InitialIncumbentSearchOptions duration_options = options;
    duration_options.milp_mode = InitialIncumbentMilpMode::Duration;
    return solve_iterative_initial_incumbent(data, graph, duration_options);
}

}  // namespace spdp
