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
    FixedKSolveOutcome outcome
) {
    if (state.finished) {
        throw std::runtime_error(
            "Cannot advance a finished initial incumbent search."
        );
    }

    InitialIncumbentSearchState next = state;
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

    if (outcome == FixedKSolveOutcome::Unknown &&
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

}  // namespace

InitialIncumbentSolveResult solve_fixed_k_duration_initial_incumbent(
    const SPDPData& data,
    const MultiDiGraph& graph,
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

    GRBLinExpr actual_departures = 0.0;
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        if (edge.u == 0 &&
            graph.node(edge.v).kind == NodeSpec::Kind::Pickup) {
            actual_departures += core.y_vars[edge_id];
        }
    }
    core.model->addConstr(
        actual_departures == static_cast<double>(options.vehicle_count),
        "initial_incumbent_fixed_vehicle_count"
    );
    core.model->set(GRB_IntParam_SolutionLimit, 1);
    core.model->update();
    core.model->optimize();

    InitialIncumbentSolveResult result;
    result.status = core.model->get(GRB_IntAttr_Status);
    result.raw_status = result.status;
    result.status_name = std::to_string(result.status);
    result.hit_time_limit = result.status == GRB_TIME_LIMIT;
    result.hit_solution_limit = result.status == GRB_SOLUTION_LIMIT;
    result.infeasible = result.status == GRB_INFEASIBLE;
    result.runtime_seconds = core.model->get(GRB_DoubleAttr_Runtime);
    result.has_feasible_solution = core.model->get(GRB_IntAttr_SolCount) > 0;
    result.outcome = result.has_feasible_solution
        ? FixedKSolveOutcome::Feasible
        : (result.infeasible ? FixedKSolveOutcome::ProvenInfeasible
                             : FixedKSolveOutcome::Unknown);
    if (!result.has_feasible_solution) {
        return result;
    }

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
    return result;
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

        const FixedKSolveOutcome outcome = attempt.solve_result.outcome;
        search_result.attempts.push_back(std::move(attempt));
        state = advance_initial_incumbent_search(state, outcome);
    }

    search_result.certified_k = state.certified_k;
    search_result.stopped_on_time_limit = state.stopped_on_time_limit;
    search_result.exhausted_increment_limit =
        state.exhausted_increment_limit;
    search_result.exhausted_candidate_limit =
        state.exhausted_candidate_limit;
    return search_result;
}

InitialIncumbentSearchResult solve_iterative_duration_initial_incumbent(
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

    const std::filesystem::path base_log_path(options.gurobi_log_base_path);
    auto log_path_for_k = [&base_log_path, &options](int vehicle_count) {
        if (options.gurobi_log_base_path.empty()) {
            return std::string{};
        }
        std::string stem = base_log_path.stem().string();
        const std::string suffix = "_gurobi";
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
                const CpFixedKSolveResult cp_result =
                    solve_fixed_k_cp_sat(data, cp_options);
                InitialIncumbentSolveResult converted;
                converted.backend = InitialIncumbentBackend::CpSat;
                converted.status = cp_result.raw_status;
                converted.raw_status = cp_result.raw_status;
                converted.status_name = cp_result.status_name;
                converted.runtime_seconds = cp_result.wall_time_seconds;
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
                case CpSolveOutcome::Unknown:
                    converted.outcome = FixedKSolveOutcome::Unknown;
                    converted.hit_time_limit = options.per_attempt_time_limit > 0.0;
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
            attempt_options.gurobi_log_path = log_path_for_k(vehicle_count);
            return solve_fixed_k_duration_initial_incumbent(
                data, graph, attempt_options
            );
        }
    );
    for (InitialIncumbentAttemptResult& attempt : result.attempts) {
        if (options.backend == InitialIncumbentBackend::TwoIndexMilp) {
            attempt.gurobi_log_path = log_path_for_k(attempt.vehicle_count);
        }
    }
    return result;
}

}  // namespace spdp
