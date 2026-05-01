#include "PstepBnP.h"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace spdp {
namespace {

constexpr double kTolerance = 1e-9;

bool double_less_or_equal(double lhs, double rhs) {
    return lhs <= rhs + kTolerance;
}

std::string format_double(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(6) << value;
    return out.str();
}

bool is_lp_status_usable(int status) {
    return status == GRB_OPTIMAL || status == GRB_SUBOPTIMAL;
}

double remaining_time_seconds(
    const std::chrono::steady_clock::time_point start_time,
    double total_time_limit_seconds
) {
    if (total_time_limit_seconds < 0.0) {
        return total_time_limit_seconds;
    }
    const double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();
    return std::max(0.0, total_time_limit_seconds - elapsed);
}

void write_iteration_log(
    std::ostream& out,
    const NodeCGIterationLog& log
) {
    out << "[node-cg] phase="
        << (log.phase == CGPhase::PhaseI ? "I" : "II")
        << " iter=" << log.iteration_index
        << " lp_obj=" << format_double(log.lp_objective_value)
        << " artificial_sum=" << format_double(log.artificial_sum)
        << " best_rc=" << format_double(log.best_reduced_cost)
        << " pricing_status=" << to_string(log.pricing_status)
        << " pricing_seconds=" << format_double(log.pricing_runtime_seconds)
        << " added_columns=" << log.added_column_count
        << " complete_labels=" << log.complete_label_count
        << " start_labels=" << log.start_label_count
        << " generated_labels=" << log.generated_label_count
        << " surviving_labels=" << log.surviving_label_count
        << " dominated_labels=" << log.dominated_label_count
        << '\n';
}

void capture_column_values_if_available(
    const CGMasterProblem& master_problem,
    NodeCGResult& result
) {
    result.columns = master_problem.columns;
    result.column_values.clear();

    if (master_problem.model->get(GRB_IntAttr_SolCount) <= 0) {
        return;
    }

    result.column_values.resize(master_problem.x_vars.size(), 0.0);
    for (std::size_t idx = 0; idx < master_problem.x_vars.size(); ++idx) {
        result.column_values[idx] = master_problem.x_vars[idx].get(GRB_DoubleAttr_X);
    }
}

bool run_node_pricing_phase(
    const MultiDiGraph& graph,
    const ForwardPricingContext& pricing_context,
    CGMasterProblem& master_problem,
    const NodeCGOptions& options,
    CGPhase phase,
    const std::chrono::steady_clock::time_point global_start_time,
    NodeCGResult& result,
    std::ostream* log_stream
) {
    const ForwardPricingOptions pricing_options{
        options.max_columns_per_start,
        options.max_total_columns_per_round,
        options.reduced_cost_tolerance,
        options.prune_pickup_symmetry_43,
        options.prune_delivery_symmetry_43,
    };
    CGLPSnapshot last_solved_snapshot;
    bool has_last_solved_snapshot = false;

    for (std::size_t iteration = 1; iteration <= options.max_iterations_per_phase; ++iteration) {
        const double remaining = remaining_time_seconds(global_start_time, options.solver_time_limit);
        if (remaining >= 0.0) {
            master_problem.model->set(GRB_DoubleParam_TimeLimit, remaining);
        }

        solve_cg_master_lp(master_problem);
        const CGLPSnapshot snapshot = capture_cg_master_snapshot(master_problem);
        last_solved_snapshot = snapshot;
        has_last_solved_snapshot = true;

        NodeCGIterationLog iteration_log;
        iteration_log.phase = phase;
        iteration_log.iteration_index = iteration;
        iteration_log.lp_objective_value = snapshot.objective_value;
        iteration_log.artificial_sum = snapshot.artificial_sum;

        if (log_stream != nullptr) {
            write_cg_master_snapshot(*log_stream, snapshot);
        }

        if (!is_lp_status_usable(snapshot.gurobi_status)) {
            result.final_snapshot = snapshot;
            result.iteration_logs.push_back(iteration_log);
            if (phase == CGPhase::PhaseI) {
                result.phase_one_objective = snapshot.objective_value;
                result.phase_one_iterations = iteration;
            } else {
                result.phase_two_objective = snapshot.objective_value;
                result.phase_two_iterations = iteration;
            }
            if (log_stream != nullptr) {
                *log_stream << "[node-cg] stopping because LP status is not usable for pricing.\n";
            }
            return false;
        }

        if (phase == CGPhase::PhaseI && snapshot.artificial_sum <= 1e-6) {
            result.phase_one_feasible = true;
            result.phase_one_objective = snapshot.objective_value;
            result.phase_one_iterations = iteration;
            result.final_snapshot = snapshot;
            result.iteration_logs.push_back(iteration_log);
            if (log_stream != nullptr) {
                *log_stream << "[node-cg] Phase I reached zero artificial sum.\n";
            }
            return true;
        }

        const CGDualSolution dual_solution = extract_cg_master_duals(master_problem);
        const auto pricing_start_time = std::chrono::steady_clock::now();
        const ForwardPricingResult pricing_result = run_forward_pricing(
            graph,
            pricing_context,
            dual_solution,
            phase,
            pricing_options
        );
        const auto pricing_end_time = std::chrono::steady_clock::now();
        iteration_log.best_reduced_cost = pricing_result.best_reduced_cost;
        iteration_log.pricing_status = pricing_result.status;
        iteration_log.pricing_runtime_seconds =
            std::chrono::duration<double>(pricing_end_time - pricing_start_time).count();
        iteration_log.complete_label_count = pricing_result.complete_label_count;
        iteration_log.start_label_count = pricing_result.start_label_count;
        iteration_log.generated_label_count = pricing_result.generated_label_count;
        iteration_log.surviving_label_count = pricing_result.surviving_label_count;
        iteration_log.dominated_label_count = pricing_result.dominated_label_count;

        const std::size_t added_columns =
            add_columns_to_cg_master(master_problem, pricing_result.columns);
        iteration_log.added_column_count = added_columns;
        result.iteration_logs.push_back(iteration_log);
        result.total_columns_added += added_columns;

        if (log_stream != nullptr) {
            write_iteration_log(*log_stream, iteration_log);
        }

        if (added_columns == 0U) {
            result.final_snapshot = snapshot;
            if (phase == CGPhase::PhaseI) {
                result.phase_one_feasible = false;
                result.phase_one_objective = snapshot.objective_value;
                result.phase_one_iterations = iteration;
            } else {
                result.phase_two_objective = snapshot.objective_value;
                result.phase_two_iterations = iteration;
                result.solved_to_completion = true;
            }
            return phase == CGPhase::PhaseII;
        }
    }

    result.final_snapshot = has_last_solved_snapshot
                                ? last_solved_snapshot
                                : capture_cg_master_snapshot(master_problem);
    if (phase == CGPhase::PhaseI) {
        result.phase_one_objective = result.final_snapshot.objective_value;
        result.phase_one_iterations = options.max_iterations_per_phase;
    } else {
        result.phase_two_objective = result.final_snapshot.objective_value;
        result.phase_two_iterations = options.max_iterations_per_phase;
    }
    return false;
}

bool run_node_phase_one_exact_cg(
    const MultiDiGraph& graph,
    const ForwardPricingContext& pricing_context,
    CGMasterProblem& master_problem,
    const NodeCGOptions& options,
    const std::chrono::steady_clock::time_point global_start_time,
    NodeCGResult& result,
    std::ostream* log_stream
) {
    if (log_stream != nullptr) {
        *log_stream << "[node-cg] Phase I mode=exact-cg\n";
    }
    return run_node_pricing_phase(
        graph,
        pricing_context,
        master_problem,
        options,
        CGPhase::PhaseI,
        global_start_time,
        result,
        log_stream
    );
}

bool run_node_phase_one_heuristic_seed(
    NodeCGResult& result,
    std::ostream* log_stream
) {
    result.phase_one_feasible = false;
    if (log_stream != nullptr) {
        *log_stream << "[node-cg] Phase I mode=heuristic-seed\n";
    }
    throw std::runtime_error(
        "node CG phase-one mode heuristic-seed is not implemented yet."
    );
}

bool run_node_phase_one(
    const MultiDiGraph& graph,
    const ForwardPricingContext& pricing_context,
    CGMasterProblem& master_problem,
    const NodeCGOptions& options,
    const std::chrono::steady_clock::time_point global_start_time,
    NodeCGResult& result,
    std::ostream* log_stream
) {
    switch (options.phase_one_mode) {
        case NodeCGPhaseOneMode::ExactCG:
            return run_node_phase_one_exact_cg(
                graph,
                pricing_context,
                master_problem,
                options,
                global_start_time,
                result,
                log_stream
            );
        case NodeCGPhaseOneMode::HeuristicSeed:
            return run_node_phase_one_heuristic_seed(result, log_stream);
    }
    throw std::runtime_error("Unsupported node CG phase-one mode.");
}

bool run_node_phase_two(
    const MultiDiGraph& graph,
    const ForwardPricingContext& pricing_context,
    CGMasterProblem& master_problem,
    const NodeCGOptions& options,
    const std::chrono::steady_clock::time_point global_start_time,
    NodeCGResult& result,
    std::ostream* log_stream
) {
    switch_cg_master_to_phase_two(master_problem);
    result.phase_two_reached = true;
    if (log_stream != nullptr) {
        *log_stream << "[node-cg] Entering Phase II.\n";
    }
    return run_node_pricing_phase(
        graph,
        pricing_context,
        master_problem,
        options,
        CGPhase::PhaseII,
        global_start_time,
        result,
        log_stream
    );
}

}  // namespace

NodeCGResult solve_node_column_generation(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const NodeCGOptions& options,
    std::ostream* log_stream
) {
    if (options.p < 1) {
        throw std::runtime_error("Node CG requires p >= 1.");
    }
    if (!double_less_or_equal(0.0, options.time_limit)) {
        throw std::runtime_error("Node CG requires a nonnegative time limit.");
    }

    NodeCGResult result;
    result.phase_one_mode = options.phase_one_mode;

    const ForwardPricingContext pricing_context =
        build_forward_pricing_context(graph, options.p, options.time_limit);

    const CGMasterOptions master_options{
        options.gurobi_log_path,
        options.gurobi_threads,
        options.solver_time_limit,
    };
    CGMasterProblem master_problem =
        build_cg_master_problem(data, graph, master_options);

    const auto global_start_time = std::chrono::steady_clock::now();

    const bool phase_one_success = run_node_phase_one(
        graph,
        pricing_context,
        master_problem,
        options,
        global_start_time,
        result,
        log_stream
    );
    if (!phase_one_success || !result.phase_one_feasible) {
        capture_column_values_if_available(master_problem, result);
        return result;
    }

    const bool phase_two_success = run_node_phase_two(
        graph,
        pricing_context,
        master_problem,
        options,
        global_start_time,
        result,
        log_stream
    );
    capture_column_values_if_available(master_problem, result);
    if (!phase_two_success) {
        return result;
    }

    if (master_problem.model->get(GRB_IntAttr_SolCount) > 0) {
        result.phase_two_objective = master_problem.model->get(GRB_DoubleAttr_ObjVal);
    }
    return result;
}

void write_node_cg_summary(
    std::ostream& out,
    const NodeCGResult& result
) {
    out << "[node-cg-summary] phase_one_mode=" << to_string(result.phase_one_mode) << '\n';
    out << "[node-cg-summary] phase_one_feasible=" << (result.phase_one_feasible ? 1 : 0) << '\n';
    out << "[node-cg-summary] phase_two_reached=" << (result.phase_two_reached ? 1 : 0) << '\n';
    out << "[node-cg-summary] solved_to_completion=" << (result.solved_to_completion ? 1 : 0) << '\n';
    out << "[node-cg-summary] phase_one_objective=" << format_double(result.phase_one_objective) << '\n';
    out << "[node-cg-summary] phase_two_objective=" << format_double(result.phase_two_objective) << '\n';
    out << "[node-cg-summary] phase_one_iterations=" << result.phase_one_iterations << '\n';
    out << "[node-cg-summary] phase_two_iterations=" << result.phase_two_iterations << '\n';
    out << "[node-cg-summary] total_columns_added=" << result.total_columns_added << '\n';
    write_cg_master_snapshot(out, result.final_snapshot);
}

void write_node_cg_solution(
    std::ostream& out,
    const NodeCGResult& result
) {
    out << "Solution:\n";

    if (!result.phase_one_feasible) {
        out << "  Node CG stopped before obtaining a Phase-I feasible RMP\n";
        out << "Solution done \n";
        return;
    }

    if (!result.solved_to_completion) {
        out << "  Node CG terminated before proving LP completion\n";
        out << "Solution done \n";
        return;
    }

    for (std::size_t idx = 0; idx < result.columns.size() && idx < result.column_values.size(); ++idx) {
        if (result.column_values[idx] <= 1e-6) {
            continue;
        }

        const CGColumn& column = result.columns[idx];
        out << "  x_" << column.id
            << " = " << format_double(result.column_values[idx])
            << " | rc(gen)=" << format_double(column.reduced_cost)
            << " | tau=" << format_double(column.tau)
            << " | time=" << format_double(column.total_time)
            << " | cost=" << format_double(column.total_cost)
            << " | edges=";

        out << "[";
        for (std::size_t edge_idx = 0; edge_idx < column.edge_ids.size(); ++edge_idx) {
            if (edge_idx > 0U) {
                out << ", ";
            }
            out << column.edge_ids[edge_idx];
        }
        out << "]\n";
    }

    out << "Solution done \n";
}

const char* to_string(NodeCGPhaseOneMode mode) {
    switch (mode) {
        case NodeCGPhaseOneMode::ExactCG:
            return "exact-cg";
        case NodeCGPhaseOneMode::HeuristicSeed:
            return "heuristic-seed";
    }
    return "unknown";
}

}  // namespace spdp
