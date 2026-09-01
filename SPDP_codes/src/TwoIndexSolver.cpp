#include "TwoIndexSolver.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "TwoIndexModelCore.h"
#include "gurobi_c++.h"

namespace spdp {
namespace {

using detail::TwoIndexCoreModel;
using detail::TwoIndexCoreObjective;
using detail::TwoIndexCoreOptions;
using detail::build_two_index_model_core;

}  // namespace

DirectTwoIndexResult solve_direct_two_index_ip(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const DirectTwoIndexOptions& options
) {
    TwoIndexCoreOptions core_options;
    switch (options.objective) {
        case DirectTwoIndexObjective::OriginalCost:
            core_options.objective = TwoIndexCoreObjective::OriginalCost;
            break;
        case DirectTwoIndexObjective::TravelCost:
            core_options.objective = TwoIndexCoreObjective::TravelCost;
            break;
        case DirectTwoIndexObjective::Duration:
            core_options.objective = TwoIndexCoreObjective::Duration;
            break;
        case DirectTwoIndexObjective::DurationPlusFixed:
            core_options.objective = TwoIndexCoreObjective::DurationPlusFixed;
            break;
    }
    core_options.binary_y = true;
    core_options.add_time_constraints = options.add_time_constraints;
    core_options.solver_time_limit = options.solver_time_limit;
    core_options.gurobi_threads = options.gurobi_threads;
    core_options.output_enabled = true;
    core_options.gurobi_log_path = options.gurobi_log_path;
    core_options.name_prefix = "direct_two_index";
    TwoIndexCoreModel core = build_two_index_model_core(data, graph, core_options);

    if (!options.initial_edge_start.empty()) {
        if (options.initial_edge_start.size() != core.y_vars.size()) {
            throw std::runtime_error(
                "The direct two-index MIP start size does not match the edge count."
            );
        }
        for (std::size_t edge_id = 0; edge_id < core.y_vars.size(); ++edge_id) {
            core.y_vars[edge_id].set(
                GRB_DoubleAttr_Start,
                options.initial_edge_start[edge_id]
            );
        }
    }

    const std::vector<PstepValidInequalityRow> vi_rows =
        build_pstep_valid_inequality_rows(
            data,
            graph,
            options.valid_inequalities
        );
    for (const PstepValidInequalityRow& row : vi_rows) {
        GRBLinExpr expression = 0.0;
        for (const auto& term : row.edge_terms) {
            if (term.first >= core.y_vars.size()) {
                throw std::runtime_error(
                    "A direct two-index valid inequality references an invalid edge."
                );
            }
            expression += term.second * core.y_vars[term.first];
        }
        if (row.sense == PstepValidInequalitySense::GreaterEqual) {
            core.model->addConstr(expression >= row.rhs, row.name);
        } else {
            core.model->addConstr(expression <= row.rhs, row.name);
        }
    }
    core.model->update();

    DirectTwoIndexResult result;
    result.variable_count = core.model->get(GRB_IntAttr_NumVars);
    result.constraint_count = core.model->get(GRB_IntAttr_NumConstrs);
    result.valid_inequality_count = static_cast<int>(vi_rows.size());

    core.model->optimize();
    result.status = core.model->get(GRB_IntAttr_Status);
    result.hit_time_limit = result.status == GRB_TIME_LIMIT;
    result.solved_to_optimality = result.status == GRB_OPTIMAL;
    result.runtime_seconds = core.model->get(GRB_DoubleAttr_Runtime);
    result.has_feasible_solution = core.model->get(GRB_IntAttr_SolCount) > 0;

    try {
        const double objective_bound = core.model->get(GRB_DoubleAttr_ObjBound);
        if (std::isfinite(objective_bound) &&
            std::abs(objective_bound) < 0.5 * GRB_INFINITY) {
            result.objective_bound = objective_bound;
            result.has_certified_bound = true;
        }
    } catch (const GRBException& error) {
        if (error.getErrorCode() != GRB_ERROR_DATA_NOT_AVAILABLE) {
            throw;
        }
    }

    if (!result.has_feasible_solution) {
        return result;
    }

    result.objective_value = core.model->get(GRB_DoubleAttr_ObjVal);
    if (result.has_certified_bound &&
        std::abs(result.objective_value) > detail::kGurobiSolverTolerance) {
        result.gap_percent =
            100.0 * std::abs(result.objective_value - result.objective_bound) /
            std::abs(result.objective_value);
    } else if (result.has_certified_bound) {
        result.gap_percent = 0.0;
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
    }
    result.total_duration_plus_fixed =
        result.total_duration +
        data.fixed_vehicle_cost * static_cast<double>(result.vehicle_count);
    result.total_travel_cost =
        result.total_original_cost -
        data.fixed_vehicle_cost * static_cast<double>(result.vehicle_count);

    double recomputed_objective = 0.0;
    switch (options.objective) {
        case DirectTwoIndexObjective::OriginalCost:
            recomputed_objective = result.total_original_cost;
            break;
        case DirectTwoIndexObjective::TravelCost:
            recomputed_objective = result.total_travel_cost;
            break;
        case DirectTwoIndexObjective::Duration:
            recomputed_objective = result.total_duration;
            break;
        case DirectTwoIndexObjective::DurationPlusFixed:
            recomputed_objective = result.total_duration_plus_fixed;
            break;
    }
    const double objective_scale =
        std::max({1.0, std::abs(result.objective_value), std::abs(recomputed_objective)});
    if (std::abs(result.objective_value - recomputed_objective) >
        detail::kGurobiSolverTolerance * objective_scale) {
        throw std::runtime_error(
            "The direct two-index objective value is inconsistent with its components."
        );
    }
    return result;
}

}  // namespace spdp
