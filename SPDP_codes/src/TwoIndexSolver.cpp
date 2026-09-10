#include "TwoIndexSolver.h"

#include <algorithm>
#include <chrono>
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

// User-cut callback: separates connectivity cuts on node LP relaxations.
class ConnectivityCutCallback : public GRBCallback {
public:
    ConnectivityCutCallback(
        const ConnectivitySeparator& separator,
        const ConnectivityCutOptions& options,
        std::vector<GRBVar>& y_vars,
        ConnectivityCutStats& stats
    )
        : separator_(separator),
          options_(options),
          y_vars_(y_vars),
          stats_(stats) {}

protected:
    void callback() override {
        if (where != GRB_CB_MIPNODE) {
            return;
        }
        if (getIntInfo(GRB_CB_MIPNODE_STATUS) != GRB_OPTIMAL) {
            return;
        }
        const double node_count = getDoubleInfo(GRB_CB_MIPNODE_NODCNT);
        const bool is_root = node_count < 0.5;
        if (options_.root_only && !is_root) {
            return;
        }
        const auto start = std::chrono::steady_clock::now();
        std::vector<double> values(y_vars_.size(), 0.0);
        double* raw = getNodeRel(y_vars_.data(), static_cast<int>(y_vars_.size()));
        for (std::size_t i = 0; i < values.size(); ++i) {
            values[i] = raw[i];
        }
        delete[] raw;
        const std::vector<PstepValidInequalityRow> cuts =
            separate_connectivity_cuts(separator_, values, options_);
        for (const PstepValidInequalityRow& row : cuts) {
            GRBLinExpr expression = 0.0;
            for (const auto& term : row.edge_terms) {
                expression += term.second * y_vars_[term.first];
            }
            addCut(expression >= row.rhs);
            ++stats_.cuts_added;
            if (is_root) {
                ++stats_.root_cuts_added;
            }
            if (row.rhs >= 2.0 - 1e-9) {
                ++stats_.rhs_two_or_more_cuts;
            }
        }
        ++stats_.rounds;
        stats_.separation_seconds +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }

private:
    const ConnectivitySeparator& separator_;
    const ConnectivityCutOptions& options_;
    std::vector<GRBVar>& y_vars_;
    ConnectivityCutStats& stats_;
};

}  // namespace

DirectTwoIndexResult solve_direct_two_index_model(
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
    core_options.binary_y = options.model_type == DirectTwoIndexModelType::IP;
    core_options.add_time_constraints = options.add_time_constraints;
    core_options.add_time_flow_formulation = options.add_time_flow_formulation;
    core_options.time_flow_state_disaggregated = options.time_flow_state_disaggregated;
    core_options.solver_time_limit = options.solver_time_limit;
    core_options.gurobi_threads = options.gurobi_threads;
    core_options.output_enabled = true;
    core_options.gurobi_log_path = options.gurobi_log_path;
    core_options.name_prefix = "direct_two_index";
    TwoIndexCoreModel core = build_two_index_model_core(data, graph, core_options);

    if (!options.initial_edge_start.empty()) {
        if (options.model_type != DirectTwoIndexModelType::IP) {
            throw std::runtime_error(
                "A direct two-index MIP start cannot be applied to an LP model."
            );
        }
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
        switch (row.sense) {
            case PstepValidInequalitySense::GreaterEqual:
                core.model->addConstr(expression >= row.rhs, row.name);
                break;
            case PstepValidInequalitySense::LessEqual:
                core.model->addConstr(expression <= row.rhs, row.name);
                break;
            case PstepValidInequalitySense::Equal:
                core.model->addConstr(expression == row.rhs, row.name);
                break;
        }
    }
    core.model->update();

    DirectTwoIndexResult result;
    result.model_type = options.model_type;
    result.variable_count = core.model->get(GRB_IntAttr_NumVars);
    result.time_flow_infeasible_edge_count = core.time_flow_infeasible_edge_count;
    result.constraint_count = core.model->get(GRB_IntAttr_NumConstrs);
    result.fixed_vehicle_constraint_count =
        options.valid_inequalities.add_fixed_vehicle_number ? 1 : 0;
    result.valid_inequality_count =
        static_cast<int>(vi_rows.size()) - result.fixed_vehicle_constraint_count;

    ConnectivitySeparator connectivity_separator;
    std::unique_ptr<ConnectivityCutCallback> connectivity_callback;
    if (options.connectivity_cuts.enabled &&
        options.model_type == DirectTwoIndexModelType::IP) {
        connectivity_separator = build_connectivity_separator(data, graph);
        connectivity_callback = std::make_unique<ConnectivityCutCallback>(
            connectivity_separator,
            options.connectivity_cuts,
            core.y_vars,
            result.connectivity_cut_stats
        );
        core.model->set(GRB_IntParam_PreCrush, 1);
        core.model->setCallback(connectivity_callback.get());
    }

    core.model->optimize();
    if (options.connectivity_cuts.enabled &&
        options.model_type == DirectTwoIndexModelType::LP) {
        // LP cutting-plane loop: separate connectivity cuts until none is violated.
        connectivity_separator = build_connectivity_separator(data, graph);
        std::vector<double> values(core.y_vars.size(), 0.0);
        for (std::size_t round = 0; round < 500; ++round) {
            if (core.model->get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
                break;
            }
            for (std::size_t i = 0; i < values.size(); ++i) {
                values[i] = core.y_vars[i].get(GRB_DoubleAttr_X);
            }
            const auto start = std::chrono::steady_clock::now();
            const std::vector<PstepValidInequalityRow> cuts =
                separate_connectivity_cuts(
                    connectivity_separator,
                    values,
                    options.connectivity_cuts
                );
            result.connectivity_cut_stats.separation_seconds +=
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - start
                ).count();
            ++result.connectivity_cut_stats.rounds;
            if (cuts.empty()) {
                break;
            }
            for (const PstepValidInequalityRow& row : cuts) {
                GRBLinExpr expression = 0.0;
                for (const auto& term : row.edge_terms) {
                    expression += term.second * core.y_vars[term.first];
                }
                core.model->addConstr(expression >= row.rhs, row.name);
                ++result.connectivity_cut_stats.cuts_added;
                ++result.connectivity_cut_stats.root_cuts_added;
                if (row.rhs >= 2.0 - 1e-9) {
                    ++result.connectivity_cut_stats.rhs_two_or_more_cuts;
                }
            }
            core.model->optimize();
        }
    }
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
    result.departure_flow = 0.0;
    result.edge_values.resize(graph.number_of_edges(), 0.0);
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const double value = core.y_vars[edge_id].get(GRB_DoubleAttr_X);
        result.edge_values[edge_id] = value;
        result.total_duration += graph.edges()[edge_id].data.time * value;
        result.total_original_cost += graph.edges()[edge_id].data.cost * value;
        const EdgeRecord& edge = graph.edges()[edge_id];
        if (edge.u == 0 &&
            graph.node(edge.v).kind == NodeSpec::Kind::Pickup) {
            result.departure_flow += value;
            if (options.model_type == DirectTwoIndexModelType::IP && value > 0.5) {
                ++result.vehicle_count;
            }
        }
    }
    result.total_duration_plus_fixed =
        result.total_duration +
        data.fixed_vehicle_cost * result.departure_flow;
    result.total_travel_cost =
        result.total_original_cost -
        data.fixed_vehicle_cost * result.departure_flow;

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
