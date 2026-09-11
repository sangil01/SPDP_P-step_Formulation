#include "TwoIndexSolver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
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

// Shared state of the three cut families for one direct two-index solve.
struct CutContext {
    // Type-closed travel-time cover rows (screened mode only).
    const std::vector<TypeTravelTimeCoverRow>* cover_rows = nullptr;
    std::vector<bool> cover_added;
    const TypeTravelTimeCoverOptions* cover_options = nullptr;
    TypeTravelTimeCoverStats* cover_stats = nullptr;
    // Capacity blossom cuts.
    const CapacityBlossomSeparator* blossom_separator = nullptr;
    const CapacityBlossomOptions* blossom_options = nullptr;
    CapacityBlossomPool blossom_pool;
    CapacityBlossomStats* blossom_stats = nullptr;
    std::size_t blossom_root_rounds = 0;
    std::size_t blossom_root_empty_rounds = 0;
    std::size_t blossom_root_weak_rounds = 0;
    bool blossom_root_closed = false;
    bool blossom_tree_closed = false;
    // Connectivity cuts.
    const ConnectivitySeparator* connectivity_separator = nullptr;
    const ConnectivityCutOptions* connectivity_options = nullptr;
    ConnectivityCutStats* connectivity_stats = nullptr;
};

bool cover_enabled(const CutContext& context) {
    return context.cover_rows != nullptr && context.cover_options != nullptr &&
        context.cover_options->enabled &&
        context.cover_options->mode == TypeTravelTimeCoverMode::Screened &&
        !context.cover_rows->empty();
}

bool blossom_enabled(const CutContext& context) {
    return context.blossom_separator != nullptr && context.blossom_options != nullptr &&
        context.blossom_options->enabled;
}

bool connectivity_enabled(const CutContext& context) {
    return context.connectivity_separator != nullptr &&
        context.connectivity_options != nullptr && context.connectivity_options->enabled;
}

// Decides whether the capacity blossom separator runs at this node and how many
// cuts it may return; zero means skip.
std::size_t blossom_budget_at_node(
    CutContext& context,
    bool is_root,
    double node_count,
    double runtime_seconds
) {
    const CapacityBlossomOptions& options = *context.blossom_options;
    CapacityBlossomStats& stats = *context.blossom_stats;
    if (stats.cuts_added >= options.max_total_cuts) {
        stats.stopped_by_total_cap = true;
        return 0;
    }
    if (is_root) {
        if (context.blossom_root_closed) {
            return 0;
        }
        if (context.blossom_root_rounds >= options.root_max_rounds ||
            stats.root_cuts_added >= options.root_max_cuts) {
            context.blossom_root_closed = true;
            return 0;
        }
        return std::min(
            options.root_max_per_round,
            options.root_max_cuts - stats.root_cuts_added
        );
    }
    if (options.scope == CapacityBlossomScope::RootOnly || context.blossom_tree_closed) {
        return 0;
    }
    if (options.max_separation_time_fraction > 0.0 && runtime_seconds > 5.0 &&
        stats.separation_seconds > options.max_separation_time_fraction * runtime_seconds) {
        stats.stopped_by_time_fraction = true;
        context.blossom_tree_closed = true;
        return 0;
    }
    if (options.scope == CapacityBlossomScope::AdaptiveTree) {
        const auto node = static_cast<unsigned long long>(node_count + 0.5);
        const bool dense_phase = node < options.tree_dense_node_limit;
        const bool periodic = options.tree_node_frequency > 0U &&
            node % options.tree_node_frequency == 0ULL;
        if (!dense_phase && !periodic) {
            return 0;
        }
    }
    return options.tree_max_per_round;
}

// Applies one round of every enabled family to the fractional point y and
// returns the rows to add. Also updates statistics and pools.
struct RoundOutput {
    std::vector<PstepValidInequalityRow> rows;
    std::size_t cover_count = 0;
    std::size_t blossom_count = 0;
    std::size_t connectivity_count = 0;
    double blossom_best_violation = 0.0;
};

RoundOutput separate_round(
    CutContext& context,
    const std::vector<double>& y_values,
    bool is_root,
    double node_count,
    double runtime_seconds,
    bool lp_mode
) {
    RoundOutput output;

    if (cover_enabled(context)) {
        const auto start = std::chrono::steady_clock::now();
        const std::vector<std::size_t> indices = screen_type_travel_time_cover_rows(
            *context.cover_rows,
            y_values,
            context.cover_added,
            context.cover_options->violation_tolerance,
            context.cover_options->max_per_round
        );
        for (std::size_t index : indices) {
            context.cover_added[index] = true;
            output.rows.push_back((*context.cover_rows)[index].row);
            ++output.cover_count;
            ++context.cover_stats->rows_added;
            if (is_root) {
                ++context.cover_stats->root_rows_added;
            }
        }
        ++context.cover_stats->screening_rounds;
        context.cover_stats->screening_seconds +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }

    if (blossom_enabled(context)) {
        std::size_t budget = 0;
        if (lp_mode) {
            budget = context.blossom_options->root_max_per_round;
        } else {
            budget = blossom_budget_at_node(context, is_root, node_count, runtime_seconds);
        }
        if (budget > 0U) {
            const std::vector<CapacityBlossomCut> cuts = separate_capacity_blossom_cuts(
                *context.blossom_separator,
                y_values,
                *context.blossom_options,
                budget,
                context.blossom_pool,
                *context.blossom_stats
            );
            for (const CapacityBlossomCut& cut : cuts) {
                output.rows.push_back(cut.row);
                output.blossom_best_violation =
                    std::max(output.blossom_best_violation, cut.violation);
            }
            output.blossom_count = cuts.size();
            ++context.blossom_stats->rounds;
            context.blossom_stats->cuts_added += cuts.size();
            if (is_root) {
                ++context.blossom_stats->root_rounds;
                context.blossom_stats->root_cuts_added += cuts.size();
                if (!lp_mode) {
                    ++context.blossom_root_rounds;
                    if (cuts.empty()) {
                        ++context.blossom_root_empty_rounds;
                    } else {
                        context.blossom_root_empty_rounds = 0;
                    }
                    if (output.blossom_best_violation <
                        context.blossom_options->root_tailing_off_violation) {
                        ++context.blossom_root_weak_rounds;
                    } else {
                        context.blossom_root_weak_rounds = 0;
                    }
                    if (context.blossom_root_empty_rounds >= 2U ||
                        context.blossom_root_weak_rounds >= 2U) {
                        context.blossom_root_closed = true;
                    }
                }
            }
        }
    }

    if (connectivity_enabled(context) &&
        (is_root || !context.connectivity_options->root_only)) {
        const auto start = std::chrono::steady_clock::now();
        const std::vector<PstepValidInequalityRow> cuts = separate_connectivity_cuts(
            *context.connectivity_separator,
            y_values,
            *context.connectivity_options
        );
        for (const PstepValidInequalityRow& row : cuts) {
            output.rows.push_back(row);
            ++context.connectivity_stats->cuts_added;
            if (is_root) {
                ++context.connectivity_stats->root_cuts_added;
            }
            if (row.rhs >= 2.0 - 1e-9) {
                ++context.connectivity_stats->rhs_two_or_more_cuts;
            }
        }
        output.connectivity_count = cuts.size();
        ++context.connectivity_stats->rounds;
        context.connectivity_stats->separation_seconds +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }
    return output;
}

GRBTempConstr row_to_temp_constr(
    const PstepValidInequalityRow& row,
    const std::vector<GRBVar>& y_vars
) {
    GRBLinExpr expression = 0.0;
    for (const auto& term : row.edge_terms) {
        expression += term.second * y_vars[term.first];
    }
    switch (row.sense) {
        case PstepValidInequalitySense::GreaterEqual:
            return expression >= row.rhs;
        case PstepValidInequalitySense::LessEqual:
            return expression <= row.rhs;
        case PstepValidInequalitySense::Equal:
            return expression == row.rhs;
    }
    throw std::runtime_error("Unknown valid inequality sense.");
}

// User-cut callback shared by every cut family. One node relaxation is read per
// invocation and handed to all separators.
class DirectTwoIndexCutCallback : public GRBCallback {
public:
    DirectTwoIndexCutCallback(CutContext& context, std::vector<GRBVar>& y_vars)
        : context_(context), y_vars_(y_vars) {}

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
        const double runtime_seconds = getDoubleInfo(GRB_CB_RUNTIME);

        const bool want_cover = cover_enabled(context_) &&
            std::find(context_.cover_added.begin(), context_.cover_added.end(), false) !=
                context_.cover_added.end();
        const bool want_blossom = blossom_enabled(context_) &&
            (is_root ? !context_.blossom_root_closed
                     : (context_.blossom_options->scope != CapacityBlossomScope::RootOnly &&
                        !context_.blossom_tree_closed));
        const bool want_connectivity = connectivity_enabled(context_) &&
            (is_root || !context_.connectivity_options->root_only);
        if (!want_cover && !want_blossom && !want_connectivity) {
            return;
        }

        std::vector<double> values(y_vars_.size(), 0.0);
        double* raw = getNodeRel(y_vars_.data(), static_cast<int>(y_vars_.size()));
        for (std::size_t i = 0; i < values.size(); ++i) {
            values[i] = raw[i];
        }
        delete[] raw;

        const RoundOutput output =
            separate_round(context_, values, is_root, node_count, runtime_seconds, false);
        for (const PstepValidInequalityRow& row : output.rows) {
            addCut(row_to_temp_constr(row, y_vars_));
        }
    }

private:
    CutContext& context_;
    std::vector<GRBVar>& y_vars_;
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
        for (const auto& term : row.edge_terms) {
            if (term.first >= core.y_vars.size()) {
                throw std::runtime_error(
                    "A direct two-index valid inequality references an invalid edge."
                );
            }
        }
        core.model->addConstr(row_to_temp_constr(row, core.y_vars), row.name);
    }

    DirectTwoIndexResult result;
    result.type_cover_stats.rows_available = options.type_cover_rows.size();
    if (options.type_cover.enabled &&
        options.type_cover.mode == TypeTravelTimeCoverMode::Static) {
        for (const TypeTravelTimeCoverRow& row : options.type_cover_rows) {
            core.model->addConstr(row_to_temp_constr(row.row, core.y_vars), row.row.name);
            ++result.type_cover_stats.rows_added;
            ++result.type_cover_stats.root_rows_added;
        }
    }
    core.model->update();

    result.model_type = options.model_type;
    result.variable_count = core.model->get(GRB_IntAttr_NumVars);
    result.time_flow_infeasible_edge_count = core.time_flow_infeasible_edge_count;
    result.constraint_count = core.model->get(GRB_IntAttr_NumConstrs);
    result.fixed_vehicle_constraint_count =
        options.valid_inequalities.add_fixed_vehicle_number ? 1 : 0;
    result.valid_inequality_count =
        static_cast<int>(vi_rows.size()) - result.fixed_vehicle_constraint_count;

    CutContext context;
    ConnectivitySeparator connectivity_separator;
    CapacityBlossomSeparator blossom_separator;
    if (options.type_cover.enabled &&
        options.type_cover.mode == TypeTravelTimeCoverMode::Screened) {
        context.cover_rows = &options.type_cover_rows;
        context.cover_added.assign(options.type_cover_rows.size(), false);
        context.cover_options = &options.type_cover;
        context.cover_stats = &result.type_cover_stats;
    }
    if (options.capacity_blossom.enabled) {
        blossom_separator = build_capacity_blossom_separator(graph);
        context.blossom_separator = &blossom_separator;
        context.blossom_options = &options.capacity_blossom;
        context.blossom_stats = &result.capacity_blossom_stats;
    }
    if (options.connectivity_cuts.enabled) {
        connectivity_separator = build_connectivity_separator(data, graph);
        context.connectivity_separator = &connectivity_separator;
        context.connectivity_options = &options.connectivity_cuts;
        context.connectivity_stats = &result.connectivity_cut_stats;
    }
    const bool any_dynamic_family =
        cover_enabled(context) || blossom_enabled(context) || connectivity_enabled(context);

    std::unique_ptr<DirectTwoIndexCutCallback> cut_callback;
    if (any_dynamic_family && options.model_type == DirectTwoIndexModelType::IP) {
        cut_callback = std::make_unique<DirectTwoIndexCutCallback>(context, core.y_vars);
        core.model->set(GRB_IntParam_PreCrush, 1);
        core.model->setCallback(cut_callback.get());
    }

    core.model->optimize();
    if (any_dynamic_family && options.model_type == DirectTwoIndexModelType::LP) {
        // LP cutting-plane loop: every family is separated until no family
        // returns a violated row (or the round limit is reached).
        std::vector<double> values(core.y_vars.size(), 0.0);
        for (std::size_t round = 0; round < 500; ++round) {
            if (core.model->get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
                break;
            }
            for (std::size_t i = 0; i < values.size(); ++i) {
                values[i] = core.y_vars[i].get(GRB_DoubleAttr_X);
            }
            const RoundOutput output = separate_round(context, values, true, 0.0, 0.0, true);
            if (output.rows.empty()) {
                break;
            }
            for (const PstepValidInequalityRow& row : output.rows) {
                core.model->addConstr(row_to_temp_constr(row, core.y_vars), row.name);
            }
            ++result.lp_cut_rounds;
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
