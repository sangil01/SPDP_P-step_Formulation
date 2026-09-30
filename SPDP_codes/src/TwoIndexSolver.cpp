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

// Shared state of the dynamic cut families for one direct two-index solve.
// Both families run the same management policy (scope, caps, violation
// thresholds, tailing-off, time budget) through independent option and state
// instances; only candidate generation differs.
struct CutContext {
    // Type-closed travel-time cover rows (screened mode only).
    const std::vector<TypeTravelTimeCoverRow>* cover_rows = nullptr;
    std::vector<bool> cover_added;
    const TypeTravelTimeCoverOptions* cover_options = nullptr;
    TypeTravelTimeCoverStats* cover_stats = nullptr;
    CutManagementState cover_state;
    // Capacity blossom cuts.
    const CapacityBlossomSeparator* blossom_separator = nullptr;
    const CapacityBlossomOptions* blossom_options = nullptr;
    CapacityBlossomPool blossom_pool;
    CapacityBlossomStats* blossom_stats = nullptr;
    CutManagementState blossom_state;
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

// Applies one round of every enabled family to the fractional point y and
// returns the rows to add. Also updates statistics and pools.
struct RoundOutput {
    std::vector<PstepValidInequalityRow> rows;
    std::size_t cover_count = 0;
    std::size_t blossom_count = 0;
};

RoundOutput separate_round(
    CutContext& context,
    const std::vector<double>& y_values,
    bool is_root,
    double node_count,
    double runtime_seconds
) {
    RoundOutput output;

    if (cover_enabled(context)) {
        const CutManagementOptions& management = context.cover_options->management;
        CutManagementStats& stats = context.cover_stats->management;
        const std::size_t budget = cut_management_budget(
            management, context.cover_state, stats, is_root, node_count, runtime_seconds);
        if (budget > 0U) {
            const auto start = std::chrono::steady_clock::now();
            double best_violation = 0.0;
            std::size_t overlap_rejections = 0;
            const std::vector<std::size_t> indices = screen_type_travel_time_cover_rows(
                *context.cover_rows,
                y_values,
                context.cover_added,
                management.min_violation,
                management.max_overlap_jaccard,
                budget,
                best_violation,
                overlap_rejections
            );
            stats.overlap_rejections += overlap_rejections;
            for (std::size_t index : indices) {
                context.cover_added[index] = true;
                output.rows.push_back((*context.cover_rows)[index].row);
                ++output.cover_count;
            }
            cut_management_finish_round(
                management, context.cover_state, stats, is_root, indices.size(), best_violation);
            stats.separation_seconds +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        }
    }

    if (blossom_enabled(context)) {
        const CutManagementOptions& management = context.blossom_options->management;
        CutManagementStats& stats = context.blossom_stats->management;
        const std::size_t budget = cut_management_budget(
            management, context.blossom_state, stats, is_root, node_count, runtime_seconds);
        if (budget > 0U) {
            double best_violation = 0.0;
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
                best_violation = std::max(best_violation, cut.violation);
            }
            output.blossom_count = cuts.size();
            cut_management_finish_round(
                management, context.blossom_state, stats, is_root, cuts.size(), best_violation);
        }
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

        // The scope and the tree schedule are part of this test, so a node the
        // families skip is left before the node relaxation is copied.
        const bool want_cover = cover_enabled(context_) &&
            !cut_management_exhausted(
                context_.cover_options->management,
                context_.cover_state,
                context_.cover_stats->management,
                is_root,
                node_count) &&
            std::find(context_.cover_added.begin(), context_.cover_added.end(), false) !=
                context_.cover_added.end();
        const bool want_blossom = blossom_enabled(context_) &&
            !cut_management_exhausted(
                context_.blossom_options->management,
                context_.blossom_state,
                context_.blossom_stats->management,
                is_root,
                node_count);
        if (!want_cover && !want_blossom) {
            return;
        }

        std::vector<double> values(y_vars_.size(), 0.0);
        double* raw = getNodeRel(y_vars_.data(), static_cast<int>(y_vars_.size()));
        for (std::size_t i = 0; i < values.size(); ++i) {
            values[i] = raw[i];
        }
        delete[] raw;

        const RoundOutput output =
            separate_round(context_, values, is_root, node_count, runtime_seconds);
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
    core_options.lp_method = options.mip_method;
    core_options.node_method = options.mip_node_method;
    core_options.crossover = options.mip_crossover;
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
            ++result.type_cover_stats.management.cuts_added;
            ++result.type_cover_stats.management.root_cuts_added;
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
    const bool any_dynamic_family = cover_enabled(context) || blossom_enabled(context);

    result.mip_time_limit_used = options.solver_time_limit;
    if (options.root_cut_mode == RootCutMode::PreMip && any_dynamic_family &&
        options.model_type == DirectTwoIndexModelType::IP) {
        const auto phase_start = std::chrono::steady_clock::now();
        TwoIndexCoreOptions lp_options = core_options;
        lp_options.binary_y = false;
        lp_options.output_enabled = true;
        lp_options.gurobi_log_path = options.pre_mip_lp_gurobi_log_path;
        lp_options.lp_method = options.pre_mip_method;
        lp_options.node_method = -1;
        lp_options.crossover = options.pre_mip_crossover;
        lp_options.solver_time_limit = options.pre_mip_root_lp_time_limit;
        lp_options.name_prefix = "root_lp_phase";
        TwoIndexCoreModel lp = build_two_index_model_core(data, graph, lp_options);
        for (const PstepValidInequalityRow& row : vi_rows) {
            lp.model->addConstr(row_to_temp_constr(row, lp.y_vars), row.name);
        }
        if (options.type_cover.enabled &&
            options.type_cover.mode == TypeTravelTimeCoverMode::Static) {
            for (const TypeTravelTimeCoverRow& row : options.type_cover_rows) {
                lp.model->addConstr(row_to_temp_constr(row.row, lp.y_vars), row.row.name);
            }
        }
        lp.model->optimize();
        std::vector<PstepValidInequalityRow> phase_rows;
        std::vector<double> values(lp.y_vars.size(), 0.0);
        const auto elapsed = [&]() {
            return std::chrono::duration<double>(
                std::chrono::steady_clock::now() - phase_start).count();
        };
        if (lp.model->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            result.pre_mip_lp_value_before = lp.model->get(GRB_DoubleAttr_ObjVal);
        }
        // Each family owns its root round/cut caps and tailing-off rules, so the
        // loop only guards the phase budget and a hard safety bound.
        for (int round = 0; round < 100000; ++round) {
            if (lp.model->get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
                break;
            }
            if (options.pre_mip_root_lp_time_limit > 0.0 &&
                elapsed() >= options.pre_mip_root_lp_time_limit) {
                break;
            }
            for (std::size_t i = 0; i < values.size(); ++i) {
                values[i] = lp.y_vars[i].get(GRB_DoubleAttr_X);
            }
            const RoundOutput output =
                separate_round(context, values, true, 0.0, elapsed());
            if (output.rows.empty()) {
                break;
            }
            for (const PstepValidInequalityRow& row : output.rows) {
                lp.model->addConstr(row_to_temp_constr(row, lp.y_vars), row.name);
                phase_rows.push_back(row);
            }
            ++result.pre_mip_rounds;
            if (options.pre_mip_root_lp_time_limit > 0.0) {
                const double remaining = options.pre_mip_root_lp_time_limit - elapsed();
                if (remaining <= 0.0) {
                    break;
                }
                lp.model->set(GRB_DoubleParam_TimeLimit, remaining);
            }
            lp.model->optimize();
        }
        if (lp.model->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            result.pre_mip_lp_value_after = lp.model->get(GRB_DoubleAttr_ObjVal);
        }
        for (const PstepValidInequalityRow& row : phase_rows) {
            core.model->addConstr(row_to_temp_constr(row, core.y_vars), row.name);
        }
        result.pre_mip_rows = static_cast<int>(phase_rows.size());
        result.pre_mip_seconds = elapsed();
        // The root phase already consumed part of the instance budget.
        if (options.solver_time_limit > 0.0) {
            result.mip_time_limit_used = std::max(
                0.0, options.solver_time_limit - result.pre_mip_seconds);
            core.model->set(GRB_DoubleParam_TimeLimit, result.mip_time_limit_used);
        }
        // By default the pre-MIP phase owns the root: both families are closed
        // there and the callback may only separate inside the tree. When the
        // caller asks for it, the root phase is reopened so the main MIP root
        // continues the same cut loop on its own relaxation. Already added rows
        // stay recorded, so no row is separated twice.
        if (options.pre_mip_separate_main_mip_root) {
            if (context.cover_stats != nullptr) {
                cut_management_reopen_root(
                    context.cover_state, context.cover_stats->management);
            }
            if (context.blossom_stats != nullptr) {
                cut_management_reopen_root(
                    context.blossom_state, context.blossom_stats->management);
            }
        } else {
            context.cover_state.root_closed = true;
            context.blossom_state.root_closed = true;
        }
        core.model->update();
        result.constraint_count = core.model->get(GRB_IntAttr_NumConstrs);
    }

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
            const RoundOutput output = separate_round(context, values, true, 0.0, 0.0);
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

const char* root_cut_mode_name(RootCutMode mode) {
    switch (mode) {
        case RootCutMode::PreMip: return "pre-mip";
        case RootCutMode::Callback: return "callback";
    }
    return "unknown";
}

bool parse_root_cut_mode(const std::string& value, RootCutMode& mode) {
    if (value == "pre-mip") {
        mode = RootCutMode::PreMip;
        return true;
    }
    if (value == "callback") {
        mode = RootCutMode::Callback;
        return true;
    }
    return false;
}

}  // namespace spdp
