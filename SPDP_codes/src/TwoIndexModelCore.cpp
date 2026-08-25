#include "TwoIndexModelCore.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace spdp::detail {
namespace {

constexpr double kTolerance = 1e-9;
constexpr double kAuxiliaryLPSolverTolerance = 1e-9;

State canonical_auxiliary_state(State state) {
    if (state[1] < state[0]) {
        std::swap(state[0], state[1]);
    }
    return state;
}

double two_index_objective_coefficient(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const EdgeRecord& edge,
    TwoIndexCoreObjective objective
) {
    switch (objective) {
        case TwoIndexCoreObjective::OriginalCost:
            return edge.data.cost;
        case TwoIndexCoreObjective::TravelCost: {
            const bool is_vehicle_departure =
                edge.u == 0 &&
                graph.node(edge.v).kind == NodeSpec::Kind::Pickup;
            return edge.data.cost -
                (is_vehicle_departure ? data.fixed_vehicle_cost : 0.0);
        }
        case TwoIndexCoreObjective::Duration:
            return edge.data.time;
        case TwoIndexCoreObjective::DurationPlusFixed: {
            const bool is_vehicle_departure =
                edge.u == 0 &&
                graph.node(edge.v).kind == NodeSpec::Kind::Pickup;
            return edge.data.time +
                (is_vehicle_departure ? data.fixed_vehicle_cost : 0.0);
        }
    }
    throw std::runtime_error("Unsupported two-index objective.");
}

void validate_two_index_inputs(
    const SPDPData& data,
    const MultiDiGraph& graph,
    double solver_time_limit
) {
    if (data.time_limit <= 0.0) {
        throw std::runtime_error(
            "The two-index formulation requires a positive route-duration limit."
        );
    }
    if (!std::isfinite(solver_time_limit) || solver_time_limit < 0.0) {
        throw std::runtime_error(
            "The two-index formulation time limit must be finite and nonnegative."
        );
    }
    if (graph.number_of_nodes() != 2U * data.requests.size() + 2U) {
        throw std::runtime_error(
            "The two-index formulation requires exactly two action nodes per request."
        );
    }

    const NodeId start_node_id = 0;
    const NodeId end_node_id = graph.end_node_id();
    if (graph.node(start_node_id).kind != NodeSpec::Kind::Start ||
        graph.node(end_node_id).kind != NodeSpec::Kind::End) {
        throw std::runtime_error(
            "The two-index formulation requires node 0 and the last node to be the depots."
        );
    }

    std::size_t dummy_edge_count = 0U;
    for (const EdgeRecord& edge : graph.edges()) {
        if (!std::isfinite(edge.data.time) || edge.data.time < -kTolerance ||
            !std::isfinite(edge.data.cost)) {
            throw std::runtime_error(
                "The two-index formulation requires finite costs and nonnegative durations."
            );
        }
        if (edge.u == start_node_id && edge.v == end_node_id) {
            ++dummy_edge_count;
            if (std::abs(edge.data.time) > kTolerance) {
                throw std::runtime_error(
                    "The two-index formulation requires zero-duration dummy depot edges."
                );
            }
        }
        if (edge.u == start_node_id && edge.v != end_node_id &&
            graph.node(edge.v).kind != NodeSpec::Kind::Pickup) {
            throw std::runtime_error(
                "A two-index departure edge must enter a pickup node."
            );
        }
        if (edge.v == end_node_id && edge.u != start_node_id &&
            graph.node(edge.u).kind != NodeSpec::Kind::Delivery) {
            throw std::runtime_error(
                "A two-index return edge must leave a delivery node."
            );
        }
    }
    if (dummy_edge_count == 0U) {
        throw std::runtime_error(
            "The two-index formulation graph has no start-to-end dummy edge."
        );
    }
}

}  // namespace

TwoIndexCoreModel build_two_index_model_core(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const TwoIndexCoreOptions& options
) {
    validate_two_index_inputs(data, graph, options.solver_time_limit);

    TwoIndexCoreModel core;
    core.environment = std::make_unique<GRBEnv>(true);
    core.environment->set(
        GRB_IntParam_OutputFlag,
        options.output_enabled ? 1 : 0
    );
    if (!options.gurobi_log_path.empty()) {
        core.environment->set(GRB_StringParam_LogFile, options.gurobi_log_path);
    }
    core.environment->start();
    core.model = std::make_unique<GRBModel>(*core.environment);
    core.model->set(
        GRB_IntParam_OutputFlag,
        options.output_enabled ? 1 : 0
    );
    if (options.gurobi_threads >= 0) {
        core.model->set(GRB_IntParam_Threads, options.gurobi_threads);
    }
    core.model->set(GRB_DoubleParam_FeasibilityTol, kAuxiliaryLPSolverTolerance);
    core.model->set(GRB_DoubleParam_OptimalityTol, kAuxiliaryLPSolverTolerance);
    if (options.solver_time_limit > 0.0) {
        core.model->set(GRB_DoubleParam_TimeLimit, options.solver_time_limit);
    }

    const NodeId start_node_id = 0;
    const NodeId end_node_id = graph.end_node_id();
    core.y_vars.reserve(graph.number_of_edges());
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        const bool is_dummy = edge.u == start_node_id && edge.v == end_node_id;
        const double objective_coefficient = two_index_objective_coefficient(
            data,
            graph,
            edge,
            options.objective
        );
        core.y_vars.push_back(core.model->addVar(
            0.0,
            is_dummy ? 0.0 : 1.0,
            objective_coefficient,
            options.binary_y ? GRB_BINARY : GRB_CONTINUOUS,
            options.name_prefix + "_y_" + std::to_string(edge_id)
        ));
    }
    core.model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    core.model->update();

    std::vector<std::set<State>> states_by_node(graph.number_of_nodes());
    for (NodeId node_id = 1; node_id < end_node_id; ++node_id) {
        if (!graph.is_physical_service_node(node_id)) {
            continue;
        }
        const NodeSpec::Kind node_kind = graph.node(node_id).kind;
        if (node_kind != NodeSpec::Kind::Pickup &&
            node_kind != NodeSpec::Kind::Delivery) {
            throw std::runtime_error(
                "The two-index formulation found a non-service physical node."
            );
        }

        GRBLinExpr incoming_degree = 0.0;
        for (std::size_t edge_id : graph.ingoing_edge_indices(node_id)) {
            incoming_degree += core.y_vars[edge_id];
        }
        core.model->addConstr(
            incoming_degree == 1.0,
            options.name_prefix + "_in_degree_" + std::to_string(node_id)
        );

        GRBLinExpr outgoing_degree = 0.0;
        for (std::size_t edge_id : graph.outgoing_edge_indices(node_id)) {
            outgoing_degree += core.y_vars[edge_id];
        }
        core.model->addConstr(
            outgoing_degree == 1.0,
            options.name_prefix + "_out_degree_" + std::to_string(node_id)
        );

        for (std::size_t edge_id : graph.ingoing_edge_indices(node_id)) {
            states_by_node[static_cast<std::size_t>(node_id)].insert(
                canonical_auxiliary_state(graph.edges()[edge_id].data.end_state)
            );
        }
        for (std::size_t edge_id : graph.outgoing_edge_indices(node_id)) {
            states_by_node[static_cast<std::size_t>(node_id)].insert(
                canonical_auxiliary_state(graph.edges()[edge_id].data.start_state)
            );
        }

        std::size_t state_index = 0U;
        for (const State& state : states_by_node[static_cast<std::size_t>(node_id)]) {
            GRBLinExpr state_balance = 0.0;
            for (std::size_t edge_id : graph.ingoing_edge_indices(node_id)) {
                if (canonical_auxiliary_state(graph.edges()[edge_id].data.end_state) ==
                    state) {
                    state_balance += core.y_vars[edge_id];
                }
            }
            for (std::size_t edge_id : graph.outgoing_edge_indices(node_id)) {
                if (canonical_auxiliary_state(graph.edges()[edge_id].data.start_state) ==
                    state) {
                    state_balance -= core.y_vars[edge_id];
                }
            }
            core.model->addConstr(
                state_balance == 0.0,
                options.name_prefix + "_state_" + std::to_string(node_id) + "_" +
                    std::to_string(state_index++)
            );
        }
    }

    if (options.add_time_constraints) {
        std::vector<std::map<State, GRBVar>> time_vars_by_node(
            graph.number_of_nodes()
        );
        for (NodeId node_id = 1; node_id < end_node_id; ++node_id) {
            std::size_t state_index = 0U;
            for (const State& state :
                 states_by_node[static_cast<std::size_t>(node_id)]) {
                time_vars_by_node[static_cast<std::size_t>(node_id)].emplace(
                    state,
                    core.model->addVar(
                        0.0,
                        data.time_limit,
                        0.0,
                        GRB_CONTINUOUS,
                        options.name_prefix + "_B_" + std::to_string(node_id) + "_" +
                            std::to_string(state_index++)
                    )
                );
            }
        }
        core.model->update();

        const auto time_var = [&](NodeId node_id, const State& state) -> GRBVar {
            const State canonical_state = canonical_auxiliary_state(state);
            const auto& variables =
                time_vars_by_node.at(static_cast<std::size_t>(node_id));
            const auto found = variables.find(canonical_state);
            if (found == variables.end()) {
                throw std::runtime_error(
                    "The two-index formulation found no time variable for an edge "
                    "endpoint state."
                );
            }
            return found->second;
        };

        for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
            const EdgeRecord& edge = graph.edges()[edge_id];
            if (edge.u == start_node_id && edge.v == end_node_id) {
                continue;
            }

            const std::string constraint_name =
                options.name_prefix + "_time_" + std::to_string(edge_id);
            if (edge.u == start_node_id) {
                const GRBVar end_time = time_var(edge.v, edge.data.end_state);
                core.model->addConstr(
                    end_time >= edge.data.time * core.y_vars[edge_id],
                    constraint_name
                );
                continue;
            }
            if (edge.v == end_node_id) {
                const GRBVar start_time = time_var(edge.u, edge.data.start_state);
                core.model->addConstr(
                    start_time + edge.data.time * core.y_vars[edge_id] <=
                        data.time_limit,
                    constraint_name
                );
                continue;
            }

            const GRBVar start_time = time_var(edge.u, edge.data.start_state);
            const GRBVar end_time = time_var(edge.v, edge.data.end_state);
            const double big_m = data.time_limit + edge.data.time;
            core.model->addConstr(
                end_time >= start_time + edge.data.time -
                    big_m * (1.0 - core.y_vars[edge_id]),
                constraint_name
            );
        }
    }

    GRBLinExpr actual_departures = 0.0;
    GRBLinExpr actual_returns = 0.0;
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        if (edge.u == start_node_id &&
            graph.node(edge.v).kind == NodeSpec::Kind::Pickup) {
            actual_departures += core.y_vars[edge_id];
        }
        if (graph.node(edge.u).kind == NodeSpec::Kind::Delivery &&
            edge.v == end_node_id) {
            actual_returns += core.y_vars[edge_id];
        }
    }
    core.model->addConstr(
        actual_departures == actual_returns,
        options.name_prefix + "_depot_balance"
    );
    core.model->update();
    return core;
}

}  // namespace spdp::detail
