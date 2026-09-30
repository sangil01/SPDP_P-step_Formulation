#include "TwoIndexModelCore.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <utility>

namespace spdp::detail {
namespace {

constexpr double kTolerance = 1e-9;

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
        case TwoIndexCoreObjective::None:
            return 0.0;
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

    const double time_horizon = options.time_horizon > 0.0
        ? options.time_horizon
        : data.time_limit;
    if (options.add_time_constraints &&
        (!std::isfinite(time_horizon) || time_horizon <= 0.0)) {
        throw std::runtime_error(
            "The two-index time horizon must be finite and positive."
        );
    }

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
    core.model->set(GRB_DoubleParam_FeasibilityTol, kGurobiSolverTolerance);
    core.model->set(GRB_DoubleParam_OptimalityTol, kGurobiSolverTolerance);
    if (options.lp_method >= -1) {
        core.model->set(GRB_IntParam_Method, options.lp_method);
    }
    if (options.node_method >= -1) {
        core.model->set(GRB_IntParam_NodeMethod, options.node_method);
    }
    if (options.crossover >= 0) {
        core.model->set(GRB_IntParam_Crossover, options.crossover);
    }
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
        core.time_vars_by_node.assign(graph.number_of_nodes(), {});
        for (NodeId node_id = 1; node_id < end_node_id; ++node_id) {
            std::size_t state_index = 0U;
            for (const State& state :
                 states_by_node[static_cast<std::size_t>(node_id)]) {
                core.time_vars_by_node[static_cast<std::size_t>(node_id)].emplace(
                    state,
                    core.model->addVar(
                        0.0,
                        time_horizon,
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
                core.time_vars_by_node.at(static_cast<std::size_t>(node_id));
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
                if (!options.enforce_route_duration_limit) {
                    continue;
                }
                const GRBVar start_time = time_var(edge.u, edge.data.start_state);
                core.model->addConstr(
                    start_time + edge.data.time * core.y_vars[edge_id] <=
                        time_horizon,
                    constraint_name
                );
                continue;
            }

            const GRBVar start_time = time_var(edge.u, edge.data.start_state);
            const GRBVar end_time = time_var(edge.v, edge.data.end_state);
            const double big_m = time_horizon + edge.data.time;
            core.model->addConstr(
                end_time >= start_time + edge.data.time -
                    big_m * (1.0 - core.y_vars[edge_id]),
                constraint_name
            );
        }
    }

    if (options.add_time_flow_formulation) {
        // Node-state shortest duration labels (relaxing node uniqueness):
        // h0[(v,s)] = min duration from the depot to (v, s), including service at v;
        // hd[(v,s)] = min duration from (v, s) to the end depot.
        using Key = std::pair<NodeId, State>;
        std::map<Key, int> index_of;
        std::vector<Key> keys;
        const auto key_index = [&](NodeId node_id, const State& state) {
            const Key key{node_id, canonical_auxiliary_state(state)};
            auto found = index_of.find(key);
            if (found == index_of.end()) {
                found = index_of.emplace(key, static_cast<int>(keys.size())).first;
                keys.push_back(key);
            }
            return found->second;
        };
        struct Arc { int to; double weight; };
        std::vector<std::vector<Arc>> forward_arcs;
        std::vector<std::vector<Arc>> backward_arcs;
        std::vector<int> edge_from(graph.number_of_edges(), -1);
        std::vector<int> edge_to(graph.number_of_edges(), -1);
        for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
            const EdgeRecord& edge = graph.edges()[edge_id];
            if (edge.u == start_node_id && edge.v == end_node_id) {
                continue;
            }
            const int from = key_index(edge.u, edge.data.start_state);
            const int to = key_index(edge.v, edge.data.end_state);
            edge_from[edge_id] = from;
            edge_to[edge_id] = to;
            if (static_cast<int>(forward_arcs.size()) <= std::max(from, to)) {
                forward_arcs.resize(static_cast<std::size_t>(std::max(from, to)) + 1);
                backward_arcs.resize(forward_arcs.size());
            }
            forward_arcs[static_cast<std::size_t>(from)].push_back({to, edge.data.time});
            backward_arcs[static_cast<std::size_t>(to)].push_back({from, edge.data.time});
        }
        forward_arcs.resize(keys.size());
        backward_arcs.resize(keys.size());
        const double infinity = std::numeric_limits<double>::infinity();
        const auto dijkstra = [&](const std::vector<std::vector<Arc>>& arcs,
                                  const std::vector<int>& sources) {
            std::vector<double> dist(keys.size(), infinity);
            using Item = std::pair<double, int>;
            std::priority_queue<Item, std::vector<Item>, std::greater<Item>> heap;
            for (int s : sources) {
                dist[static_cast<std::size_t>(s)] = 0.0;
                heap.emplace(0.0, s);
            }
            while (!heap.empty()) {
                const auto [d, u] = heap.top();
                heap.pop();
                if (d > dist[static_cast<std::size_t>(u)] + 1e-12) {
                    continue;
                }
                for (const Arc& arc : arcs[static_cast<std::size_t>(u)]) {
                    const double nd = d + arc.weight;
                    if (nd + 1e-12 < dist[static_cast<std::size_t>(arc.to)]) {
                        dist[static_cast<std::size_t>(arc.to)] = nd;
                        heap.emplace(nd, arc.to);
                    }
                }
            }
            return dist;
        };
        std::vector<int> depot_keys;
        std::vector<int> end_keys;
        for (std::size_t i = 0; i < keys.size(); ++i) {
            if (keys[i].first == start_node_id) depot_keys.push_back(static_cast<int>(i));
            if (keys[i].first == end_node_id) end_keys.push_back(static_cast<int>(i));
        }
        const std::vector<double> h0 = dijkstra(forward_arcs, depot_keys);
        const std::vector<double> hd = dijkstra(backward_arcs, end_keys);

        core.time_flow_vars.reserve(graph.number_of_edges());
        std::vector<double> lower_coefficient(graph.number_of_edges(), 0.0);
        std::vector<double> upper_coefficient(graph.number_of_edges(), 0.0);
        for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
            const EdgeRecord& edge = graph.edges()[edge_id];
            const bool is_dummy = edge.u == start_node_id && edge.v == end_node_id;
            double lower = 0.0;
            double upper = 0.0;
            if (!is_dummy) {
                const double prefix = h0[static_cast<std::size_t>(edge_from[edge_id])];
                const double suffix = hd[static_cast<std::size_t>(edge_to[edge_id])];
                lower = (std::isfinite(prefix) ? prefix : 0.0) + edge.data.time;
                upper = time_horizon - (std::isfinite(suffix) ? suffix : 0.0);
                if (lower > upper + 1e-9) {
                    // No duration-feasible route can use this edge.
                    core.y_vars[edge_id].set(GRB_DoubleAttr_UB, 0.0);
                    ++core.time_flow_infeasible_edge_count;
                    upper = std::max(upper, lower);
                }
            }
            lower_coefficient[edge_id] = lower;
            upper_coefficient[edge_id] = upper;
            core.time_flow_vars.push_back(core.model->addVar(
                0.0,
                is_dummy ? 0.0 : std::max(upper, 0.0),
                0.0,
                GRB_CONTINUOUS,
                options.name_prefix + "_g_" + std::to_string(edge_id)
            ));
        }
        core.model->update();
        for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
            const EdgeRecord& edge = graph.edges()[edge_id];
            if (edge.u == start_node_id && edge.v == end_node_id) {
                continue;
            }
            const std::string name =
                options.name_prefix + "_tf_" + std::to_string(edge_id);
            if (edge.u == start_node_id) {
                core.model->addConstr(
                    core.time_flow_vars[edge_id] == edge.data.time * core.y_vars[edge_id],
                    name + "_depot"
                );
                continue;
            }
            core.model->addConstr(
                core.time_flow_vars[edge_id] <=
                    upper_coefficient[edge_id] * core.y_vars[edge_id],
                name + "_ub"
            );
            core.model->addConstr(
                core.time_flow_vars[edge_id] >=
                    lower_coefficient[edge_id] * core.y_vars[edge_id],
                name + "_lb"
            );
        }
        for (NodeId node_id = 1; node_id < end_node_id; ++node_id) {
            if (!graph.is_physical_service_node(node_id)) {
                continue;
            }
            if (options.time_flow_state_disaggregated) {
                std::size_t state_index = 0U;
                for (const State& state :
                     states_by_node[static_cast<std::size_t>(node_id)]) {
                    GRBLinExpr outgoing = 0.0;
                    GRBLinExpr incoming = 0.0;
                    for (std::size_t edge_id : graph.outgoing_edge_indices(node_id)) {
                        if (canonical_auxiliary_state(
                                graph.edges()[edge_id].data.start_state) != state) {
                            continue;
                        }
                        outgoing += core.time_flow_vars[edge_id];
                        outgoing -= graph.edges()[edge_id].data.time * core.y_vars[edge_id];
                    }
                    for (std::size_t edge_id : graph.ingoing_edge_indices(node_id)) {
                        if (canonical_auxiliary_state(
                                graph.edges()[edge_id].data.end_state) != state) {
                            continue;
                        }
                        incoming += core.time_flow_vars[edge_id];
                    }
                    core.model->addConstr(
                        outgoing == incoming,
                        options.name_prefix + "_tf_balance_" + std::to_string(node_id) +
                            "_" + std::to_string(state_index++)
                    );
                }
                continue;
            }
            GRBLinExpr outgoing = 0.0;
            GRBLinExpr incoming = 0.0;
            for (std::size_t edge_id : graph.outgoing_edge_indices(node_id)) {
                outgoing += core.time_flow_vars[edge_id];
                outgoing -= graph.edges()[edge_id].data.time * core.y_vars[edge_id];
            }
            for (std::size_t edge_id : graph.ingoing_edge_indices(node_id)) {
                incoming += core.time_flow_vars[edge_id];
            }
            core.model->addConstr(
                outgoing == incoming,
                options.name_prefix + "_tf_balance_" + std::to_string(node_id)
            );
        }
        core.model->update();
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

GRBVar two_index_time_var(
    const TwoIndexCoreModel& core,
    NodeId node_id,
    const State& state
) {
    const State canonical_state = canonical_auxiliary_state(state);
    const auto& variables =
        core.time_vars_by_node.at(static_cast<std::size_t>(node_id));
    const auto found = variables.find(canonical_state);
    if (found == variables.end()) {
        throw std::runtime_error(
            "The two-index formulation found no time variable for an endpoint state."
        );
    }
    return found->second;
}

}  // namespace spdp::detail
