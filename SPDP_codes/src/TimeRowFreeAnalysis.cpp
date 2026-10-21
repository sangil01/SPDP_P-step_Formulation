#include "TimeRowFreeAnalysis.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace spdp {
namespace {

constexpr double kTolerance = 1e-6;

State canonical_state(State state) {
    if (state[1] < state[0]) {
        std::swap(state[0], state[1]);
    }
    return state;
}

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

std::string node_degree_message(
    NodeId node_id,
    int incoming,
    int outgoing
) {
    std::ostringstream out;
    out << "Time-row-free solution has invalid selected degree at node "
        << node_id << ": incoming=" << incoming
        << ", outgoing=" << outgoing << '.';
    return out.str();
}

}  // namespace

TimeRowFreeAnalysis analyze_time_row_free_solution(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const std::vector<int>& selected_edge_ids
) {
    require(data.time_limit > 0.0, "Time-row-free analysis requires a positive time limit.");
    require(graph.number_of_nodes() >= 2U, "Time-row-free analysis requires both depots.");

    const NodeId start_node_id = 0;
    const NodeId end_node_id = graph.end_node_id();
    const std::vector<EdgeRecord>& edges = graph.edges();
    std::vector<bool> selected(edges.size(), false);
    std::vector<bool> consumed(edges.size(), false);
    std::vector<int> incoming(graph.number_of_nodes(), 0);
    std::vector<int> outgoing(graph.number_of_nodes(), 0);
    std::vector<int> selected_incoming(graph.number_of_nodes(), -1);
    std::vector<int> selected_outgoing(graph.number_of_nodes(), -1);
    std::vector<int> depot_departures;

    for (int edge_id : selected_edge_ids) {
        require(edge_id >= 0, "A selected time-row-free edge id is negative.");
        require(
            static_cast<std::size_t>(edge_id) < edges.size(),
            "A selected time-row-free edge id is out of range."
        );
        require(
            !selected[static_cast<std::size_t>(edge_id)],
            "A time-row-free edge was selected more than once."
        );
        selected[static_cast<std::size_t>(edge_id)] = true;

        const EdgeRecord& edge = edges[static_cast<std::size_t>(edge_id)];
        require(
            edge.u >= start_node_id && edge.u <= end_node_id &&
                edge.v >= start_node_id && edge.v <= end_node_id,
            "A selected time-row-free edge references an invalid node."
        );
        require(
            !(edge.u == start_node_id && edge.v == end_node_id),
            "The disabled depot dummy edge cannot belong to a time-row-free solution."
        );

        ++outgoing[static_cast<std::size_t>(edge.u)];
        ++incoming[static_cast<std::size_t>(edge.v)];
        if (edge.u == start_node_id) {
            depot_departures.push_back(edge_id);
        } else {
            selected_outgoing[static_cast<std::size_t>(edge.u)] = edge_id;
        }
        if (edge.v != end_node_id) {
            selected_incoming[static_cast<std::size_t>(edge.v)] = edge_id;
        }
    }

    require(incoming.front() == 0, "A time-row-free solution enters the start depot.");
    require(
        outgoing[static_cast<std::size_t>(end_node_id)] == 0,
        "A time-row-free solution leaves the end depot."
    );
    require(
        outgoing.front() == incoming[static_cast<std::size_t>(end_node_id)],
        "Time-row-free depot departure and return counts differ."
    );

    for (NodeId node_id = 1; node_id < end_node_id; ++node_id) {
        const int in_degree = incoming[static_cast<std::size_t>(node_id)];
        const int out_degree = outgoing[static_cast<std::size_t>(node_id)];
        require(
            in_degree == 1 && out_degree == 1,
            node_degree_message(node_id, in_degree, out_degree)
        );
        const EdgeRecord& in_edge = edges[static_cast<std::size_t>(
            selected_incoming[static_cast<std::size_t>(node_id)])];
        const EdgeRecord& out_edge = edges[static_cast<std::size_t>(
            selected_outgoing[static_cast<std::size_t>(node_id)])];
        require(
            canonical_state(in_edge.data.end_state) ==
                canonical_state(out_edge.data.start_state),
            "Time-row-free solution violates load-state continuity."
        );
    }

    std::sort(depot_departures.begin(), depot_departures.end());
    TimeRowFreeAnalysis analysis;
    analysis.selected_edge_count = selected_edge_ids.size();

    const auto consume_component = [&](int first_edge_id, bool depot_route) {
        TimeRowFreeComponent component;
        int current_edge_id = first_edge_id;
        const NodeId cycle_start = edges[static_cast<std::size_t>(first_edge_id)].u;

        for (std::size_t step = 0; step <= selected_edge_ids.size(); ++step) {
            require(
                selected[static_cast<std::size_t>(current_edge_id)],
                "Time-row-free decomposition reached an unselected edge."
            );
            require(
                !consumed[static_cast<std::size_t>(current_edge_id)],
                "Time-row-free decomposition assigned an edge more than once."
            );
            consumed[static_cast<std::size_t>(current_edge_id)] = true;
            component.edge_ids.push_back(current_edge_id);

            const EdgeRecord& edge = edges[static_cast<std::size_t>(current_edge_id)];
            component.duration += edge.data.time;
            if (depot_route && edge.v == end_node_id) {
                return component;
            }
            if (!depot_route && edge.v == cycle_start) {
                return component;
            }
            require(
                edge.v != start_node_id && edge.v != end_node_id,
                depot_route
                    ? "A depot-connected time-row-free walk does not end at the end depot."
                    : "A disconnected time-row-free component reaches a depot."
            );
            current_edge_id = selected_outgoing[static_cast<std::size_t>(edge.v)];
            require(
                current_edge_id >= 0,
                "Time-row-free decomposition found no selected continuation."
            );
        }
        throw std::runtime_error("Time-row-free decomposition did not terminate.");
    };

    for (int edge_id : depot_departures) {
        analysis.depot_routes.push_back(consume_component(edge_id, true));
    }
    for (int edge_id : selected_edge_ids) {
        if (!consumed[static_cast<std::size_t>(edge_id)]) {
            analysis.disconnected_subtours.push_back(
                consume_component(edge_id, false)
            );
        }
    }

    std::size_t consumed_count = 0U;
    for (std::size_t edge_id = 0; edge_id < selected.size(); ++edge_id) {
        if (selected[edge_id]) {
            require(
                consumed[edge_id],
                "A selected time-row-free edge was not assigned to a component."
            );
            ++consumed_count;
        }
    }
    require(
        consumed_count == selected_edge_ids.size(),
        "Time-row-free component edge count does not match the selected support."
    );

    for (const TimeRowFreeComponent& route : analysis.depot_routes) {
        analysis.maximum_depot_route_duration =
            std::max(analysis.maximum_depot_route_duration, route.duration);
        const double excess = std::max(0.0, route.duration - data.time_limit);
        analysis.maximum_duration_excess =
            std::max(analysis.maximum_duration_excess, excess);
        if (excess > kTolerance) {
            ++analysis.duration_violating_route_count;
        }
    }
    return analysis;
}

}  // namespace spdp
