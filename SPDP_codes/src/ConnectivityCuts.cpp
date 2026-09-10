#include "ConnectivityCuts.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <set>
#include <stdexcept>

namespace spdp {
namespace {

// Edmonds-Karp on a dense capacity matrix. Node counts are small (2n + 2).
double max_flow_dense(
    std::vector<std::vector<double>>& capacity,
    int source,
    int sink
) {
    const int n = static_cast<int>(capacity.size());
    double flow = 0.0;
    std::vector<int> parent(n);
    while (true) {
        std::fill(parent.begin(), parent.end(), -2);
        parent[source] = -1;
        std::deque<int> queue{source};
        while (!queue.empty() && parent[sink] == -2) {
            const int u = queue.front();
            queue.pop_front();
            for (int v = 0; v < n; ++v) {
                if (parent[v] == -2 && capacity[u][v] > 1e-9) {
                    parent[v] = u;
                    queue.push_back(v);
                }
            }
        }
        if (parent[sink] == -2) {
            return flow;
        }
        double bottleneck = 1e300;
        for (int v = sink; parent[v] != -1; v = parent[v]) {
            bottleneck = std::min(bottleneck, capacity[parent[v]][v]);
        }
        for (int v = sink; parent[v] != -1; v = parent[v]) {
            capacity[parent[v]][v] -= bottleneck;
            capacity[v][parent[v]] += bottleneck;
        }
        flow += bottleneck;
    }
}

}  // namespace

const char* connectivity_cut_rhs_mode_name(ConnectivityCutRhsMode mode) {
    switch (mode) {
        case ConnectivityCutRhsMode::One: return "one";
        case ConnectivityCutRhsMode::Duration: return "duration";
    }
    return "unknown";
}

ConnectivitySeparator build_connectivity_separator(
    const SPDPData& data,
    const MultiDiGraph& graph
) {
    ConnectivitySeparator separator;
    separator.end_node_id = graph.end_node_id();
    separator.node_count = graph.number_of_nodes();
    separator.route_time_limit = data.time_limit;
    separator.min_in_edge_time_by_node.assign(separator.node_count, 0.0);
    separator.in_edge_ids_by_node.assign(separator.node_count, {});
    separator.edge_tail.resize(graph.number_of_edges());
    separator.edge_head.resize(graph.number_of_edges());
    std::vector<bool> seen(separator.node_count, false);
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        separator.edge_tail[edge_id] = edge.u;
        separator.edge_head[edge_id] = edge.v;
        const auto v = static_cast<std::size_t>(edge.v);
        separator.in_edge_ids_by_node[v].push_back(edge_id);
        if (edge.v == separator.end_node_id) {
            continue;
        }
        if (!seen[v] || edge.data.time < separator.min_in_edge_time_by_node[v]) {
            separator.min_in_edge_time_by_node[v] = edge.data.time;
            seen[v] = true;
        }
    }
    return separator;
}

std::vector<PstepValidInequalityRow> separate_connectivity_cuts(
    const ConnectivitySeparator& separator,
    const std::vector<double>& y_values,
    const ConnectivityCutOptions& options
) {
    std::vector<PstepValidInequalityRow> cuts;
    if (y_values.size() != separator.edge_tail.size()) {
        throw std::runtime_error(
            "Connectivity separation received a y vector of the wrong size."
        );
    }
    const int n = static_cast<int>(separator.node_count);
    const int depot = 0;
    const int end_node = static_cast<int>(separator.end_node_id);

    // Aggregate support capacities over parallel edges.
    std::vector<std::vector<double>> base_capacity(n, std::vector<double>(n, 0.0));
    for (std::size_t edge_id = 0; edge_id < y_values.size(); ++edge_id) {
        const double value = y_values[edge_id];
        if (value <= options.support_tolerance) {
            continue;
        }
        base_capacity[separator.edge_tail[edge_id]][separator.edge_head[edge_id]] += value;
    }

    std::set<std::vector<int>> emitted_sets;
    std::vector<std::pair<double, std::vector<int>>> candidates;
    for (int sink = 1; sink < end_node; ++sink) {
        std::vector<std::vector<double>> capacity = base_capacity;
        const double flow = max_flow_dense(capacity, depot, sink);

        // Sink side of the minimum cut: nodes not reachable from the depot in
        // the residual graph.
        std::vector<bool> reachable(n, false);
        reachable[depot] = true;
        std::deque<int> queue{depot};
        while (!queue.empty()) {
            const int u = queue.front();
            queue.pop_front();
            for (int v = 0; v < n; ++v) {
                if (!reachable[v] && capacity[u][v] > 1e-9) {
                    reachable[v] = true;
                    queue.push_back(v);
                }
            }
        }
        std::vector<int> sink_side;
        double h_min = 0.0;
        for (int v = 1; v < end_node; ++v) {
            if (!reachable[v]) {
                sink_side.push_back(v);
                h_min += separator.min_in_edge_time_by_node[static_cast<std::size_t>(v)];
            }
        }
        if (sink_side.empty()) {
            continue;
        }
        double rhs = 1.0;
        if (options.rhs_mode == ConnectivityCutRhsMode::Duration &&
            separator.route_time_limit > 0.0) {
            rhs = std::max(
                1.0,
                std::ceil(h_min / separator.route_time_limit - 1e-9)
            );
        }
        // The max-flow value equals the entering support flow of the sink side.
        const double violation = rhs - flow;
        if (violation <= options.violation_tolerance) {
            continue;
        }
        if (!emitted_sets.insert(sink_side).second) {
            continue;
        }
        candidates.emplace_back(violation, sink_side);
    }

    std::sort(
        candidates.begin(),
        candidates.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.first > rhs.first; }
    );
    if (candidates.size() > options.max_cuts_per_round) {
        candidates.resize(options.max_cuts_per_round);
    }

    for (const auto& [violation, sink_side] : candidates) {
        std::vector<bool> in_set(n, false);
        double h_min = 0.0;
        for (int v : sink_side) {
            in_set[v] = true;
            h_min += separator.min_in_edge_time_by_node[static_cast<std::size_t>(v)];
        }
        PstepValidInequalityRow row;
        row.sense = PstepValidInequalitySense::GreaterEqual;
        row.rhs = 1.0;
        if (options.rhs_mode == ConnectivityCutRhsMode::Duration &&
            separator.route_time_limit > 0.0) {
            row.rhs = std::max(
                1.0,
                std::ceil(h_min / separator.route_time_limit - 1e-9)
            );
        }
        row.name = "conn_cut_" + std::to_string(sink_side.size()) + "_" +
            std::to_string(static_cast<int>(row.rhs));
        for (int v : sink_side) {
            for (std::size_t edge_id :
                 separator.in_edge_ids_by_node[static_cast<std::size_t>(v)]) {
                if (!in_set[separator.edge_tail[edge_id]]) {
                    row.edge_terms.emplace_back(edge_id, 1.0);
                }
            }
        }
        cuts.push_back(std::move(row));
    }
    return cuts;
}

}  // namespace spdp
