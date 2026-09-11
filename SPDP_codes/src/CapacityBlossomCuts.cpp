#include "CapacityBlossomCuts.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace spdp {
namespace {

constexpr double kResidualEps = 1e-12;

// Dinic maximum flow on a small dense symmetric capacity matrix (undirected
// graph). After max_flow() the source side of a minimum cut is the set of
// nodes still reachable from the source in the residual network.
class DenseDinic {
public:
    explicit DenseDinic(std::vector<std::vector<double>> capacity)
        : n_(static_cast<int>(capacity.size())),
          capacity_(std::move(capacity)),
          adjacency_(static_cast<std::size_t>(n_)),
          level_(static_cast<std::size_t>(n_), -1),
          iterator_(static_cast<std::size_t>(n_), 0) {
        for (int u = 0; u < n_; ++u) {
            for (int v = 0; v < n_; ++v) {
                if (u != v && (capacity_[u][v] > kResidualEps || capacity_[v][u] > kResidualEps)) {
                    adjacency_[static_cast<std::size_t>(u)].push_back(v);
                }
            }
        }
    }

    double max_flow(int source, int sink) {
        double flow = 0.0;
        while (bfs(source, sink)) {
            std::fill(iterator_.begin(), iterator_.end(), 0);
            while (true) {
                const double pushed = dfs(source, sink, std::numeric_limits<double>::infinity());
                if (pushed <= kResidualEps) {
                    break;
                }
                flow += pushed;
            }
        }
        // The final failed BFS leaves level_ >= 0 exactly on the source side.
        return flow;
    }

    std::vector<bool> source_side() const {
        std::vector<bool> side(static_cast<std::size_t>(n_), false);
        for (int v = 0; v < n_; ++v) {
            side[static_cast<std::size_t>(v)] = level_[static_cast<std::size_t>(v)] >= 0;
        }
        return side;
    }

private:
    bool bfs(int source, int sink) {
        std::fill(level_.begin(), level_.end(), -1);
        level_[static_cast<std::size_t>(source)] = 0;
        std::deque<int> queue{source};
        while (!queue.empty()) {
            const int u = queue.front();
            queue.pop_front();
            for (int v : adjacency_[static_cast<std::size_t>(u)]) {
                if (level_[static_cast<std::size_t>(v)] < 0 && capacity_[u][v] > kResidualEps) {
                    level_[static_cast<std::size_t>(v)] = level_[static_cast<std::size_t>(u)] + 1;
                    queue.push_back(v);
                }
            }
        }
        return level_[static_cast<std::size_t>(sink)] >= 0;
    }

    double dfs(int u, int sink, double limit) {
        if (u == sink) {
            return limit;
        }
        std::size_t& it = iterator_[static_cast<std::size_t>(u)];
        const std::vector<int>& neighbors = adjacency_[static_cast<std::size_t>(u)];
        for (; it < neighbors.size(); ++it) {
            const int v = neighbors[it];
            if (capacity_[u][v] <= kResidualEps ||
                level_[static_cast<std::size_t>(v)] != level_[static_cast<std::size_t>(u)] + 1) {
                continue;
            }
            const double pushed = dfs(v, sink, std::min(limit, capacity_[u][v]));
            if (pushed > kResidualEps) {
                capacity_[u][v] -= pushed;
                capacity_[v][u] += pushed;
                return pushed;
            }
        }
        return 0.0;
    }

    int n_;
    std::vector<std::vector<double>> capacity_;
    std::vector<std::vector<int>> adjacency_;
    std::vector<int> level_;
    std::vector<std::size_t> iterator_;
};

struct CutTree {
    std::vector<int> parent;
    std::vector<double> weight;
};

// Gusfield's algorithm for a Gomory--Hu (cut-equivalent) tree: m - 1 maximum
// flows, no contractions. Node 0 is the root.
CutTree gusfield_cut_tree(
    const std::vector<std::vector<double>>& capacity,
    std::size_t& max_flow_calls
) {
    const int m = static_cast<int>(capacity.size());
    CutTree tree;
    tree.parent.assign(static_cast<std::size_t>(m), 0);
    tree.weight.assign(static_cast<std::size_t>(m), 0.0);
    for (int s = 1; s < m; ++s) {
        const int t = tree.parent[static_cast<std::size_t>(s)];
        DenseDinic flow(capacity);
        const double value = flow.max_flow(s, t);
        ++max_flow_calls;
        const std::vector<bool> source_side = flow.source_side();
        tree.weight[static_cast<std::size_t>(s)] = value;
        for (int i = 0; i < m; ++i) {
            if (i != s && source_side[static_cast<std::size_t>(i)] &&
                tree.parent[static_cast<std::size_t>(i)] == t) {
                tree.parent[static_cast<std::size_t>(i)] = s;
            }
        }
        const int parent_of_t = tree.parent[static_cast<std::size_t>(t)];
        if (t != 0 && source_side[static_cast<std::size_t>(parent_of_t)]) {
            tree.parent[static_cast<std::size_t>(s)] = parent_of_t;
            tree.parent[static_cast<std::size_t>(t)] = s;
            tree.weight[static_cast<std::size_t>(s)] = tree.weight[static_cast<std::size_t>(t)];
            tree.weight[static_cast<std::size_t>(t)] = value;
        }
    }
    return tree;
}

// Component containing 'start' after deleting the tree edge (start, parent[start]).
std::vector<bool> tree_side(const CutTree& tree, int start) {
    const int m = static_cast<int>(tree.parent.size());
    std::vector<std::vector<int>> adjacency(static_cast<std::size_t>(m));
    for (int v = 1; v < m; ++v) {
        if (v == start) {
            continue;
        }
        const int p = tree.parent[static_cast<std::size_t>(v)];
        adjacency[static_cast<std::size_t>(v)].push_back(p);
        adjacency[static_cast<std::size_t>(p)].push_back(v);
    }
    std::vector<bool> side(static_cast<std::size_t>(m), false);
    side[static_cast<std::size_t>(start)] = true;
    std::deque<int> queue{start};
    while (!queue.empty()) {
        const int u = queue.front();
        queue.pop_front();
        for (int v : adjacency[static_cast<std::size_t>(u)]) {
            if (!side[static_cast<std::size_t>(v)]) {
                side[static_cast<std::size_t>(v)] = true;
                queue.push_back(v);
            }
        }
    }
    return side;
}

std::size_t pair_key(std::size_t i, std::size_t j, std::size_t size) {
    if (i > j) {
        std::swap(i, j);
    }
    return i * size + j;
}

void build_family(
    const MultiDiGraph& graph,
    NodeSpec::Kind kind,
    CapacityFamilyData& family
) {
    family.compact_index_by_node.assign(graph.number_of_nodes(), -1);
    for (NodeId node_id = 1; node_id < graph.end_node_id(); ++node_id) {
        if (!graph.is_physical_service_node(node_id) || graph.node(node_id).kind != kind) {
            continue;
        }
        family.compact_index_by_node[static_cast<std::size_t>(node_id)] =
            static_cast<int>(family.nodes.size());
        family.nodes.push_back(node_id);
    }
    const std::size_t size = family.nodes.size();
    family.in_edges_by_compact.assign(size, {});
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        const int head = family.compact_index_by_node[static_cast<std::size_t>(edge.v)];
        if (head < 0) {
            continue;
        }
        family.in_edges_by_compact[static_cast<std::size_t>(head)].push_back(edge_id);
        const int tail = family.compact_index_by_node[static_cast<std::size_t>(edge.u)];
        if (tail < 0 || tail == head) {
            continue;
        }
        family.pair_edges[pair_key(
            static_cast<std::size_t>(tail), static_cast<std::size_t>(head), size
        )].push_back(edge_id);
    }
}

struct FamilyCandidate {
    std::vector<int> compact_set;
    double violation = 0.0;
};

std::vector<FamilyCandidate> separate_family(
    const CapacityFamilyData& family,
    const std::vector<double>& y_values,
    const CapacityBlossomOptions& options,
    CapacityBlossomStats& stats
) {
    std::vector<FamilyCandidate> candidates;
    const std::size_t size = family.nodes.size();
    if (size < 3U) {
        return candidates;
    }

    // Aggregate parallel edges into z_ij and node degrees.
    std::vector<double> degree(size, 0.0);
    std::vector<std::pair<std::pair<std::size_t, std::size_t>, double>> weights;
    for (const auto& entry : family.pair_edges) {
        double weight = 0.0;
        for (std::size_t edge_id : entry.second) {
            weight += y_values[edge_id];
        }
        if (weight <= options.support_tolerance) {
            continue;
        }
        const std::size_t i = entry.first / size;
        const std::size_t j = entry.first % size;
        degree[i] += weight;
        degree[j] += weight;
        weights.push_back({{i, j}, weight});
    }

    std::vector<int> active_of_compact(size, -1);
    std::vector<std::size_t> compact_of_active;
    for (std::size_t i = 0; i < size; ++i) {
        if (degree[i] > options.support_tolerance) {
            active_of_compact[i] = static_cast<int>(compact_of_active.size());
            compact_of_active.push_back(i);
        }
    }
    const std::size_t active_count = compact_of_active.size();
    if (active_count < 3U) {
        return candidates;
    }

    // Padberg--Rao graph: active nodes plus a root that absorbs the degree slack.
    const std::size_t m = active_count + 1U;
    const int root = static_cast<int>(active_count);
    std::vector<std::vector<double>> capacity(m, std::vector<double>(m, 0.0));
    for (const auto& [pair, weight] : weights) {
        const int a = active_of_compact[pair.first];
        const int b = active_of_compact[pair.second];
        capacity[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] += weight;
        capacity[static_cast<std::size_t>(b)][static_cast<std::size_t>(a)] += weight;
    }
    for (std::size_t a = 0; a < active_count; ++a) {
        const double slack = std::max(0.0, 1.0 - degree[compact_of_active[a]]);
        capacity[a][static_cast<std::size_t>(root)] = slack;
        capacity[static_cast<std::size_t>(root)][a] = slack;
    }

    // Root as node 0 of the cut tree keeps the parent pointers simple.
    // Reindex: tree node 0 = root, tree node k (k >= 1) = active node k - 1.
    std::vector<std::vector<double>> tree_capacity(m, std::vector<double>(m, 0.0));
    for (std::size_t a = 0; a < m; ++a) {
        for (std::size_t b = 0; b < m; ++b) {
            const std::size_t ta = a == static_cast<std::size_t>(root) ? 0U : a + 1U;
            const std::size_t tb = b == static_cast<std::size_t>(root) ? 0U : b + 1U;
            tree_capacity[ta][tb] = capacity[a][b];
        }
    }
    const CutTree tree = gusfield_cut_tree(tree_capacity, stats.max_flow_calls);

    std::set<std::vector<int>> seen;
    for (int s = 1; s < static_cast<int>(m); ++s) {
        const double weight = tree.weight[static_cast<std::size_t>(s)];
        if (weight >= 1.0 - options.min_violation) {
            continue;
        }
        const std::vector<bool> side = tree_side(tree, s);
        const bool root_on_side = side[0];
        std::vector<int> compact_set;
        for (std::size_t k = 1; k < m; ++k) {
            const bool in_shore = side[k] != root_on_side;
            if (in_shore) {
                compact_set.push_back(static_cast<int>(compact_of_active[k - 1U]));
            }
        }
        if (compact_set.size() < 3U || compact_set.size() % 2U == 0U) {
            continue;
        }
        std::sort(compact_set.begin(), compact_set.end());
        if (!seen.insert(compact_set).second) {
            continue;
        }
        const double internal = capacity_internal_flow(family, compact_set, y_values);
        const double violation =
            internal - static_cast<double>(compact_set.size() - 1U) / 2.0;
        if (violation <= options.min_violation) {
            continue;
        }
        candidates.push_back(FamilyCandidate{std::move(compact_set), violation});
    }
    return candidates;
}

double jaccard(const std::vector<int>& a, const std::vector<int>& b) {
    std::size_t common = 0;
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] == b[j]) {
            ++common;
            ++i;
            ++j;
        } else if (a[i] < b[j]) {
            ++i;
        } else {
            ++j;
        }
    }
    const std::size_t union_size = a.size() + b.size() - common;
    return union_size == 0U ? 0.0 : static_cast<double>(common) / static_cast<double>(union_size);
}

PstepValidInequalityRow build_row(
    const CapacityBlossomSeparator& separator,
    const CapacityFamilyData& family,
    CapacityFamily family_kind,
    const std::vector<int>& compact_set,
    CapacityBlossomRowForm form
) {
    PstepValidInequalityRow row;
    const std::size_t size = family.nodes.size();
    const std::size_t k = compact_set.size();
    row.name = std::string("capblossom_") + capacity_family_name(family_kind) + "_" +
        std::to_string(k);
    if (form == CapacityBlossomRowForm::Internal) {
        row.sense = PstepValidInequalitySense::LessEqual;
        row.rhs = static_cast<double>(k / 2U);
        for (std::size_t a = 0; a < k; ++a) {
            for (std::size_t b = a + 1U; b < k; ++b) {
                const auto found = family.pair_edges.find(pair_key(
                    static_cast<std::size_t>(compact_set[a]),
                    static_cast<std::size_t>(compact_set[b]),
                    size
                ));
                if (found == family.pair_edges.end()) {
                    continue;
                }
                for (std::size_t edge_id : found->second) {
                    row.edge_terms.emplace_back(edge_id, 1.0);
                }
            }
        }
        return row;
    }
    row.sense = PstepValidInequalitySense::GreaterEqual;
    row.rhs = static_cast<double>((k + 1U) / 2U);
    std::vector<bool> in_set(size, false);
    for (int c : compact_set) {
        in_set[static_cast<std::size_t>(c)] = true;
    }
    for (int c : compact_set) {
        for (std::size_t edge_id : family.in_edges_by_compact[static_cast<std::size_t>(c)]) {
            const int tail = family.compact_index_by_node[static_cast<std::size_t>(separator.edge_tail[edge_id])];
            if (tail >= 0 && in_set[static_cast<std::size_t>(tail)]) {
                continue;
            }
            row.edge_terms.emplace_back(edge_id, 1.0);
        }
    }
    return row;
}

}  // namespace

const char* capacity_family_name(CapacityFamily family) {
    switch (family) {
        case CapacityFamily::Pickup: return "pickup";
        case CapacityFamily::Delivery: return "delivery";
    }
    return "unknown";
}

const char* capacity_blossom_scope_name(CapacityBlossomScope scope) {
    switch (scope) {
        case CapacityBlossomScope::RootOnly: return "root-only";
        case CapacityBlossomScope::AdaptiveTree: return "adaptive-tree";
        case CapacityBlossomScope::FullTree: return "full-tree";
    }
    return "unknown";
}

const char* capacity_blossom_row_form_name(CapacityBlossomRowForm form) {
    switch (form) {
        case CapacityBlossomRowForm::Internal: return "internal";
        case CapacityBlossomRowForm::Inbound: return "inbound";
    }
    return "unknown";
}

double capacity_internal_flow(
    const CapacityFamilyData& family,
    const std::vector<int>& compact_set,
    const std::vector<double>& y_values
) {
    const std::size_t size = family.nodes.size();
    double total = 0.0;
    for (std::size_t a = 0; a < compact_set.size(); ++a) {
        for (std::size_t b = a + 1U; b < compact_set.size(); ++b) {
            const auto found = family.pair_edges.find(pair_key(
                static_cast<std::size_t>(compact_set[a]),
                static_cast<std::size_t>(compact_set[b]),
                size
            ));
            if (found == family.pair_edges.end()) {
                continue;
            }
            for (std::size_t edge_id : found->second) {
                total += y_values[edge_id];
            }
        }
    }
    return total;
}

CapacityBlossomSeparator build_capacity_blossom_separator(const MultiDiGraph& graph) {
    CapacityBlossomSeparator separator;
    separator.edge_count = graph.number_of_edges();
    separator.edge_tail.resize(separator.edge_count);
    separator.edge_head.resize(separator.edge_count);
    for (std::size_t edge_id = 0; edge_id < separator.edge_count; ++edge_id) {
        separator.edge_tail[edge_id] = graph.edges()[edge_id].u;
        separator.edge_head[edge_id] = graph.edges()[edge_id].v;
    }
    build_family(graph, NodeSpec::Kind::Pickup, separator.pickup);
    build_family(graph, NodeSpec::Kind::Delivery, separator.delivery);
    return separator;
}

std::vector<CapacityBlossomCut> separate_capacity_blossom_cuts(
    const CapacityBlossomSeparator& separator,
    const std::vector<double>& y_values,
    const CapacityBlossomOptions& options,
    std::size_t max_cuts,
    CapacityBlossomPool& pool,
    CapacityBlossomStats& stats
) {
    if (y_values.size() != separator.edge_count) {
        throw std::runtime_error(
            "Capacity blossom separation received a y vector of the wrong size."
        );
    }
    const auto start = std::chrono::steady_clock::now();
    std::vector<CapacityBlossomCut> cuts;
    if (max_cuts == 0U) {
        return cuts;
    }

    struct Tagged {
        CapacityFamily family;
        FamilyCandidate candidate;
    };
    std::vector<Tagged> tagged;
    for (FamilyCandidate& candidate :
         separate_family(separator.pickup, y_values, options, stats)) {
        tagged.push_back(Tagged{CapacityFamily::Pickup, std::move(candidate)});
    }
    for (FamilyCandidate& candidate :
         separate_family(separator.delivery, y_values, options, stats)) {
        tagged.push_back(Tagged{CapacityFamily::Delivery, std::move(candidate)});
    }
    stats.candidates += tagged.size();
    std::sort(tagged.begin(), tagged.end(), [](const Tagged& lhs, const Tagged& rhs) {
        if (lhs.candidate.violation != rhs.candidate.violation) {
            return lhs.candidate.violation > rhs.candidate.violation;
        }
        return lhs.candidate.compact_set.size() < rhs.candidate.compact_set.size();
    });

    std::vector<const Tagged*> accepted;
    for (const Tagged& item : tagged) {
        if (cuts.size() >= max_cuts) {
            break;
        }
        std::set<std::vector<int>>& emitted = item.family == CapacityFamily::Pickup
            ? pool.emitted_pickup_sets
            : pool.emitted_delivery_sets;
        if (emitted.count(item.candidate.compact_set) != 0U) {
            ++stats.duplicate_rejections;
            continue;
        }
        bool overlaps = false;
        for (const Tagged* other : accepted) {
            if (other->family == item.family &&
                jaccard(other->candidate.compact_set, item.candidate.compact_set) >
                    options.max_overlap_jaccard) {
                overlaps = true;
                break;
            }
        }
        if (overlaps) {
            ++stats.overlap_rejections;
            continue;
        }
        accepted.push_back(&item);
        emitted.insert(item.candidate.compact_set);

        const CapacityFamilyData& family = item.family == CapacityFamily::Pickup
            ? separator.pickup
            : separator.delivery;
        CapacityBlossomCut cut;
        cut.family = item.family;
        cut.compact_set = item.candidate.compact_set;
        cut.violation = item.candidate.violation;
        cut.row = build_row(separator, family, item.family, cut.compact_set, options.row_form);
        stats.largest_set_size = std::max(stats.largest_set_size, cut.compact_set.size());
        if (item.family == CapacityFamily::Pickup) {
            ++stats.pickup_cuts;
        } else {
            ++stats.delivery_cuts;
        }
        cuts.push_back(std::move(cut));
    }
    stats.separation_seconds +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return cuts;
}

double max_capacity_blossom_violation_by_enumeration(
    const CapacityBlossomSeparator& separator,
    CapacityFamily family_kind,
    const std::vector<double>& y_values,
    std::size_t max_family_size
) {
    const CapacityFamilyData& family = family_kind == CapacityFamily::Pickup
        ? separator.pickup
        : separator.delivery;
    const std::size_t size = family.nodes.size();
    if (size > max_family_size || size >= 63U) {
        throw std::runtime_error("Capacity blossom enumeration family is too large.");
    }
    double best = -std::numeric_limits<double>::infinity();
    const std::uint64_t limit = std::uint64_t{1} << size;
    std::vector<int> compact_set;
    for (std::uint64_t mask = 1; mask < limit; ++mask) {
        const int popcount = __builtin_popcountll(mask);
        if (popcount < 3 || popcount % 2 == 0) {
            continue;
        }
        compact_set.clear();
        for (std::size_t i = 0; i < size; ++i) {
            if ((mask >> i) & 1ULL) {
                compact_set.push_back(static_cast<int>(i));
            }
        }
        const double violation = capacity_internal_flow(family, compact_set, y_values) -
            static_cast<double>(popcount - 1) / 2.0;
        best = std::max(best, violation);
    }
    return best;
}

}  // namespace spdp

