#ifndef SPDP_CAPACITY_BLOSSOM_CUTS_H
#define SPDP_CAPACITY_BLOSSOM_CUTS_H

#include <cstddef>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "GenMultiGraph.h"
#include "PstepValidInequality.h"

namespace spdp {

// Capacity blossom cuts.
//
// A vehicle holds at most Q = 2 skips, a pickup occupies a slot, emptying keeps
// the slot occupied (F -> E) and only a delivery frees it. Hence any route
// visits at most two pickups consecutively and at most two deliveries
// consecutively, whatever treatment happens on the connecting arcs. For every
// pickup set S (or delivery set S) the integer solutions satisfy
//
//     sum_{e in delta^-(S)} y_e >= ceil(|S| / 2)          (inbound form)
//     sum_{e in E(S)}       y_e <= floor(|S| / 2)         (internal form)
//
// which are equivalent under the unit in/out-degree rows. Aggregating parallel
// edges between two pickups into z_ij turns the internal form into the odd-set
// (blossom) inequality of the fractional matching polytope, so the exact
// separation is a minimum odd cut (Padberg--Rao) on the support graph.

enum class CapacityFamily {
    Pickup,
    Delivery,
};

enum class CapacityBlossomScope {
    RootOnly,
    AdaptiveTree,
    FullTree,
};

enum class CapacityBlossomRowForm {
    Internal,  // sum of edges inside S <= floor(|S|/2)   (sparse for small S)
    Inbound,   // sum of edges entering S >= ceil(|S|/2)
};

struct CapacityBlossomOptions {
    bool enabled = false;
    CapacityBlossomScope scope = CapacityBlossomScope::AdaptiveTree;
    CapacityBlossomRowForm row_form = CapacityBlossomRowForm::Internal;
    std::size_t root_max_rounds = 20;
    std::size_t root_max_cuts = 256;
    std::size_t root_max_per_round = 32;
    std::size_t tree_max_per_round = 8;
    // Adaptive tree: separate at every node while node_count < tree_dense_node_limit,
    // afterwards only when node_count % tree_node_frequency == 0.
    std::size_t tree_dense_node_limit = 50;
    std::size_t tree_node_frequency = 100;
    std::size_t max_total_cuts = 2000;
    double min_violation = 1e-4;
    double support_tolerance = 1e-6;
    // Stop tree separation when the accumulated separation time exceeds this
    // fraction of the solver runtime (0 disables the check).
    double max_separation_time_fraction = 0.15;
    // Candidates of one round whose node set overlaps an already accepted
    // candidate (same family) by more than this Jaccard index are skipped.
    double max_overlap_jaccard = 0.9;
    // Root rounds stop when the best violation of two consecutive rounds is
    // below this value (tailing off on violation, not on the bound).
    double root_tailing_off_violation = 0.02;
};

struct CapacityBlossomStats {
    std::size_t rounds = 0;
    std::size_t root_rounds = 0;
    std::size_t candidates = 0;
    std::size_t cuts_added = 0;
    std::size_t root_cuts_added = 0;
    std::size_t pickup_cuts = 0;
    std::size_t delivery_cuts = 0;
    std::size_t duplicate_rejections = 0;
    std::size_t overlap_rejections = 0;
    std::size_t largest_set_size = 0;
    std::size_t max_flow_calls = 0;
    bool stopped_by_total_cap = false;
    bool stopped_by_time_fraction = false;
    double separation_seconds = 0.0;
};

struct CapacityBlossomCut {
    CapacityFamily family = CapacityFamily::Pickup;
    std::vector<int> compact_set;  // sorted compact indices inside the family
    double violation = 0.0;
    PstepValidInequalityRow row;
};

struct CapacityFamilyData {
    std::vector<NodeId> nodes;               // compact index -> node id
    std::vector<int> compact_index_by_node;  // node id -> compact index or -1
    // Edges between two family nodes (both directions, all parallel copies),
    // keyed by min(i,j) * size + max(i,j).
    std::unordered_map<std::size_t, std::vector<std::size_t>> pair_edges;
    // Incoming edge ids for every family node (all sources).
    std::vector<std::vector<std::size_t>> in_edges_by_compact;
};

struct CapacityBlossomSeparator {
    std::size_t edge_count = 0;
    std::vector<NodeId> edge_tail;
    std::vector<NodeId> edge_head;
    CapacityFamilyData pickup;
    CapacityFamilyData delivery;
};

// Sets already handed to the solver; the solver cannot remove user cuts, so the
// separator never emits the same set twice.
struct CapacityBlossomPool {
    std::set<std::vector<int>> emitted_pickup_sets;
    std::set<std::vector<int>> emitted_delivery_sets;
};

CapacityBlossomSeparator build_capacity_blossom_separator(const MultiDiGraph& graph);

// Exact odd-set separation for both families on the fractional point y.
// Returns at most max_cuts cuts sorted by decreasing violation; every returned
// cut is re-evaluated on the original y vector and is new with respect to the
// pool (the pool is updated with the returned sets).
std::vector<CapacityBlossomCut> separate_capacity_blossom_cuts(
    const CapacityBlossomSeparator& separator,
    const std::vector<double>& y_values,
    const CapacityBlossomOptions& options,
    std::size_t max_cuts,
    CapacityBlossomPool& pool,
    CapacityBlossomStats& stats
);

// Exact reference used by tests: enumerates every odd subset of one family
// (feasible only for small families) and returns the maximum violation.
double max_capacity_blossom_violation_by_enumeration(
    const CapacityBlossomSeparator& separator,
    CapacityFamily family,
    const std::vector<double>& y_values,
    std::size_t max_family_size
);

// Evaluates the internal-form left-hand side sum_{e in E(S)} y_e for a set S.
double capacity_internal_flow(
    const CapacityFamilyData& family,
    const std::vector<int>& compact_set,
    const std::vector<double>& y_values
);

const char* capacity_family_name(CapacityFamily family);
const char* capacity_blossom_scope_name(CapacityBlossomScope scope);
const char* capacity_blossom_row_form_name(CapacityBlossomRowForm form);

}  // namespace spdp

#endif

