#ifndef SPDP_CONNECTIVITY_CUTS_H
#define SPDP_CONNECTIVITY_CUTS_H

#include <cstddef>
#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "PstepValidInequality.h"
#include "ReadData.h"

namespace spdp {

enum class ConnectivityCutRhsMode {
    One,       // classic connectivity: every service-node set is entered at least once
    Duration,  // max(1, ceil(H_min(S) / T)) with H_min(S) = sum of minimum in-edge durations
};

struct ConnectivityCutOptions {
    bool enabled = false;
    bool root_only = false;
    ConnectivityCutRhsMode rhs_mode = ConnectivityCutRhsMode::Duration;
    std::size_t max_cuts_per_round = 32;
    double violation_tolerance = 1e-4;
    double support_tolerance = 1e-6;
};

struct ConnectivityCutStats {
    std::size_t rounds = 0;
    std::size_t cuts_added = 0;
    std::size_t root_cuts_added = 0;
    std::size_t rhs_two_or_more_cuts = 0;
    double separation_seconds = 0.0;
};

// Precomputed per-node data used by the separator.
struct ConnectivitySeparator {
    NodeId end_node_id = 0;
    std::size_t node_count = 0;
    double route_time_limit = 0.0;
    std::vector<double> min_in_edge_time_by_node;
    // Incoming graph edge ids per node (all parallel edges).
    std::vector<std::vector<std::size_t>> in_edge_ids_by_node;
    std::vector<NodeId> edge_tail;
    std::vector<NodeId> edge_head;
};

ConnectivitySeparator build_connectivity_separator(
    const SPDPData& data,
    const MultiDiGraph& graph
);

// Exact separation on the fractional support: for every service node v,
// compute the maximum depot->v flow with capacities y. If it is below the
// required entry count of the sink-side set S, return the cut
//     sum_{e in delta^-(S)} y_e >= rhs(S).
std::vector<PstepValidInequalityRow> separate_connectivity_cuts(
    const ConnectivitySeparator& separator,
    const std::vector<double>& y_values,
    const ConnectivityCutOptions& options
);

const char* connectivity_cut_rhs_mode_name(ConnectivityCutRhsMode mode);

}  // namespace spdp

#endif
