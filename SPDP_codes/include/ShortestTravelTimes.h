#ifndef SPDP_SHORTEST_TRAVEL_TIMES_H
#define SPDP_SHORTEST_TRAVEL_TIMES_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include "ReadData.h"

namespace spdp {

// All-pairs shortest travel times on a directed, nonnegative travel-time
// matrix. Computed with one Dijkstra search per source. Independent of the
// SPDP request structure so that every consumer (metric closure for restricted
// instances, duration labels, diagnostics) uses the same numbers.
struct ShortestTravelTimes {
    std::size_t location_count = 0;
    // distance[a][b] = length of the shortest directed path a -> b
    // (0 on the diagonal, +infinity when b is unreachable from a).
    std::vector<std::vector<double>> distance;
    // predecessor[a][b] = node preceding b on the shortest path from a
    // (-1 for b == a or when unreachable).
    std::vector<std::vector<int>> predecessor;
    // Number of ordered pairs (a, b), a != b, whose direct travel time is
    // strictly longer than the shortest path (non-metric pairs).
    std::size_t non_metric_pair_count = 0;
    // Largest reduction raw[a][b] - distance[a][b] over all pairs.
    double max_reduction = 0.0;

    double between(int from, int to) const;
    bool reachable(int from, int to) const;
    // Node sequence from 'from' to 'to', both endpoints included. Empty when
    // unreachable; {from} when from == to.
    std::vector<int> path(int from, int to) const;
    std::uint64_t fingerprint() const;
};

// Directed Dijkstra from every source. The matrix must be square with
// nonnegative finite entries; the diagonal is ignored.
ShortestTravelTimes compute_shortest_travel_times(
    const std::vector<std::vector<double>>& travel_time
);

// Shortest travel times between the physical locations 0..data.locations-1 of
// an instance. The virtual depot row/column (index data.locations) is excluded,
// because its zero entries would otherwise connect every pair at cost zero.
ShortestTravelTimes compute_physical_shortest_travel_times(const SPDPData& data);

// Returns a copy of data.time whose physical block is replaced by the shortest
// travel times. The virtual depot row/column keeps its original (zero) values,
// so open routes still start and end without depot travel.
std::vector<std::vector<double>> metric_closure_travel_time_matrix(
    const SPDPData& data,
    const ShortestTravelTimes& shortest
);

}  // namespace spdp

#endif

