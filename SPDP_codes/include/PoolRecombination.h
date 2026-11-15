#ifndef SPDP_POOL_RECOMBINATION_H
#define SPDP_POOL_RECOMBINATION_H

#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "ReadData.h"

namespace spdp {

// Recombination step after PALNS version 2 (plan 4.4 of
// SPDP_test_research/05_arc_elimination_palns/notes/00_plan.md).
//
// During Phase II the search cuts its good feasible candidates into blocks (the
// part of a route between two moments with an empty vehicle) and keeps, per
// block, the cost of the cheapest solution that contained it. After PALNS the
// real model (two-index model on the multigraph) is solved on a subset of the
// edges, starting from the PALNS solution:
//   - every edge of the kinds used inside a block whose cheapest solution costs
//     at most (best cost in the pool) + delta, where a kind is (kind and
//     location of the tail, kind, location and container type of the head,
//     facility sequence in between), so every request of a class is included;
//   - every edge that leaves the depot or a delivery with an empty vehicle
//     (links between blocks, start and end of a route);
//   - the edges of the PALNS solution itself.
// The restricted model contains the PALNS solution, so its result is never worse.

// One operation of a block. kind 'P': pickup at `location`; 'E': emptying at the
// facility `location`; 'D': delivery at `location`.
struct PoolOp {
    char kind = 'P';
    int location = -1;
    int container_type = -1;
};

struct PoolBlock {
    std::vector<PoolOp> ops;
    double min_cost = 0.0;  // original cost of the cheapest solution containing the block
};

struct PoolRecombinationOptions {
    double delta = 3.0;        // blocks with min_cost <= best + delta are used
    double time_limit = 120.0; // seconds for the restricted model
    int threads = 8;
    int mip_focus = 1;
};

struct PoolRecombinationResult {
    bool ran = false;       // the restricted model was built and optimized
    bool improved = false;  // a solution cheaper than the start was found and validated
    std::vector<int> active_edge_ids;  // edges of the improved solution (ids of `graph`)
    double total_duration = 0.0;
    double total_original_cost = 0.0;
    double seconds = 0.0;   // wall time of the whole step
    std::string log;        // one line, "[palns-recomb] key=value ..."
};

// `start_edge_ids`: the PALNS solution mapped to `graph` (ids of graph.edges()).
// The number of vehicles is fixed to `vehicle_count`, as in PALNS.
PoolRecombinationResult recombine_pool_edges(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const std::vector<PoolBlock>& blocks,
    const std::vector<int>& start_edge_ids,
    int vehicle_count,
    const PoolRecombinationOptions& options
);

}  // namespace spdp

#endif
