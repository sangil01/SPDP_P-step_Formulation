#ifndef SPDP_PALNS_V2_INCUMBENT_H
#define SPDP_PALNS_V2_INCUMBENT_H

#include <memory>
#include <string>

#include "PalnsIncumbent.h"
#include "PoolRecombination.h"

namespace spdp {

// Version 2 of the PALNS incumbent heuristic (stage 2 of
// SPDP_test_research/05_arc_elimination_palns/notes/00_plan.md).
// It lives next to PalnsIncumbentSolver so that both can be run under the same
// conditions; PalnsIncumbent.cpp is not touched. The source started as a copy
// of version 1: with every switch below at its default the search is identical
// to version 1 for the same PalnsOptions (same iterations, same solutions).
// Only the time differs: the insertion evaluation, the regret/greedy insertion
// and the trajectory relocate are computed in a faster but equivalent way.
struct PalnsV2Options {
    PalnsOptions base;
    // Phase II penalty weight (lambda) control. The update itself is that of
    // version 1: every 100 iterations, times 1.5 if fewer than 40% of the
    // candidates were feasible, divided by 1.5 if more than 60%, within
    // [1e-3, 1e6].
    int lam_back = 0;       // > 0: the current solution has been infeasible for this
                            // many iterations in a row -> go back to a feasible solution
    int lam_to = 0;         // where lam_back goes: 0 thread best, 1 the last feasible
                            // current solution (the one the excursion started from)
    // Threads that have no feasible solution of their own.
    //     p1_join >= 0: an independent (odd) thread still in Phase I takes the
    //     shared best and starts Phase II once a shared feasible solution has
    //     existed for p1_join iterations (-1: it stays in Phase I, version 1).
    //     rescue_p2: after the rescue found a feasible solution every thread
    //     runs Phase II from it (version 1 stops there).
    int p1_join = -1;
    bool rescue_p2 = false;
    // Recombination after PALNS (PoolRecombination.h), run by the caller of
    // solve_fixed_k. recomb_tl = seconds for the restricted model (0: off);
    // recomb_delta: blocks whose cheapest solution costs at most the best cost
    // + recomb_delta are used. The model uses as many Gurobi threads as PALNS.
    double recomb_tl = 0.0;
    double recomb_delta = 3.0;
    // Block pool kept in memory for the recombination step (does not change the
    // search): a feasible Phase II candidate is cut into blocks (the part of a
    // route between two moments with an empty vehicle) when its travel cost is
    // at most the thread best + pool_delta. Set to recomb_delta by the
    // constructor when recomb_tl > 0; not a configuration key.
    double pool_delta = 0.0;
    // Text of the configuration these options were parsed from (for the log).
    std::string config_text;
};

// Parses "key=value,key=value" into the version 2 switches. Keys: lam_back,
// lam_to, p1_join, rescue_p2, recomb_tl, recomb_delta. Unknown keys and
// malformed values are errors. An empty text leaves every switch off.
bool parse_palns_v2_config(const std::string& text, PalnsV2Options& options, std::string& error);

class PalnsV2IncumbentSolver {
public:
    PalnsV2IncumbentSolver(const SPDPData& data, const PalnsV2Options& options);
    ~PalnsV2IncumbentSolver();
    PalnsV2IncumbentSolver(const PalnsV2IncumbentSolver&) = delete;
    PalnsV2IncumbentSolver& operator=(const PalnsV2IncumbentSolver&) = delete;

    double pattern_generation_seconds() const;
    std::string pattern_summary() const;
    PalnsFixedKResult solve_fixed_k(int vehicle_count, double time_limit_seconds);
    // Blocks recorded by the last solve_fixed_k call, merged over the threads
    // (empty unless the pool is on).
    const std::vector<PoolBlock>& pool_blocks() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace spdp

#endif
