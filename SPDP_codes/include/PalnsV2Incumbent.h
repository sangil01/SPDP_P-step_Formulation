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
    // (1) Phase II trip selection (composite choice and pickup-delivery
    //     assignment) by cost instead of duration.
    bool sel_cost = false;
    // (2) Focus cycles in Phase II. Every focus_period-th simulated-annealing
    //     cycle is a focus cycle (0 disables; 2 = every other cycle). A focus
    //     cycle is a sequence of episodes of focus_len iterations. An episode
    //     fixes all routes except focus_routes of them (a seed route, taken in
    //     turn, and routes close to it), runs its own annealing on those routes
    //     starting from the thread's best solution, and puts the best feasible
    //     sub-solution it found back into the solution.
    int focus_period = 0;
    int focus_routes = 3;
    int focus_len = 4000;
    double focus_kmax_ratio = 0.3;    // k <= max(kmax_min, ceil(ratio * trajectories of the focus routes))
    double focus_rebuild = 0.1;       // probability of removing every trajectory of the focus routes
    int focus_deep = 2;               // trajectory relocate in an episode: 0 never, 1 always, 2 near the episode best
    double focus_start_ratio = -1.0;  // annealing start ratio of an episode (negative: sa_start_ratio)
    // (3) Trajectory relocate (the deep local search) evaluated with a cache of
    //     best insertions per (trajectory, unchanged route). Same moves in the
    //     same order, so the search itself does not change; only time. On by
    //     default; fast_relocate=0 uses the routine of version 1.
    bool fast_relocate = true;
    // (4) Portfolio diversity. temp_spread = f > 1: thread i uses the annealing
    //     start ratio times {1, f, 1/f}[i mod 3]. sel_cost_threads = m > 0:
    //     every m-th thread selects by cost (as sel_cost), the others by duration.
    double temp_spread = 0.0;
    int sel_cost_threads = 0;
    // (5) Phase II penalty weight (lambda) control. lambda0 is the start value
    //     max(1, 2C/D). With the defaults the control is that of version 1:
    //     every 100 iterations, times 1.5 if fewer than 40% of the candidates
    //     were feasible, divided by 1.5 if more than 60%, within [1e-3, 1e6].
    double lam_lo = 0.0;    // > 0: lower bound lam_lo * lambda0 (0: 1e-3)
    double lam_hi = 0.0;    // > 0: upper bound lam_hi * lambda0 (0: 1e6)
    double lam_step = 1.5;  // factor of one update
    double lam_dstep = 0.0; // > 1: factor of an update downwards (0: lam_step)
    int lam_win = 100;      // iterations per update
    int lam_sig = 0;        // share measured over the window: 0 feasible candidates,
                            // 1 iterations whose current solution is feasible
    int lam_back = 0;       // > 0: the current solution has been infeasible for this
                            // many iterations in a row -> go back to a feasible solution
    int lam_to = 0;         // where lam_back goes: 0 thread best, 1 the last feasible
                            // current solution (the one the excursion started from)
    double lam_tlo = 0.4;   // share below which lambda goes up
    double lam_thi = 0.6;   // share above which lambda goes down
    // (6) Threads that have no feasible solution of their own.
    //     p1_join >= 0: an independent (odd) thread still in Phase I takes the
    //     shared best and starts Phase II once a shared feasible solution has
    //     existed for p1_join iterations (-1: it stays in Phase I, version 1).
    //     rescue_p2: after the rescue found a feasible solution every thread
    //     runs Phase II from it (version 1 stops there).
    int p1_join = -1;
    bool rescue_p2 = false;
    //     indep_adopt = m > 0: an independent (odd) thread in Phase II whose own
    //     best did not improve during m annealing cycles in a row takes the
    //     shared best at the end of the cycle, if that is better (0: never).
    int indep_adopt = 0;
    // (7) Adaptive annealing temperature in Phase II. The start temperature of a
    //     cycle is multiplied by m (1 at the beginning). At the end of a cycle:
    //     m = 1 if the thread best improved in it, otherwise after reheat_after
    //     such cycles m = min(reheat_max, m + reheat_step). reheat_max <= 1: off.
    double reheat_max = 1.0;
    double reheat_step = 1.0;
    int reheat_after = 1;
    // Diagnostics: time per part of an iteration in the log (does not change the search).
    bool profile = false;
    // Diagnostics: Phase II counters in the log, per thread (operators, removal
    // size, trips inserted, local-search moves, one line per annealing cycle).
    // Does not change the search.
    bool diag = false;
    // Diagnostics: block pool (does not change the search). pool_ratio = r > 0:
    // every feasible Phase II candidate whose travel cost is at most (1 + r)
    // times the thread best is cut into blocks (the part of a route between two
    // moments with an empty vehicle). pool_file receives, per thread, every
    // block with the iteration at which it was first seen and the iterations at
    // which the cheapest solution containing it got cheaper, together with the
    // request classes, the shared-best trace and the iteration clock.
    double pool_ratio = 0.0;
    std::string pool_file;
    // Block pool kept in memory for the recombination step (does not change the
    // search). pool_delta = d > 0: a feasible Phase II candidate is cut into
    // blocks when its travel cost is at most the thread best + d. With both
    // pool_ratio and pool_delta set, the larger of the two limits applies.
    double pool_delta = 0.0;
    // Recombination after PALNS (PoolRecombination.h), run by the caller of
    // solve_fixed_k. recomb_tl = seconds for the restricted model (0: off);
    // recomb_delta: blocks whose cheapest solution costs at most the best cost
    // + recomb_delta are used; recomb_threads: Gurobi threads (0: PALNS threads).
    // recomb_tl > 0 switches the pool on with pool_delta >= recomb_delta.
    double recomb_tl = 0.0;
    double recomb_delta = 3.0;
    int recomb_threads = 0;
    // Text of the configuration these options were parsed from (for the log).
    std::string config_text;
};

// Parses "key=value,key=value" into the version 2 switches. Unknown keys and
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
