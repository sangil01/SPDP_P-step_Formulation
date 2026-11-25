#ifndef SPDP_ARC_ELIMINATION_H
#define SPDP_ARC_ELIMINATION_H

#include <array>
#include <ostream>
#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "ReadData.h"
#include "TwoIndexSolver.h"

namespace spdp {

// Gurobi logs of the P=0 direct solve, "<instance>_p0_<name>_gurobi.log",
// named by the edge set and the order of the solves on it: E0 the full edge
// set, E1 the edges left by the root LP deletion (Step 1 below), E2 the edges
// collected in Step 2. Without arc elimination only E0_1 and E0_2 are written.
inline constexpr const char* kLogE0RootLp = "E0_1_rootlp";            // Step 1, or the pre-MIP root LP without arc elimination
inline constexpr const char* kLogE0Milp = "E0_2_milp";                // main MIP without arc elimination
inline constexpr const char* kLogE1InfeasibleCheck = "E1_1_infeasible_check";       // Step 2a
inline constexpr const char* kLogE1FeasibleArcCollect = "E1_2_feasible_arc_collect"; // Step 2b
inline constexpr const char* kLogE1RootLp = "E1_3_rootlp";            // fallback: pre-MIP root LP on E1
inline constexpr const char* kLogE1Milp = "E1_4_milp";                // fallback: MIP on E1
inline constexpr const char* kLogE2RootLp = "E2_1_rootlp";            // Step 3: pre-MIP root LP on E2
inline constexpr const char* kLogE2Milp = "E2_2_milp";                // Step 3: MIP on E2
inline constexpr std::array<const char*, 8> kP0GurobiLogNames = {
    kLogE0RootLp, kLogE0Milp, kLogE1InfeasibleCheck, kLogE1FeasibleArcCollect,
    kLogE1RootLp, kLogE1Milp, kLogE2RootLp, kLogE2Milp};

// prefix + name + "_gurobi.log"; empty prefix: no log (empty path).
inline std::string p0_gurobi_log_path(const std::string& prefix, const char* name) {
    return prefix.empty() ? std::string() : prefix + name + "_gurobi.log";
}

// Arc elimination with the time-free model (P=0 direct two-index IP, original
// cost objective, integer edge costs). U is the cost of the known solution; the
// aim is an edge set E2 that contains every edge of every solution of cost
// <= U - 1, so that the optimum is min(U, optimum of the main model on E2).
//
//   Step 1  Root LP of the main model (with its cuts). An edge e with
//           (lower bound of the LP with y_e = 1) > U - 1 is deleted: E -> E1.
//   Step 2  Merged time-free model on E1: the main model without its time rows,
//           with equal requests merged (a relaxation of the main model).
//     2a    With "cost <= U - 1", stopped at its first solution. No solution:
//           U is optimal. A solution: its classes join E2, which starts with
//           the classes of the incumbent edges.
//     2b    One branch-and-bound with "cost <= U - 1" and "at least one class
//           outside E2": every integer solution found puts its classes into E2
//           and is rejected by a lazy constraint. When the search ends without
//           a solution, no time-free solution of cost <= U - 1 uses an edge
//           outside E2, so E2 is complete.
//   Step 3  Main model on E2 with "cost <= U - 1". No solution: U is optimal.
//           A solution: it is optimal.
//   Fallback  If Step 2 does not finish, the main model is solved on E1, with the
//           known solution as MIP start when all its edges are in E1.
//
// A class is the set of edges that are the same move for different requests:
// same kinds, locations, head container type, vehicle states and facility
// sequence. Adding more edges to E2 never hurts the argument.
struct ArcEliminationOptions {
    bool enabled = false;
    // Step 1. The root LP always gives the bound check (LP bound >= U: U is
    // optimal). Edge deletion: 0 none (E1 = E), 1 barrier with crossover,
    // 2 barrier without crossover (0 also runs without crossover).
    // 1 and 2 use a bound computed from the dual values that is valid for any
    // dual vector, so the choice only changes the strength and the time.
    int rc_mode = 2;
    // Seconds for the root LP of Step 1 (0: no limit of its own). It took at
    // most 49 s in the 67-instance run with rc_mode = 2.
    double rc_time_limit = 300.0;
    // > 0: everything before arc elimination (PALNS included) plus Steps 1 and 2
    // must end within this many seconds of the instance start; Step 2 is cut at
    // that point, and Steps 1 and 2 are skipped when less than two minutes are
    // left. PALNS itself is never cut by this limit. 0: no such limit (Step 2
    // may then use the whole budget). Default 600 s: in the 67-instance run
    // (2026-11-17, limit 1200 s) every instance whose Step 2 finished did so
    // within 428 s of the instance start (D6).
    double total_time_limit = 600.0;
    // Seconds of the instance used before arc elimination starts (set by the caller).
    double elapsed_before = 0.0;
    // "<output dir>/<instance>_p0_" (set by the caller): the Gurobi logs of the
    // steps are p0_gurobi_log_path(log_prefix, kLog...). Empty: no logs.
    std::string log_prefix;
};

struct ArcEliminationResult {
    bool ran = false;
    // result holds the final answer and the main solve is not needed.
    bool finished = false;
    bool proved_optimal = false;
    double incumbent_cost = -1.0;  // cost of the best solution at the end
    double lower_bound = -1.0;     // valid lower bound on the optimal cost (-1: none)
    double seconds = 0.0;
    DirectTwoIndexResult result;
};

// Keys: rc, rc_tl, total_tl. Returns false and fills error for an unknown key
// or a bad value.
bool parse_arc_elimination_config(
    const std::string& text,
    ArcEliminationOptions& options,
    std::string& error
);

// main_options are the options of the main model (temporal rows, valid
// inequalities, cut families, Gurobi settings); its solver_time_limit is the
// budget of the whole procedure. incumbent_edge_values is the 0/1 edge vector
// of the known solution and incumbent_cost its original cost.
ArcEliminationResult run_arc_elimination(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const DirectTwoIndexOptions& main_options,
    const std::vector<double>& incumbent_edge_values,
    double incumbent_cost,
    const ArcEliminationOptions& options,
    std::ostream* log
);

}  // namespace spdp

#endif
