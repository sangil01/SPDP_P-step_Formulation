#ifndef SPDP_INITIAL_INCUMBENT_H
#define SPDP_INITIAL_INCUMBENT_H

#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "ReadData.h"

namespace spdp {

struct InitialIncumbentSolveOptions {
    int vehicle_count = 0;
    double solver_time_limit = 0.0;  // 0 means no time limit.
    int gurobi_threads = -1;
    std::string gurobi_log_path;
};

struct InitialIncumbentSolveResult {
    int status = 0;
    bool hit_time_limit = false;
    bool hit_solution_limit = false;
    bool infeasible = false;
    bool has_feasible_solution = false;
    double runtime_seconds = 0.0;
    double total_duration = -1.0;
    double total_original_cost = -1.0;
    int vehicle_count = 0;
    std::vector<double> edge_values;
};

enum class InitialIncumbentTimeoutAction {
    Stop,
    Advance,
};

struct InitialIncumbentSearchOptions {
    int initial_vehicle_count = 0;
    int max_k_increments = 0;
    InitialIncumbentTimeoutAction timeout_action =
        InitialIncumbentTimeoutAction::Stop;
    double per_attempt_time_limit = 0.0;
    int gurobi_threads = -1;
    std::string gurobi_log_base_path;
};

struct InitialIncumbentAttemptResult {
    int attempt_index = 0;
    int vehicle_count = 0;
    std::string gurobi_log_path;
    InitialIncumbentSolveResult solve_result;
};

struct InitialIncumbentSearchResult {
    int initial_k = 0;
    int certified_k = 0;
    int last_candidate_k = 0;
    int incumbent_k = 0;
    bool has_feasible_solution = false;
    bool stopped_on_time_limit = false;
    bool exhausted_increment_limit = false;
    bool exhausted_candidate_limit = false;
    double total_runtime_seconds = 0.0;
    std::vector<InitialIncumbentAttemptResult> attempts;
    InitialIncumbentSolveResult incumbent_result;
};

InitialIncumbentSolveResult solve_fixed_k_duration_initial_incumbent(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const InitialIncumbentSolveOptions& options
);

InitialIncumbentSearchResult solve_iterative_duration_initial_incumbent(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const InitialIncumbentSearchOptions& options
);

}  // namespace spdp

#endif
