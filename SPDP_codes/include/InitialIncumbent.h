#ifndef SPDP_INITIAL_INCUMBENT_H
#define SPDP_INITIAL_INCUMBENT_H

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "CpSatSpdpSolver.h"
#include "GenMultiGraph.h"
#include "ReadData.h"

namespace spdp {

enum class InitialIncumbentBackend {
    TwoIndexMilp,
    CpSat,
};

enum class InitialIncumbentMilpMode {
    Duration,
    Makespan,
};

const char* initial_incumbent_milp_mode_name(InitialIncumbentMilpMode mode);

enum class FixedKSolveOutcome {
    Feasible,
    ProvenInfeasible,
    TimedOutUnknown,
    EarlyUnknown,
    ModelInvalid,
    AdapterError,
};

struct InitialIncumbentSolveOptions {
    int vehicle_count = 0;
    double solver_time_limit = 0.0;  // 0 means no time limit.
    int gurobi_threads = -1;
    std::string gurobi_log_path;
    InitialIncumbentMilpMode milp_mode = InitialIncumbentMilpMode::Duration;
    double milp_makespan_horizon_factor = 1.5;
};

struct InitialIncumbentSolveResult {
    InitialIncumbentBackend backend = InitialIncumbentBackend::TwoIndexMilp;
    FixedKSolveOutcome outcome = FixedKSolveOutcome::EarlyUnknown;
    int status = 0;
    int raw_status = 0;
    std::string status_name;
    std::string termination_name;
    bool hit_time_limit = false;
    bool early_unknown = false;
    bool hit_solution_limit = false;
    bool infeasible = false;
    bool has_feasible_solution = false;
    double runtime_seconds = 0.0;
    double configured_time_limit_seconds = 0.0;
    InitialIncumbentMilpMode milp_mode = InitialIncumbentMilpMode::Duration;
    double milp_makespan_horizon_factor = 0.0;
    double milp_model_horizon = 0.0;
    bool milp_has_objective_value = false;
    double milp_objective_value = 0.0;
    bool milp_has_objective_bound = false;
    double milp_best_objective_bound = 0.0;
    bool milp_stopped_by_feasible_callback = false;
    bool milp_stopped_by_bound_callback = false;
    CpGraphMode cp_graph_mode = CpGraphMode::OriginalGraph;
    CpSolveMode cp_solve_mode = CpSolveMode::Satisfaction;
    double cp_threshold_horizon_factor = 0.0;
    std::int64_t cp_model_horizon = 0;
    bool cp_has_objective_value = false;
    double cp_objective_value = 0.0;
    bool cp_has_objective_bound = false;
    double cp_best_objective_bound = 0.0;
    bool cp_stopped_by_feasible_observer = false;
    bool cp_stopped_by_bound_callback = false;
    double total_duration = -1.0;
    double total_original_cost = -1.0;
    int vehicle_count = 0;
    std::vector<double> edge_values;
    std::vector<int> active_edge_ids;
    std::int64_t conflicts = 0;
    std::int64_t branches = 0;
    CpModelBuildStats cp_build_stats;
    std::string adapter_status;
    std::string error_message;
    std::string solution_info;
};

enum class InitialIncumbentTimeoutAction {
    Stop,
    Advance,
};

struct InitialIncumbentSearchOptions {
    InitialIncumbentBackend backend = InitialIncumbentBackend::TwoIndexMilp;
    InitialIncumbentMilpMode milp_mode = InitialIncumbentMilpMode::Duration;
    double milp_makespan_horizon_factor = 1.5;
    int initial_vehicle_count = 0;
    int max_k_increments = 0;
    InitialIncumbentTimeoutAction timeout_action =
        InitialIncumbentTimeoutAction::Stop;
    double per_attempt_time_limit = 0.0;
    int gurobi_threads = -1;
    std::string gurobi_log_base_path;
    std::string cp_sat_log_base_path;
    CpSatSolveOptions cp_sat;
};

using FixedKSolveCallback =
    std::function<InitialIncumbentSolveResult(int vehicle_count)>;

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
    bool stopped_on_early_unknown = false;
    int early_unknown_k = 0;
    bool exhausted_increment_limit = false;
    bool exhausted_candidate_limit = false;
    double total_runtime_seconds = 0.0;
    std::vector<InitialIncumbentAttemptResult> attempts;
    InitialIncumbentSolveResult incumbent_result;
};

InitialIncumbentSolveResult solve_fixed_k_initial_incumbent(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const InitialIncumbentSolveOptions& options
);

// Backward-compatible duration-only entry point.
InitialIncumbentSolveResult solve_fixed_k_duration_initial_incumbent(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const InitialIncumbentSolveOptions& options
);

InitialIncumbentSearchResult solve_iterative_initial_incumbent(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const InitialIncumbentSearchOptions& options
);

// Backward-compatible duration-only entry point.
InitialIncumbentSearchResult solve_iterative_duration_initial_incumbent(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const InitialIncumbentSearchOptions& options
);

InitialIncumbentSearchResult run_iterative_initial_incumbent_search(
    const InitialIncumbentSearchOptions& options,
    int max_candidate_k,
    const FixedKSolveCallback& solve_fixed_k
);

}  // namespace spdp

#endif
