#ifndef SPDP_CP_SAT_SPDP_SOLVER_H
#define SPDP_CP_SAT_SPDP_SOLVER_H

#include <cstdint>
#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "ReadData.h"

namespace spdp {

enum class CpActionKind { Pickup, Treatment, Delivery };
enum class CpSolveMode {
    Satisfaction,
    ThresholdOptimization,
};
enum class CpGraphMode {
    OriginalGraph,
    Multigraph,
};
enum class CpSolveOutcome {
    Feasible,
    ProvenInfeasible,
    TimedOutUnknown,
    EarlyUnknown,
    ModelInvalid,
};

struct CpSymmetryOptions {
    bool first_pickup_vehicle_ordering = false;
    bool cor_43_identical_pickup_time_ordering = false;
};

struct CpRedundantOptions {
    bool terminal_balance = true;
    bool full_skip_reservoir = true;
    bool container_workload = false;
    bool aggregate_duration = true;
};

struct CpSatSolveOptions {
    int vehicle_count = 0;
    double time_limit_seconds = 0.0;
    int workers = 0;
    CpGraphMode graph_mode = CpGraphMode::OriginalGraph;
    CpSolveMode solve_mode = CpSolveMode::Satisfaction;
    double threshold_horizon_factor = 1.5;
    std::string log_file_path;
    CpRedundantOptions redundant;
    CpSymmetryOptions symmetry;
};

struct CpActionVisit {
    CpActionKind kind = CpActionKind::Pickup;
    int request_index = -1;
    int position = -1;
    std::int64_t start_time = -1;
};

struct CpActionRoute {
    int vehicle_index = -1;
    std::vector<CpActionVisit> actions;
    std::int64_t completion_time = -1;
};

struct CpModelBuildStats {
    std::int64_t candidate_action_arcs = 0;
    std::int64_t created_action_arcs = 0;
    std::int64_t core_precedence_pruned_arcs = 0;
    std::int64_t cor_40_pruned_arcs = 0;
    std::int64_t cor_41_pruned_arcs = 0;
    std::int64_t cor_43_constraint_count = 0;
    std::int64_t input_multigraph_edges = 0;
    std::int64_t created_multigraph_internal_arcs = 0;
    std::int64_t created_multigraph_start_arc_copies = 0;
    std::int64_t created_multigraph_end_arc_copies = 0;
    std::int64_t fixed_multigraph_connector_arcs = 0;
    std::int64_t state_continuity_constraint_count = 0;
    std::int64_t terminal_balance_constraint_count = 0;
    std::int64_t container_workload_constraint_count = 0;
    std::int64_t aggregate_duration_constraint_count = 0;
    std::int64_t first_pickup_symmetry_constraint_count = 0;
    bool full_skip_state_embedded = false;
};

struct CpFixedKSolveResult {
    CpSolveOutcome outcome = CpSolveOutcome::EarlyUnknown;
    CpGraphMode graph_mode = CpGraphMode::OriginalGraph;
    CpSolveMode solve_mode = CpSolveMode::Satisfaction;
    int raw_status = 0;
    std::string status_name;
    std::string termination_name;
    std::string error_message;
    std::string solution_info;
    double wall_time_seconds = 0.0;
    double configured_time_limit_seconds = 0.0;
    double threshold_horizon_factor = 0.0;
    std::int64_t model_horizon = 0;
    bool has_objective_value = false;
    double objective_value = 0.0;
    bool has_objective_bound = false;
    double best_objective_bound = 0.0;
    bool stopped_by_feasible_observer = false;
    bool stopped_by_bound_callback = false;
    bool hit_time_limit = false;
    bool early_unknown = false;
    std::int64_t conflicts = 0;
    std::int64_t branches = 0;
    CpModelBuildStats build_stats;
    std::vector<CpActionRoute> routes;
    std::vector<int> active_edge_ids;
};

const char* cp_solve_mode_name(CpSolveMode mode);
const char* cp_graph_mode_name(CpGraphMode mode);

CpSolveOutcome classify_cp_unknown_outcome(
    double configured_time_limit_seconds,
    double solver_wall_time_seconds
);

CpFixedKSolveResult solve_fixed_k_cp_sat(
    const SPDPData& data,
    const CpSatSolveOptions& options
);

CpFixedKSolveResult solve_fixed_k_cp_sat(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const CpSatSolveOptions& options
);

}  // namespace spdp

#endif
