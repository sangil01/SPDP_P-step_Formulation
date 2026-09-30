#ifndef SPDP_K_MIN_COMPUTATION_H
#define SPDP_K_MIN_COMPUTATION_H

#include <optional>
#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "ReadData.h"

namespace spdp {

enum class VI44SubproblemType {
    LP,
    IP,
};

struct VI44SubproblemOptions {
    VI44SubproblemType type = VI44SubproblemType::LP;
    int gurobi_threads = 1;
    double time_limit = 0.0;  // 0 means no time limit.
    bool rounded_bound_stop = true;
    bool dff_fs_enabled = false;
    std::vector<double> dff_fs_lambdas;
};

struct VI44VehicleAssignmentOptions {
    bool add_tsp_bound = false;
    bool add_container_bound = false;
    double time_limit = 0.0;  // 0 means no time limit.
};

struct VI44KMinResult;

struct VI44KMinOptions {
    bool use_cor = true;
    bool use_subproblem = false;
    bool use_vehicle_assignment = false;
    VI44SubproblemOptions subproblem;
    VI44VehicleAssignmentOptions vehicle_assignment;
    // Non-owning cache used by one solver run. When present, VI rows reuse the
    // k_min computation already performed for initial-incumbent generation.
    const VI44KMinResult* precomputed_result = nullptr;
};

struct VI44VehicleAssignmentVehicleResult {
    int vehicle_index = -1;
    int pickup_count = 0;
    int delivery_count = 0;
    double service_time = 0.0;
    double container_lb = -1.0;
    double tsp_lb = -1.0;
    double active_duration_lb = 0.0;
};

struct VI44DffSubproblemResult {
    std::string family;
    double parameter = 0.0;
    int status = 0;
    bool hit_time_limit = false;
    bool has_certified_bound = false;
    double objective_value = -1.0;
    double objective_bound = -1.0;
    double safe_lower_bound = -1.0;
    double numerical_tolerance = 0.0;
    double runtime_seconds = 0.0;
    int k_min = 0;
    bool stopped_by_rounded_bound = false;
    bool rounded_bound_certified = false;
    int certified_rounded_k = 0;
    double callback_objective_ub = -1.0;
    double callback_safe_objective_lb = -1.0;
    bool skipped_duplicate = false;
    double duplicate_of_parameter = -1.0;
};

struct VI44KMinResult {
    bool cor_enabled = false;
    int cor_k_min = 0;
    bool subproblem_enabled = false;
    int subproblem_k_min = 0;
    int subproblem_status = 0;
    bool subproblem_hit_time_limit = false;
    bool subproblem_has_certified_bound = false;
    double subproblem_objective_value = -1.0;
    double subproblem_objective_bound = -1.0;
    double subproblem_safe_lower_bound = -1.0;
    double subproblem_numerical_tolerance = 0.0;
    double subproblem_runtime_seconds = 0.0;
    bool subproblem_stopped_by_rounded_bound = false;
    bool subproblem_rounded_bound_certified = false;
    int subproblem_certified_rounded_k = 0;
    double subproblem_callback_objective_ub = -1.0;
    double subproblem_callback_safe_objective_lb = -1.0;
    bool dff_subproblem_fs_enabled = false;
    int dff_subproblem_k_min = 0;
    std::vector<VI44DffSubproblemResult> dff_subproblem_results;
    bool vehicle_assignment_enabled = false;
    int vehicle_assignment_k_min = 0;
    int vehicle_assignment_status = 0;
    bool vehicle_assignment_hit_time_limit = false;
    bool vehicle_assignment_has_certified_bound = false;
    bool vehicle_assignment_has_feasible_solution = false;
    bool vehicle_assignment_optimal = false;
    int vehicle_assignment_tested_vehicle_count = 0;
    double vehicle_assignment_runtime_seconds = 0.0;
    std::vector<VI44VehicleAssignmentVehicleResult> vehicle_assignment_vehicles;
    int selected_k_min = 0;
};

VI44KMinResult compute_vi44_k_min(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const VI44KMinOptions& options
);

// Minimum-total-duration IP with the rounded-bound early stop: the solve is
// aborted as soon as ceil(safe_LB / T) == ceil(exact_UB / T), because the
// route-count bound ceil(eta / T) is then certified without an optimal eta.
struct DurationRoundedBoundOptions {
    double solver_time_limit = 0.0;  // 0 means no time limit.
    bool rounded_bound_stop = true;
    int gurobi_threads = 1;
    std::string name_prefix = "duration_bound";
};

struct DurationRoundedBoundResult {
    int status = 0;
    bool hit_time_limit = false;
    bool stopped_by_rounded_bound = false;
    // Exact duration of the best integer solution (sum of selected edge times),
    // -1 when none was found.
    double exact_upper_bound = -1.0;
    // Solver objective bound and the tolerance-corrected value <= eta.
    double objective_bound = -1.0;
    double safe_lower_bound = -1.0;
    bool has_lower_bound = false;
    // ceil(safe_lower_bound / T): always a valid route-count bound.
    int rounded_lower_bound = 0;
    // True when ceil(safe_lower_bound / T) == ceil(exact_upper_bound / T).
    bool certified = false;
    double runtime_seconds = 0.0;
};

DurationRoundedBoundResult solve_duration_rounded_bound(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const DurationRoundedBoundOptions& options
);

// Returns the common rounded duration bound only when the safe lower bound
// and an exact feasible incumbent objective certify the same value.
std::optional<int> certified_rounded_duration_bound(
    double exact_incumbent_upper_bound,
    double solver_objective_lower_bound,
    double route_time_limit
);

}  // namespace spdp

#endif
