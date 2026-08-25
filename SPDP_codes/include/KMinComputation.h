#ifndef SPDP_K_MIN_COMPUTATION_H
#define SPDP_K_MIN_COMPUTATION_H

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
    bool add_time_constraints = false;
    double time_limit = 0.0;  // 0 means no time limit.
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

}  // namespace spdp

#endif
