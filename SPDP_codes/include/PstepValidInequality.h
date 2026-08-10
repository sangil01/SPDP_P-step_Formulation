#ifndef SPDP_PSTEP_VALID_INEQUALITY_H
#define SPDP_PSTEP_VALID_INEQUALITY_H

#include <cstddef>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

#include "GenMultiGraph.h"
#include "ReadData.h"

namespace spdp {

struct VI44KMinResult;

enum class PstepValidInequalitySense {
    GreaterEqual,
    LessEqual,
};

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

struct PstepValidInequalityRow {
    std::string name;
    PstepValidInequalitySense sense = PstepValidInequalitySense::GreaterEqual;
    double rhs = 0.0;
    std::vector<std::pair<std::size_t, double>> edge_terms;
};

struct PstepValidInequalityOptions {
    bool add_vi_35 = false;
    bool add_vi_36_combined = false;
    std::size_t vi_36_subset_max_size = 0;
    bool add_vi_request_block_sec = false;
    std::size_t vi_request_block_sec_max_size = 0;
    bool add_vi_44 = false;
    VI44KMinOptions vi_44_k_min_options;
    std::ostream* log_stream = nullptr;
};

enum class DirectTwoIndexObjective {
    OriginalCost,       // variable travel cost + fixed vehicle cost
    TravelCost,         // variable travel cost only
    Duration,           // total route duration only
    DurationPlusFixed,  // total route duration + fixed vehicle cost
};

struct DirectTwoIndexOptions {
    DirectTwoIndexObjective objective = DirectTwoIndexObjective::OriginalCost;
    bool add_time_constraints = false;
    double solver_time_limit = 0.0;  // 0 means no time limit.
    int gurobi_threads = -1;         // Negative means the Gurobi default.
    std::string gurobi_log_path;
    std::vector<double> initial_edge_start;
    PstepValidInequalityOptions valid_inequalities;
};

struct DirectTwoIndexResult {
    int status = 0;
    bool hit_time_limit = false;
    bool solved_to_optimality = false;
    bool has_feasible_solution = false;
    bool has_certified_bound = false;
    double objective_value = -1.0;
    double objective_bound = -1.0;
    double gap_percent = -1.0;
    double runtime_seconds = 0.0;
    double total_duration = -1.0;
    double total_duration_plus_fixed = -1.0;
    double total_travel_cost = -1.0;
    double total_original_cost = -1.0;
    int vehicle_count = 0;
    int variable_count = 0;
    int constraint_count = 0;
    int valid_inequality_count = 0;
    std::vector<double> edge_values;
};

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

VI44KMinResult compute_vi44_k_min(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const VI44KMinOptions& options
);

// Builds every enabled p-step valid inequality as a sparse edge row.
// A maximum subset size of zero disables the corresponding subset family.
std::vector<PstepValidInequalityRow> build_pstep_valid_inequality_rows(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const PstepValidInequalityOptions& options
);

// Solves the binary two-index formulation directly. This is the P=0
// enumeration path; its core constraints are shared with the VI-44 auxiliary
// LP/IP, while the enabled valid inequalities are added only to this main model.
DirectTwoIndexResult solve_direct_two_index_ip(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const DirectTwoIndexOptions& options
);

// Finds the first feasible duration-model solution using exactly k vehicles.
// Fixed-K and SolutionLimit are intentionally internal to this preprocessing
// routine rather than exposed as general enumeration CLI controls.
InitialIncumbentSolveResult solve_fixed_k_duration_initial_incumbent(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const InitialIncumbentSolveOptions& options
);

}  // namespace spdp

#endif
