#ifndef SPDP_TWO_INDEX_SOLVER_H
#define SPDP_TWO_INDEX_SOLVER_H

#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "PstepValidInequality.h"
#include "ReadData.h"

namespace spdp {

enum class DirectTwoIndexModelType {
    IP,
    LP,
};

enum class DirectTwoIndexObjective {
    OriginalCost,       // variable travel cost + fixed vehicle cost
    TravelCost,         // variable travel cost only
    Duration,           // total route duration only
    DurationPlusFixed,  // total route duration + fixed vehicle cost
};

struct DirectTwoIndexOptions {
    DirectTwoIndexModelType model_type = DirectTwoIndexModelType::IP;
    DirectTwoIndexObjective objective = DirectTwoIndexObjective::OriginalCost;
    bool add_time_constraints = false;
    double solver_time_limit = 0.0;  // 0 means no time limit.
    int gurobi_threads = -1;         // Negative means the Gurobi default.
    std::string gurobi_log_path;
    std::vector<double> initial_edge_start;
    PstepValidInequalityOptions valid_inequalities;
};

struct DirectTwoIndexResult {
    DirectTwoIndexModelType model_type = DirectTwoIndexModelType::IP;
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
    double departure_flow = -1.0;
    int vehicle_count = 0;
    int variable_count = 0;
    int constraint_count = 0;
    int valid_inequality_count = 0;
    int fixed_vehicle_constraint_count = 0;
    std::vector<double> edge_values;
};

DirectTwoIndexResult solve_direct_two_index_model(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const DirectTwoIndexOptions& options
);

}  // namespace spdp

#endif
