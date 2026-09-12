#ifndef SPDP_TWO_INDEX_SOLVER_H
#define SPDP_TWO_INDEX_SOLVER_H

#include <string>
#include <vector>

#include "CapacityBlossomCuts.h"
#include "ConnectivityCuts.h"
#include "GenMultiGraph.h"
#include "PstepValidInequality.h"
#include "ReadData.h"
#include "TypeTravelTimeCoverCuts.h"

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
    bool add_time_flow_formulation = false;
    bool time_flow_state_disaggregated = false;
    double solver_time_limit = 0.0;  // 0 means no time limit.
    int gurobi_threads = -1;         // Negative means the Gurobi default.
    std::string gurobi_log_path;
    std::vector<double> initial_edge_start;
    PstepValidInequalityOptions valid_inequalities;
    ConnectivityCutOptions connectivity_cuts;
    CapacityBlossomOptions capacity_blossom;
    // Rows are precomputed by the caller (they need the restricted duration
    // IPs); the solver adds them statically or screens them as user cuts.
    TypeTravelTimeCoverOptions type_cover;
    std::vector<TypeTravelTimeCoverRow> type_cover_rows;
    // Root LP cutting-plane phase: before the MIP starts, separate every enabled
    // family to closure on the LP relaxation and add the resulting rows to the
    // MIP as constraints. On large instances Gurobi's own root processing
    // leaves almost no MIPNODE callbacks inside the time limit, so this is
    // where the families act.
    bool root_lp_cut_phase = false;
    int root_lp_cut_phase_max_rounds = 50;
    double root_lp_cut_phase_time_limit = 30.0;  // seconds, 0 = no limit
    bool root_lp_cut_phase_deduct_time = true;    // subtract the phase time from the MIP limit
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
    int time_flow_infeasible_edge_count = 0;
    ConnectivityCutStats connectivity_cut_stats;
    CapacityBlossomStats capacity_blossom_stats;
    TypeTravelTimeCoverStats type_cover_stats;
    int lp_cut_rounds = 0;
    int root_lp_cut_phase_rounds = 0;
    int root_lp_cut_phase_rows = 0;
    double root_lp_cut_phase_seconds = 0.0;
    double root_lp_value_before = -1.0;
    double root_lp_value_after = -1.0;
    double mip_time_limit_used = 0.0;
    std::vector<double> edge_values;
};

DirectTwoIndexResult solve_direct_two_index_model(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const DirectTwoIndexOptions& options
);

}  // namespace spdp

#endif
