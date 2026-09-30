#ifndef SPDP_TWO_INDEX_SOLVER_H
#define SPDP_TWO_INDEX_SOLVER_H

#include <string>
#include <vector>

#include "CapacityBlossomCuts.h"
#include "CutManagement.h"
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

// Where the root cuts of the dynamic families are separated.
enum class RootCutMode {
    // Solve a separate continuous LP first, run the cut loop there and copy the
    // resulting rows into the MIP as ordinary constraints.
    PreMip,
    // Separate inside the MIPNODE callback of the main MIP at the root node.
    Callback,
};

struct DirectTwoIndexOptions {
    DirectTwoIndexModelType model_type = DirectTwoIndexModelType::IP;
    DirectTwoIndexObjective objective = DirectTwoIndexObjective::OriginalCost;
    // Duration modelling of the main model. The big-M rows and the time-flow
    // rows are independent; at least one of them must be active.
    bool add_time_constraints = false;
    bool add_time_flow_formulation = false;
    bool time_flow_state_disaggregated = false;
    double solver_time_limit = 0.0;  // 0 means no time limit.
    int gurobi_threads = -1;         // Negative means the Gurobi default.
    std::string gurobi_log_path;
    // Separate Gurobi log of the pre-MIP root LP (PreMip mode only).
    std::string pre_mip_lp_gurobi_log_path;
    std::vector<double> initial_edge_start;
    PstepValidInequalityOptions valid_inequalities;
    CapacityBlossomOptions capacity_blossom;
    // Rows are precomputed by the caller (they need the restricted duration
    // IPs); the solver adds them statically or screens them as user cuts.
    TypeTravelTimeCoverOptions type_cover;
    std::vector<TypeTravelTimeCoverRow> type_cover_rows;
    RootCutMode root_cut_mode = RootCutMode::Callback;
    // Wall-clock budget of the pre-MIP root LP phase; always deducted from the
    // main MIP budget afterwards. 0 means no separate limit.
    double pre_mip_root_lp_time_limit = 30.0;
    // PreMip mode only: keep separating at the root node of the main MIP after
    // the pre-MIP cut loop copied its rows into the model. Ignored in Callback
    // mode, where the root is separated by the callback from the start.
    bool pre_mip_separate_main_mip_root = false;
    // Gurobi LP algorithms. -1 automatic, 0 primal, 1 dual, 2 barrier,
    // 4 deterministic concurrent (method only). Crossover -1 keeps the default
    // and 0 disables it; 0 is valid only together with barrier.
    int pre_mip_method = -1;
    int pre_mip_crossover = -1;
    int mip_method = -1;
    int mip_node_method = -1;
    int mip_crossover = -1;
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
    CapacityBlossomStats capacity_blossom_stats;
    TypeTravelTimeCoverStats type_cover_stats;
    int lp_cut_rounds = 0;
    int pre_mip_rounds = 0;
    int pre_mip_rows = 0;
    double pre_mip_seconds = 0.0;
    double pre_mip_lp_value_before = -1.0;
    double pre_mip_lp_value_after = -1.0;
    double mip_time_limit_used = 0.0;
    std::vector<double> edge_values;
};

const char* root_cut_mode_name(RootCutMode mode);
bool parse_root_cut_mode(const std::string& value, RootCutMode& mode);

DirectTwoIndexResult solve_direct_two_index_model(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const DirectTwoIndexOptions& options
);

}  // namespace spdp

#endif
