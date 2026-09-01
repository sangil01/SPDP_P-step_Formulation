#ifndef SPDP_TWO_INDEX_MODEL_CORE_H
#define SPDP_TWO_INDEX_MODEL_CORE_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "ReadData.h"
#include "gurobi_c++.h"

namespace spdp::detail {

// Shared by every Gurobi model whose returned primal values or bounds are
// interpreted by the two-index/VI44 incumbent pipeline.  Keeping this value in
// one place prevents solver parameters and post-solve certification checks
// from silently using different numerical accuracies.
inline constexpr double kGurobiSolverTolerance = 1e-6;

enum class TwoIndexCoreObjective {
    None,
    OriginalCost,
    TravelCost,
    Duration,
    DurationPlusFixed,
};

struct TwoIndexCoreOptions {
    TwoIndexCoreObjective objective = TwoIndexCoreObjective::Duration;
    bool binary_y = false;
    bool add_time_constraints = false;
    double time_horizon = 0.0;  // 0 uses data.time_limit.
    bool enforce_route_duration_limit = true;
    double solver_time_limit = 0.0;
    int gurobi_threads = 1;
    bool output_enabled = false;
    std::string gurobi_log_path;
    std::string name_prefix = "two_index";
};

struct TwoIndexCoreModel {
    std::unique_ptr<GRBEnv> environment;
    std::unique_ptr<GRBModel> model;
    std::vector<GRBVar> y_vars;
    std::vector<std::map<State, GRBVar>> time_vars_by_node;
};

TwoIndexCoreModel build_two_index_model_core(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const TwoIndexCoreOptions& options
);

GRBVar two_index_time_var(
    const TwoIndexCoreModel& core,
    NodeId node_id,
    const State& state
);

}  // namespace spdp::detail

#endif
