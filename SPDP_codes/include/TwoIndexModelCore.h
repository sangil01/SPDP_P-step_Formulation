#ifndef SPDP_TWO_INDEX_MODEL_CORE_H
#define SPDP_TWO_INDEX_MODEL_CORE_H

#include <memory>
#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "ReadData.h"
#include "gurobi_c++.h"

namespace spdp::detail {

enum class TwoIndexCoreObjective {
    OriginalCost,
    TravelCost,
    Duration,
    DurationPlusFixed,
};

struct TwoIndexCoreOptions {
    TwoIndexCoreObjective objective = TwoIndexCoreObjective::Duration;
    bool binary_y = false;
    bool add_time_constraints = false;
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
};

TwoIndexCoreModel build_two_index_model_core(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const TwoIndexCoreOptions& options
);

}  // namespace spdp::detail

#endif
