#include "TestSupport.h"

#include "CpIncumbentAdapter.h"
#include "CpSatSpdpSolver.h"
#include "GenMultiGraph.h"

namespace {

spdp::SPDPData batching_data() {
    spdp::SPDPData data;
    data.fixed_vehicle_cost = 500;
    data.time_pickup = 2;
    data.time_empty = 3;
    data.time_delivery = 2;
    data.time_limit = 100;
    data.locations = 3;
    data.requests = {
        spdp::Request{0, 2, 4},
        spdp::Request{1, 2, 4},
    };
    data.time = {{0, 2, 5, 0}, {2, 0, 4, 0}, {5, 4, 0, 0}, {0, 0, 0, 0}};
    data.distance = data.time;
    return data;
}

}  // namespace

SPDP_TEST(cp_adapter) {
    const auto data = batching_data();
    spdp::CpSatSolveOptions options;
    options.vehicle_count = 1;
    options.time_limit_seconds = 5;
    const auto cp = spdp::solve_fixed_k_cp_sat(data, options);
    SPDP_CHECK_EQ(cp.outcome, spdp::CpSolveOutcome::Feasible);

    for (bool pruning : {false, true}) {
        spdp::GraphBuildOptions graph_options;
        graph_options.prune_symmetry_40 = pruning;
        graph_options.prune_symmetry_41 = pruning;
        const auto graph = spdp::build_multigraph(data, graph_options);
        const auto mapped = spdp::map_cp_incumbent_to_multigraph(
            data, graph, cp.routes
        );
        SPDP_CHECK(mapped.success);
        SPDP_CHECK(!mapped.active_edge_ids.empty());
    }

    options.vehicle_count = 2;
    const auto two_routes = spdp::solve_fixed_k_cp_sat(data, options);
    SPDP_CHECK_EQ(two_routes.outcome, spdp::CpSolveOutcome::Feasible);
    const auto graph = spdp::build_multigraph(data);
    const auto mapped_two_routes = spdp::map_cp_incumbent_to_multigraph(
        data, graph, two_routes.routes
    );
    SPDP_CHECK(mapped_two_routes.success);
}
