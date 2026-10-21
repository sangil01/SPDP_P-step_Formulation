#ifndef SPDP_TIME_ROW_FREE_ANALYSIS_H
#define SPDP_TIME_ROW_FREE_ANALYSIS_H

#include <cstddef>
#include <vector>

#include "GenMultiGraph.h"
#include "ReadData.h"

namespace spdp {

struct TimeRowFreeComponent {
    std::vector<int> edge_ids;
    double duration = 0.0;
};

struct TimeRowFreeAnalysis {
    std::vector<TimeRowFreeComponent> depot_routes;
    std::vector<TimeRowFreeComponent> disconnected_subtours;
    int duration_violating_route_count = 0;
    double maximum_depot_route_duration = 0.0;
    double maximum_duration_excess = 0.0;
    std::size_t selected_edge_count = 0U;
};

TimeRowFreeAnalysis analyze_time_row_free_solution(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const std::vector<int>& selected_edge_ids
);

}  // namespace spdp

#endif
