#ifndef SPDP_CP_INCUMBENT_ADAPTER_H
#define SPDP_CP_INCUMBENT_ADAPTER_H

#include <string>
#include <vector>

#include "CpSatSpdpSolver.h"
#include "GenMultiGraph.h"
#include "ReadData.h"

namespace spdp {

struct CpMappedIncumbent {
    bool success = false;
    std::string error_message;
    std::vector<int> active_edge_ids;
    double total_duration = 0.0;
    double total_original_cost = 0.0;
};

CpMappedIncumbent map_cp_incumbent_to_multigraph(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const std::vector<CpActionRoute>& routes
);

CpMappedIncumbent validate_multigraph_cp_incumbent(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const std::vector<int>& active_edge_ids,
    int expected_vehicle_count
);

}  // namespace spdp

#endif
