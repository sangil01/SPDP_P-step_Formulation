#ifndef SPDP_REQUEST_BLOCK_SEC_H
#define SPDP_REQUEST_BLOCK_SEC_H

#include <cstddef>
#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "ReadData.h"

namespace spdp {

struct RequestBlockSECRow {
    std::string name;
    double rhs = 0.0;
    std::vector<std::size_t> edge_ids;
};

// Builds request-block SEC rows for every I subseteq N satisfying
// 2 <= |I| <= min(max_subset_size, |N|).
std::vector<RequestBlockSECRow> build_request_block_sec_rows(
    const SPDPData& data,
    const MultiDiGraph& graph,
    std::size_t max_subset_size
);

}  // namespace spdp

#endif
