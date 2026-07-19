#ifndef SPDP_PSTEP_VALID_INEQUALITY_H
#define SPDP_PSTEP_VALID_INEQUALITY_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "GenMultiGraph.h"
#include "ReadData.h"

namespace spdp {

enum class PstepValidInequalitySense {
    GreaterEqual,
    LessEqual,
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
};

// Builds every enabled p-step valid inequality as a sparse edge row.
// A maximum subset size of zero disables the corresponding subset family.
std::vector<PstepValidInequalityRow> build_pstep_valid_inequality_rows(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const PstepValidInequalityOptions& options
);

}  // namespace spdp

#endif
