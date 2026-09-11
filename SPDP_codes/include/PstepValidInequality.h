#ifndef SPDP_PSTEP_VALID_INEQUALITY_H
#define SPDP_PSTEP_VALID_INEQUALITY_H

#include <cstddef>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

#include "GenMultiGraph.h"
#include "KMinComputation.h"
#include "ReadData.h"

namespace spdp {

enum class PstepValidInequalitySense {
    GreaterEqual,
    LessEqual,
    Equal,
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
    // Legacy COR projection: an edge whose embedded treatment lies between two
    // pickups (or two deliveries) of the subset is counted as leaving and
    // re-entering the subset. On the action-based multigraph emptying keeps the
    // slot occupied, so pickups/deliveries joined through a treatment still form
    // one capacity block; false counts only the real graph boundary, which
    // strengthens the pickup/delivery rows. Treatment-location rows are unaffected.
    bool vi_36_treatment_boundary = true;
    bool add_vi_request_block_sec = false;
    std::size_t vi_request_block_sec_max_size = 0;
    bool add_vi_44 = false;
    VI44KMinOptions vi_44_k_min_options;
    // This is a conditional model restriction, not a globally valid inequality.
    // It lives here to reuse the sparse departure-edge row representation.
    bool add_fixed_vehicle_number = false;
    int fixed_vehicle_number = 0;
    std::ostream* log_stream = nullptr;
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
