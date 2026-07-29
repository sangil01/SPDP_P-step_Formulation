#ifndef SPDP_PSTEP_VALID_INEQUALITY_H
#define SPDP_PSTEP_VALID_INEQUALITY_H

#include <cstddef>
#include <iosfwd>
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

enum class VI44KMinMode {
    COR,
    SubLP,
    SubIP,
};

struct VI44KMinResult {
    int cor_k_min = 0;
    int sub_lp_k_min = 0;
    int selected_k_min = 0;
    int sub_lp_status = 0;
    bool sub_lp_hit_time_limit = false;
    bool sub_lp_has_certified_bound = false;
    double sub_lp_objective_value = -1.0;
    double sub_lp_objective_bound = -1.0;
    double sub_lp_safe_lower_bound = -1.0;
    double sub_lp_numerical_tolerance = 0.0;
    double sub_lp_runtime_seconds = 0.0;
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
    VI44KMinMode vi_44_k_min_mode = VI44KMinMode::COR;
    // Applied to either auxiliary mode; 0 means no time limit.
    double vi_44_sub_lp_time_limit = 0.0;
    std::ostream* log_stream = nullptr;
};

VI44KMinResult compute_vi44_k_min(
    const SPDPData& data,
    const MultiDiGraph& graph,
    VI44KMinMode mode,
    double sub_lp_time_limit
);

// Builds every enabled p-step valid inequality as a sparse edge row.
// A maximum subset size of zero disables the corresponding subset family.
std::vector<PstepValidInequalityRow> build_pstep_valid_inequality_rows(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const PstepValidInequalityOptions& options
);

}  // namespace spdp

#endif
