#ifndef SPDP_VI36_LOCATION_SUBSET_H
#define SPDP_VI36_LOCATION_SUBSET_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "GenMultiGraph.h"

namespace spdp {

struct VI36CombinedRow {
    std::string name;
    double rhs = 0.0;
    std::vector<std::pair<std::size_t, int>> edge_terms;
};

// Builds the combined pickup/treatment/delivery location-subset rows for
// 1 <= |S| <= min(max_subset_size, |O_X|).
std::vector<VI36CombinedRow> build_vi36_combined_location_subset_rows(
    const MultiDiGraph& graph,
    std::size_t max_subset_size,
    bool skip_full_location_sets
);

}  // namespace spdp

#endif
