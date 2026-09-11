#ifndef SPDP_TYPE_TRAVEL_TIME_COVER_CUTS_H
#define SPDP_TYPE_TRAVEL_TIME_COVER_CUTS_H

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "PstepValidInequality.h"
#include "ReadData.h"
#include "ShortestTravelTimes.h"

namespace spdp {

// Type-closed travel-time cover cuts.
//
// For a set H of skip types let S_H be the pickup and delivery nodes of every
// request whose type lies in H. Deleting all other actions from a feasible
// route yields a feasible route of the restricted instance that contains only
// the H requests, provided travel times satisfy the triangle inequality; the
// restricted instance therefore uses the metric closure (shortest travel
// times). With eta(H) the minimum total route duration of that instance,
// every integer solution satisfies
//
//     sum_{e in delta^-(S_H)} y_e >= rho(H) = ceil(eta(H) / T).
//
// eta(H) is not solved to optimality: the duration IP stops as soon as the
// rounded safe lower bound and the rounded exact incumbent agree, which
// certifies rho(H). Uncertified masks are never turned into rows.

enum class TypeTravelTimeCoverMode {
    Static,    // add every certified row to the model before solving
    Screened,  // add certified rows as user cuts only when violated
};

struct TypeTravelTimeCoverOptions {
    bool enabled = false;
    TypeTravelTimeCoverMode mode = TypeTravelTimeCoverMode::Screened;
    std::size_t max_type_set_size = 0;  // 0 = every non-empty union of types
    double subproblem_time_limit = 20.0;
    double total_time_limit = 120.0;
    std::size_t max_per_round = 8;
    bool add_time_flow_formulation = true;
    bool time_flow_state_disaggregated = true;
    bool add_time_constraints = false;  // big-M rows in the duration IP
    int gurobi_threads = 1;
    // Rows with rho(H) below this value are not added (rho = 1 is an ordinary
    // connectivity requirement).
    int min_rho = 1;
    double violation_tolerance = 1e-4;
    std::string cache_dir;  // empty disables the certificate cache
    std::ostream* log_stream = nullptr;
};

struct TypeTravelTimeCoverEntry {
    std::uint64_t type_mask = 0;
    std::vector<int> types;
    std::size_t request_count = 0;
    std::size_t restricted_edge_count = 0;
    std::uint64_t restricted_graph_fingerprint = 0;
    double duration_lower_bound = -1.0;  // safe LB on eta(H)
    double duration_upper_bound = -1.0;  // exact incumbent duration
    int rho = 0;                          // ceil(LB / T), valid whenever LB is
    bool has_lower_bound = false;
    bool certified = false;               // rounded LB == rounded UB
    bool stopped_by_rounded_bound = false;
    bool hit_time_limit = false;
    bool skipped = false;
    std::string skip_reason;
    bool from_cache = false;
    int status = 0;
    double runtime_seconds = 0.0;
};

struct TypeTravelTimeCoverResult {
    std::vector<int> types;  // bit i of a mask <-> types[i]
    std::uint64_t instance_fingerprint = 0;
    std::uint64_t shortest_time_fingerprint = 0;
    std::size_t non_metric_pair_count = 0;
    double max_travel_time_reduction = 0.0;
    std::vector<TypeTravelTimeCoverEntry> entries;
    double total_seconds = 0.0;
    bool total_time_limit_hit = false;
    std::size_t cache_hits = 0;
};

struct TypeTravelTimeCoverRow {
    std::uint64_t type_mask = 0;
    int rho = 0;
    std::vector<NodeId> nodes;  // S_H
    PstepValidInequalityRow row;
};

struct TypeTravelTimeCoverStats {
    std::size_t rows_available = 0;
    std::size_t rows_added = 0;
    std::size_t root_rows_added = 0;
    std::size_t screening_rounds = 0;
    double screening_seconds = 0.0;
};

// Sorted distinct container types of the instance.
std::vector<int> instance_container_types(const SPDPData& data);

// Restricted instance: requests with a type in the mask, travel times replaced
// by the metric closure, every other field unchanged.
SPDPData build_type_restricted_instance(
    const SPDPData& data,
    const std::vector<std::vector<double>>& closure_travel_time,
    const std::vector<int>& types,
    std::uint64_t type_mask
);

TypeTravelTimeCoverResult compute_type_travel_time_cover_bounds(
    const SPDPData& data,
    const TypeTravelTimeCoverOptions& options
);

// One row per certified entry with rho >= options.min_rho, skipping the full
// type set (its row is the departure-count bound already covered by VI44).
std::vector<TypeTravelTimeCoverRow> build_type_travel_time_cover_rows(
    const SPDPData& data,
    const MultiDiGraph& main_graph,
    const TypeTravelTimeCoverResult& result,
    int min_rho
);

// Indices of rows not yet added that are violated by y, sorted by decreasing
// violation and truncated to max_rows.
std::vector<std::size_t> screen_type_travel_time_cover_rows(
    const std::vector<TypeTravelTimeCoverRow>& rows,
    const std::vector<double>& y_values,
    const std::vector<bool>& already_added,
    double violation_tolerance,
    std::size_t max_rows
);

void write_type_travel_time_cover_log(
    std::ostream& out,
    const TypeTravelTimeCoverResult& result
);

const char* type_travel_time_cover_mode_name(TypeTravelTimeCoverMode mode);

}  // namespace spdp

#endif

