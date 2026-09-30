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
#include "CutManagement.h"

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
// rounded safe lower bound and the rounded exact incumbent agree, which pins
// rho(H) exactly. When the two disagree the row still uses the rounded safe
// lower bound, because safe_LB <= eta(H) gives
//
//     ceil(safe_LB / T) <= ceil(eta(H) / T) = rho(H),
//
// so the row is a weaker but still valid member of the same family. Only the
// solver objective bound ever feeds a right-hand side; an incumbent duration
// is used to detect the exact case and never as a right-hand side itself. The
// duration IP always uses the big-M time formulation and is always resolved
// from scratch; no certificate is ever read from disk.

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
    int gurobi_threads = 1;
    // Rows with rho(H) below this value are not added (rho = 1 is an ordinary
    // connectivity requirement).
    int min_rho = 1;
    // Scope, caps, violation thresholds and tailing-off rules. Screened mode
    // only; a static family is added to the model before the solve starts.
    CutManagementOptions management;
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
};

struct TypeTravelTimeCoverRow {
    std::uint64_t type_mask = 0;
    int rho = 0;
    // True when rho was pinned to the exact minimum route count of the
    // restricted instance; false when it is the rounded safe lower bound of a
    // duration IP that ran out of time. Both cases are valid rows.
    bool certified = false;
    std::vector<NodeId> nodes;  // S_H
    PstepValidInequalityRow row;
};

struct TypeTravelTimeCoverStats {
    CutManagementStats management;
    std::size_t rows_available = 0;
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

// One row per entry that produced a solver objective bound and reaches
// rho >= options.min_rho, skipping the full type set (its row is the
// departure-count bound already covered by VI44). Entries whose duration IP
// hit its time limit still contribute the rounded safe lower bound.
std::vector<TypeTravelTimeCoverRow> build_type_travel_time_cover_rows(
    const SPDPData& data,
    const MultiDiGraph& main_graph,
    const TypeTravelTimeCoverResult& result,
    int min_rho
);

// Indices of rows not yet added whose violation exceeds min_violation, sorted
// by decreasing violation, filtered against overlapping selections of the same
// round and truncated to max_rows.
std::vector<std::size_t> screen_type_travel_time_cover_rows(
    const std::vector<TypeTravelTimeCoverRow>& rows,
    const std::vector<double>& y_values,
    const std::vector<bool>& already_added,
    double min_violation,
    double max_overlap_jaccard,
    std::size_t max_rows,
    double& best_violation,
    std::size_t& overlap_rejections
);

void write_type_travel_time_cover_log(
    std::ostream& out,
    const TypeTravelTimeCoverResult& result
);

const char* type_travel_time_cover_mode_name(TypeTravelTimeCoverMode mode);

}  // namespace spdp

#endif
