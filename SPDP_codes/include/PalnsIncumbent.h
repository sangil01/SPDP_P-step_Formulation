#ifndef SPDP_PALNS_INCUMBENT_H
#define SPDP_PALNS_INCUMBENT_H

#include <memory>
#include <string>
#include <vector>

#include "CpSatSpdpSolver.h"
#include "GenMultiGraph.h"
#include "ReadData.h"

namespace spdp {

// Pattern-based ALNS for the fixed-K incumbent problem (note 15).
enum class PalnsSkeletonSet { All, S3S5 };
enum class PalnsCompositeFilter { DeltaPositive, All };
enum class PalnsConstructOrder { Lpt, Random, Regret, Mixed };
enum class PalnsLsFrequency { Every, OnImprove, EveryM };

struct PalnsOptions {
    unsigned int seed = 1;
    // 1. Pattern generation
    int composite_kappa = 1;
    PalnsSkeletonSet composite_skeletons = PalnsSkeletonSet::All;
    double sim_s1 = 1.0;
    double sim_s2 = 1.0;
    double sim_s3 = 1.0;
    // 2. Construction
    int construct_starts = 10;
    PalnsCompositeFilter construct_composite_filter =
        PalnsCompositeFilter::DeltaPositive;
    PalnsConstructOrder construct_order = PalnsConstructOrder::Mixed;
    // 3. PALNS search
    double eta = 0.01;
    double kmax_ratio = 0.1;
    int kmax_min = 2;
    // 4. Local search
    PalnsLsFrequency ls_frequency = PalnsLsFrequency::Every;
    int ls_every_m = 10;
    bool ls_first_improvement = true;
    // 6. Phase II and termination
    bool phase2_enable = true;
    int phase2_max_noimp = 2000;
    // Search-mode stop (a): stop as soon as a feasible solution with
    // original cost <= target_cost is found. Negative disables it (mode b).
    double target_cost = -1.0;
    // Simulated-annealing parameters (tuning knobs, not part of note 4.2).
    // Theta0 = sa_start_ratio * f(phase start) / ln 2; geometric cooling to
    // sa_end_ratio * Theta0 within each cycle of sa_cycle_iters iterations;
    // at a cycle end the search restarts from the best solution when
    // sa_restart_best is set.
    double sa_start_ratio = 0.05;
    double sa_end_ratio = 1e-3;
    long long sa_cycle_iters = 5000;
    bool sa_restart_best = true;
    // Trajectory relocate is applied only to candidates whose objective after
    // the light local search is within ls_deep_ratio of the best objective.
    double ls_deep_ratio = 0.01;
    // Independent portfolio threads with different seeds. The first thread
    // that terminates (target, Phase II stop, or time) stops the others; the
    // best feasible solution over all threads is returned.
    int threads = 1;
};

struct PalnsFixedKResult {
    bool feasible = false;
    bool hit_time_limit = false;
    bool reached_target = false;
    bool construction_feasible = false;
    double runtime_seconds = 0.0;
    double construction_seconds = 0.0;
    double time_to_feasible = -1.0;
    double time_to_target = -1.0;
    double best_violation = 0.0;
    double best_duration = 0.0;
    double best_travel_cost = 0.0;
    double best_original_cost = 0.0;
    long long phase1_iterations = 0;
    long long phase2_iterations = 0;
    int construction_starts_done = 0;
    std::string stop_reason;
    std::vector<CpActionRoute> routes;
    std::string log;
};

class PalnsIncumbentSolver {
public:
    PalnsIncumbentSolver(const SPDPData& data, const PalnsOptions& options);
    ~PalnsIncumbentSolver();
    PalnsIncumbentSolver(const PalnsIncumbentSolver&) = delete;
    PalnsIncumbentSolver& operator=(const PalnsIncumbentSolver&) = delete;

    double pattern_generation_seconds() const;
    std::string pattern_summary() const;
    PalnsFixedKResult solve_fixed_k(int vehicle_count, double time_limit_seconds);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

struct PalnsMappedIncumbent {
    bool success = false;
    std::string error_message;
    std::vector<int> active_edge_ids;
    double total_duration = 0.0;
    double total_original_cost = 0.0;
};

// Maps canonicalized PALNS action routes to multigraph edges. Among parallel
// edges with the same endpoint states, an exact match is preferred; otherwise
// an edge that is no slower and no more expensive, otherwise the fastest one.
PalnsMappedIncumbent map_palns_routes_to_multigraph(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const std::vector<CpActionRoute>& routes,
    int expected_vehicle_count
);

const char* palns_skeleton_set_name(PalnsSkeletonSet value);
const char* palns_composite_filter_name(PalnsCompositeFilter value);
const char* palns_construct_order_name(PalnsConstructOrder value);
const char* palns_ls_frequency_name(PalnsLsFrequency value);
bool parse_palns_skeleton_set(const std::string& text, PalnsSkeletonSet& value);
bool parse_palns_composite_filter(const std::string& text, PalnsCompositeFilter& value);
bool parse_palns_construct_order(const std::string& text, PalnsConstructOrder& value);
bool parse_palns_ls_frequency(const std::string& text, PalnsLsFrequency& value);

}  // namespace spdp

#endif
