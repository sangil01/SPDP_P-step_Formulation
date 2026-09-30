#ifndef SPDP_CUT_MANAGEMENT_H
#define SPDP_CUT_MANAGEMENT_H

#include <cstddef>
#include <string>
#include <vector>

namespace spdp {

// Where a dynamically separated cut family may run.
enum class CutSeparationScope {
    RootOnly,
    AdaptiveTree,
    FullTree,
};

// Management policy shared by every dynamically separated cut family.
//
// The families differ only in how one round produces candidates: the capacity
// blossom family solves a fresh odd-cut separation problem, the type-closed
// travel-time cover family screens the certified rows built during
// preprocessing. Everything after that (violation filter, duplicate removal,
// overlap filter, ordering by violation, per-round and total caps, tailing-off
// detection, separation time budget) is identical, so each family owns its own
// independent instance of this policy and its own stop state.
struct CutManagementOptions {
    CutSeparationScope scope = CutSeparationScope::AdaptiveTree;
    std::size_t root_max_rounds = 20;
    std::size_t root_max_cuts = 256;
    std::size_t root_max_per_round = 32;
    std::size_t tree_max_per_round = 8;
    // Adaptive tree: separate at every node while node_count is below
    // tree_dense_node_limit, afterwards only when the node index is a multiple
    // of tree_node_frequency.
    std::size_t tree_dense_node_limit = 50;
    std::size_t tree_node_frequency = 100;
    std::size_t max_total_cuts = 2000;
    // A candidate is usable only when its violation exceeds this value.
    double min_violation = 1e-4;
    // Stop tree separation once the accumulated separation time exceeds this
    // fraction of the solver runtime (0 disables the check).
    double max_separation_time_fraction = 0.15;
    // Inside one round, a candidate overlapping an already accepted candidate
    // of the same family by more than this Jaccard index is skipped.
    double max_overlap_jaccard = 0.9;
    // Root separation stops after this many consecutive rounds without a cut
    // (0 disables the check).
    std::size_t root_no_cut_round_limit = 2;
    // Root separation stops after this many consecutive rounds whose best
    // violation stays below root_low_violation (0 disables the check).
    double root_low_violation = 0.02;
    std::size_t root_low_violation_round_limit = 2;
};

// Per-family stop state of one solve.
struct CutManagementState {
    std::size_t root_rounds = 0;
    std::size_t no_cut_streak = 0;
    std::size_t low_violation_streak = 0;
    bool root_closed = false;
    bool tree_closed = false;
};

struct CutManagementStats {
    std::size_t rounds = 0;
    std::size_t root_rounds = 0;
    std::size_t cuts_added = 0;
    std::size_t root_cuts_added = 0;
    std::size_t duplicate_rejections = 0;
    std::size_t overlap_rejections = 0;
    bool stopped_by_total_cap = false;
    bool stopped_by_time_fraction = false;
    bool stopped_by_root_rounds = false;
    bool stopped_by_root_cuts = false;
    bool stopped_by_no_cut_streak = false;
    bool stopped_by_low_violation = false;
    double separation_seconds = 0.0;
};

// Cuts the family may return at this node; zero means skip the node.
std::size_t cut_management_budget(
    const CutManagementOptions& options,
    CutManagementState& state,
    CutManagementStats& stats,
    bool is_root,
    double node_count,
    double runtime_seconds
);

// Round counters and tailing-off streaks after one separation round.
void cut_management_finish_round(
    const CutManagementOptions& options,
    CutManagementState& state,
    CutManagementStats& stats,
    bool is_root,
    std::size_t cuts_added,
    double best_violation
);

// True when the adaptive tree schedule covers this node index. Root-only and
// full-tree scopes always return true; the adaptive scope separates while the
// node index stays below tree_dense_node_limit and afterwards only on multiples
// of tree_node_frequency.
bool cut_management_node_scheduled(
    const CutManagementOptions& options,
    double node_count
);

// True when the family cannot contribute a cut at this node, either because its
// stop state is reached or because the scope and the tree schedule skip the
// node. The node index is part of the test so a callback can decide before it
// copies the node relaxation.
bool cut_management_exhausted(
    const CutManagementOptions& options,
    const CutManagementState& state,
    const CutManagementStats& stats,
    bool is_root,
    double node_count
);

// Reopens the root phase of one family after the pre-MIP cut loop handed its
// rows to the main MIP. The round counters and the tailing-off streaks of the
// finished phase are reset, while the accumulated cut counts are kept so the
// root and total caps keep applying across both phases.
void cut_management_reopen_root(
    CutManagementState& state,
    CutManagementStats& stats
);

// Jaccard index of two sorted index sets.
double cut_set_overlap_jaccard(
    const std::vector<int>& left,
    const std::vector<int>& right
);

const char* cut_separation_scope_name(CutSeparationScope scope);
bool parse_cut_separation_scope(const std::string& value, CutSeparationScope& scope);
std::string cut_management_stop_reason(const CutManagementStats& stats);

}  // namespace spdp

#endif
