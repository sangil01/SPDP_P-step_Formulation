#include "CutManagement.h"

#include <algorithm>
#include <sstream>

namespace spdp {

bool cut_management_node_scheduled(
    const CutManagementOptions& options,
    double node_count
) {
    if (options.scope != CutSeparationScope::AdaptiveTree) {
        return true;
    }
    const auto node = static_cast<unsigned long long>(node_count + 0.5);
    const bool dense_phase = node < options.tree_dense_node_limit;
    const bool periodic = options.tree_node_frequency > 0U &&
        node % options.tree_node_frequency == 0ULL;
    return dense_phase || periodic;
}

std::size_t cut_management_budget(
    const CutManagementOptions& options,
    CutManagementState& state,
    CutManagementStats& stats,
    bool is_root,
    double node_count,
    double runtime_seconds
) {
    if (stats.cuts_added >= options.max_total_cuts) {
        stats.stopped_by_total_cap = true;
        state.root_closed = true;
        state.tree_closed = true;
        return 0;
    }
    const std::size_t total_remaining = options.max_total_cuts - stats.cuts_added;

    if (is_root) {
        if (state.root_closed) {
            return 0;
        }
        if (state.root_rounds >= options.root_max_rounds) {
            state.root_closed = true;
            stats.stopped_by_root_rounds = true;
            return 0;
        }
        if (stats.root_cuts_added >= options.root_max_cuts) {
            state.root_closed = true;
            stats.stopped_by_root_cuts = true;
            return 0;
        }
        return std::min(
            std::min(options.root_max_per_round, options.root_max_cuts - stats.root_cuts_added),
            total_remaining
        );
    }

    if (options.scope == CutSeparationScope::RootOnly || state.tree_closed) {
        return 0;
    }
    if (options.max_separation_time_fraction > 0.0 && runtime_seconds > 5.0 &&
        stats.separation_seconds > options.max_separation_time_fraction * runtime_seconds) {
        stats.stopped_by_time_fraction = true;
        state.tree_closed = true;
        return 0;
    }
    if (!cut_management_node_scheduled(options, node_count)) {
        return 0;
    }
    return std::min(options.tree_max_per_round, total_remaining);
}

void cut_management_finish_round(
    const CutManagementOptions& options,
    CutManagementState& state,
    CutManagementStats& stats,
    bool is_root,
    std::size_t cuts_added,
    double best_violation
) {
    ++stats.rounds;
    stats.cuts_added += cuts_added;
    if (!is_root) {
        return;
    }
    ++stats.root_rounds;
    ++state.root_rounds;
    stats.root_cuts_added += cuts_added;

    if (cuts_added == 0U) {
        ++state.no_cut_streak;
    } else {
        state.no_cut_streak = 0;
    }
    if (best_violation < options.root_low_violation) {
        ++state.low_violation_streak;
    } else {
        state.low_violation_streak = 0;
    }

    if (options.root_no_cut_round_limit > 0U &&
        state.no_cut_streak >= options.root_no_cut_round_limit) {
        state.root_closed = true;
        stats.stopped_by_no_cut_streak = true;
    }
    if (options.root_low_violation_round_limit > 0U &&
        state.low_violation_streak >= options.root_low_violation_round_limit) {
        state.root_closed = true;
        stats.stopped_by_low_violation = true;
    }
}

bool cut_management_exhausted(
    const CutManagementOptions& options,
    const CutManagementState& state,
    const CutManagementStats& stats,
    bool is_root,
    double node_count
) {
    if (stats.cuts_added >= options.max_total_cuts) {
        return true;
    }
    if (is_root) {
        return state.root_closed;
    }
    if (options.scope == CutSeparationScope::RootOnly || state.tree_closed) {
        return true;
    }
    return !cut_management_node_scheduled(options, node_count);
}

void cut_management_reopen_root(
    CutManagementState& state,
    CutManagementStats& stats
) {
    state.root_closed = false;
    state.root_rounds = 0;
    state.no_cut_streak = 0;
    state.low_violation_streak = 0;
    stats.stopped_by_root_rounds = false;
    stats.stopped_by_root_cuts = false;
    stats.stopped_by_no_cut_streak = false;
    stats.stopped_by_low_violation = false;
}

double cut_set_overlap_jaccard(
    const std::vector<int>& left,
    const std::vector<int>& right
) {
    if (left.empty() || right.empty()) {
        return 0.0;
    }
    std::vector<int> intersection;
    std::set_intersection(
        left.begin(), left.end(),
        right.begin(), right.end(),
        std::back_inserter(intersection)
    );
    const std::size_t union_size =
        left.size() + right.size() - intersection.size();
    if (union_size == 0U) {
        return 0.0;
    }
    return static_cast<double>(intersection.size()) / static_cast<double>(union_size);
}

const char* cut_separation_scope_name(CutSeparationScope scope) {
    switch (scope) {
        case CutSeparationScope::RootOnly: return "root-only";
        case CutSeparationScope::AdaptiveTree: return "adaptive-tree";
        case CutSeparationScope::FullTree: return "full-tree";
    }
    return "unknown";
}

bool parse_cut_separation_scope(const std::string& value, CutSeparationScope& scope) {
    if (value == "root-only") {
        scope = CutSeparationScope::RootOnly;
        return true;
    }
    if (value == "adaptive-tree") {
        scope = CutSeparationScope::AdaptiveTree;
        return true;
    }
    if (value == "full-tree") {
        scope = CutSeparationScope::FullTree;
        return true;
    }
    return false;
}

std::string cut_management_stop_reason(const CutManagementStats& stats) {
    std::ostringstream reason;
    const char* separator = "";
    const auto append = [&](const char* name, bool active) {
        if (!active) {
            return;
        }
        reason << separator << name;
        separator = "+";
    };
    append("total-cap", stats.stopped_by_total_cap);
    append("time-fraction", stats.stopped_by_time_fraction);
    append("root-rounds", stats.stopped_by_root_rounds);
    append("root-cuts", stats.stopped_by_root_cuts);
    append("no-cut-streak", stats.stopped_by_no_cut_streak);
    append("low-violation", stats.stopped_by_low_violation);
    const std::string text = reason.str();
    return text.empty() ? std::string("none") : text;
}

}  // namespace spdp
