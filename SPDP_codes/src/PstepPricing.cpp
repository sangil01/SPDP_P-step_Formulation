#include "PstepPricing.h"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <map>
#include <optional>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace spdp {
namespace {

constexpr double kTolerance = 1e-9;

bool double_equal(double lhs, double rhs) {
    return std::fabs(lhs - rhs) <= kTolerance;
}

bool double_less_or_equal(double lhs, double rhs) {
    return lhs <= rhs + kTolerance;
}

std::string format_double(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(6) << value;
    return out.str();
}

State canonicalize_state(State state) {
    if (state[1] < state[0]) {
        std::swap(state[0], state[1]);
    }
    return state;
}

double lookup_dual(
    const std::map<NodeId, double>& duals,
    NodeId node_id
) {
    const auto found = duals.find(node_id);
    if (found == duals.end()) {
        return 0.0;
    }
    return found->second;
}

double lookup_dual(
    const std::map<NodeStateKey, double>& duals,
    const NodeStateKey& key
) {
    const auto found = duals.find(key);
    if (found == duals.end()) {
        return 0.0;
    }
    return found->second;
}

// symmetry-43 pickup ordering용 class:
// 같은 (pickup location, treatment location, type)의 pickup만 같은 class로 본다.
std::optional<PricingRequestClassKey> pickup_request_class_of_node(
    const MultiDiGraph& graph,
    NodeId node_id
) {
    if (!graph.is_physical_service_node(node_id)) {
        return std::nullopt;
    }

    const NodeSpec& node = graph.node(node_id);
    if (node.kind != NodeSpec::Kind::Pickup || !node.container_type.has_value() ||
        !node.landfill_location.has_value()) {
        return std::nullopt;
    }

    return PricingRequestClassKey{
        node.location,
        node.landfill_location.value(),
        node.container_type.value(),
    };
}

// delivery ordering symmetry용 class:
// delivery service는 E(h) -> N만 수행하므로 landfill identity는 쓰지 않고
// 같은 (delivery location, type)의 delivery만 같은 class로 본다.
std::optional<PricingDeliveryClassKey> delivery_request_class_of_node(
    const MultiDiGraph& graph,
    NodeId node_id
) {
    if (!graph.is_physical_service_node(node_id)) {
        return std::nullopt;
    }

    const NodeSpec& node = graph.node(node_id);
    if (node.kind != NodeSpec::Kind::Delivery || !node.container_type.has_value()) {
        return std::nullopt;
    }

    return PricingDeliveryClassKey{
        node.location,
        node.container_type.value(),
    };
}

std::size_t bit_word_count(std::size_t bit_count) {
    return (bit_count + 63U) / 64U;
}

struct Label {
    int node_state_index = -1;
    int depth = 0;
    double base_cost = 0.0;
    double total_time = 0.0;
    std::vector<std::uint64_t> visited_physical_mask;
    std::vector<int> pickup_max_request_idx_by_class;
    std::vector<int> delivery_max_request_idx_by_class;
    int parent_label_index = -1;
    int incoming_pricing_edge_index = -1;
};

struct ShallowSearchState {
    int node_state_index = -1;
    int depth = 0;
    double total_time = 0.0;
    std::vector<std::uint64_t> visited_physical_mask;
    std::vector<int> pickup_max_request_idx_by_class;
    std::vector<int> delivery_max_request_idx_by_class;
    std::vector<int> edge_ids;
};

bool mask_contains(
    const std::vector<std::uint64_t>& mask,
    int bit_index
) {
    if (bit_index < 0) {
        return false;
    }
    const std::size_t word = static_cast<std::size_t>(bit_index) / 64U;
    const std::size_t offset = static_cast<std::size_t>(bit_index) % 64U;
    return (mask[word] & (std::uint64_t{1} << offset)) != 0U;
}

void mask_add(
    std::vector<std::uint64_t>& mask,
    int bit_index
) {
    if (bit_index < 0) {
        return;
    }
    const std::size_t word = static_cast<std::size_t>(bit_index) / 64U;
    const std::size_t offset = static_cast<std::size_t>(bit_index) % 64U;
    mask[word] |= (std::uint64_t{1} << offset);
}

bool mask_subset_of(
    const std::vector<std::uint64_t>& lhs,
    const std::vector<std::uint64_t>& rhs
) {
    for (std::size_t idx = 0; idx < lhs.size(); ++idx) {
        if ((lhs[idx] & ~rhs[idx]) != 0U) {
            return false;
        }
    }
    return true;
}

bool symmetry_state_less_or_equal(
    const std::vector<int>& lhs,
    const std::vector<int>& rhs
) {
    for (std::size_t idx = 0; idx < lhs.size(); ++idx) {
        if (lhs[idx] > rhs[idx]) {
            return false;
        }
    }
    return true;
}

bool is_complete_label(
    NodeId start_node_id,
    int depth,
    int p
) {
    if (start_node_id == 0) {
        return depth >= 1 && depth <= p;
    }
    return depth == p;
}

double phase_edge_objective_cost(
    CGPhase phase,
    double original_edge_cost
) {
    if (phase == CGPhase::PhaseII) {
        return original_edge_cost;
    }
    return 0.0;
}

std::size_t compute_heuristic_search_limit(
    std::size_t output_limit,
    double search_ratio
) {
    const double scaled_limit = std::ceil(search_ratio * static_cast<double>(output_limit));
    const double capped_limit = std::min(
        scaled_limit,
        static_cast<double>(std::numeric_limits<std::size_t>::max())
    );
    return std::max<std::size_t>(
        output_limit,
        static_cast<std::size_t>(capped_limit)
    );
}

double compute_heuristic_start_score(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    HeuristicStartScoreMode mode,
    int start_node_state_index
) {
    switch (mode) {
        case HeuristicStartScoreMode::OneStepMin: {
            double best_score = std::numeric_limits<double>::infinity();
            for (int pricing_edge_index :
                 context.outgoing_edges_by_node_state[static_cast<std::size_t>(start_node_state_index)]) {
                const ForwardPricingContext::EdgeInfo& pricing_edge =
                    context.edges[static_cast<std::size_t>(pricing_edge_index)];
                const EdgeRecord& graph_edge =
                    graph.edges()[static_cast<std::size_t>(pricing_edge.graph_edge_id)];
                const double alpha_head =
                    graph.is_physical_service_node(graph_edge.v)
                        ? lookup_dual(dual_solution.visit_duals, graph_edge.v)
                        : 0.0;
                const double delta =
                    (pricing_edge.graph_edge_id >= 0 &&
                     static_cast<std::size_t>(pricing_edge.graph_edge_id) <
                         dual_solution.edge_duals.size())
                        ? dual_solution.edge_duals
                              [static_cast<std::size_t>(pricing_edge.graph_edge_id)]
                        : 0.0;
                const double score =
                    phase_edge_objective_cost(phase, graph_edge.data.cost) -
                    delta - 2.0 * alpha_head;
                best_score = std::min(best_score, score);
            }
            return best_score;
        }
    }
    return std::numeric_limits<double>::infinity();
}

std::vector<int> build_heuristic_start_order_impl(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    HeuristicStartScoreMode mode
) {
    std::vector<int> ordered_starts = context.start_node_state_indices;

    std::vector<std::pair<double, int>> scored_starts;
    scored_starts.reserve(ordered_starts.size());
    for (int start_node_state_index : ordered_starts) {
        scored_starts.push_back({
            compute_heuristic_start_score(
                graph,
                context,
                dual_solution,
                phase,
                mode,
                start_node_state_index
            ),
            start_node_state_index,
        });
    }

    std::stable_sort(
        scored_starts.begin(),
        scored_starts.end(),
        [](const std::pair<double, int>& lhs, const std::pair<double, int>& rhs) {
            if (double_equal(lhs.first, rhs.first)) {
                return lhs.second < rhs.second;
            }
            return lhs.first < rhs.first;
        }
    );

    ordered_starts.clear();
    ordered_starts.reserve(scored_starts.size());
    for (std::size_t idx = 0; idx < scored_starts.size(); ++idx) {
        ordered_starts.push_back(scored_starts[idx].second);
    }
    return ordered_starts;
}

std::vector<int> build_start_order(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    const ForwardPricingOptions& options
) {
    if (!options.explicit_start_node_state_indices.empty()) {
        return options.explicit_start_node_state_indices;
    }

    std::vector<int> ordered_starts = context.start_node_state_indices;
    if (!options.heuristic_pricing) {
        return ordered_starts;
    }

    ordered_starts = build_heuristic_start_order_impl(
        graph,
        context,
        dual_solution,
        phase,
        options.heuristic_start_score_mode
    );

    const std::size_t ratio_limit = static_cast<std::size_t>(
        std::ceil(options.heuristic_start_ratio * static_cast<double>(ordered_starts.size()))
    );
    const std::size_t limit = std::min(
        options.heuristic_max_starts,
        std::min(ratio_limit, ordered_starts.size())
    );
    ordered_starts.resize(limit);
    return ordered_starts;
}

double choose_best_tau(
    const MultiDiGraph& graph,
    const CGDualSolution& dual_solution,
    NodeId start_node_id,
    const State& start_state,
    NodeId last_node_id,
    const State& last_state,
    double total_time,
    double time_limit
) {
    const double latest_start = time_limit - total_time;
    if (latest_start < -kTolerance) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (start_node_id == 0) {
        return 0.0;
    }

    if (last_node_id == graph.end_node_id()) {
        return std::max(0.0, latest_start);
    }

    const NodeStateKey start_key{start_node_id, canonicalize_state(start_state)};
    const NodeStateKey last_key{last_node_id, canonicalize_state(last_state)};
    const double gamma_start = graph.is_physical_service_node(start_node_id)
                                   ? lookup_dual(dual_solution.time_duals, start_key)
                                   : 0.0;
    const double gamma_end = graph.is_physical_service_node(last_node_id)
                                 ? lookup_dual(dual_solution.time_duals, last_key)
                                 : 0.0;

    return (gamma_end - gamma_start < 0.0) ? std::max(0.0, latest_start) : 0.0;
}

std::string build_column_key(
    const std::vector<int>& edge_ids,
    double tau
) {
    // edge id sequence와 tau를 조합한 문자열 key. 서로 다른 p-step임을 파악하기 위해 사용.
    std::ostringstream out;
    for (std::size_t idx = 0; idx < edge_ids.size(); ++idx) {
        if (idx > 0U) {
            out << ',';
        }
        out << edge_ids[idx];
    }
    out << '|';
    out << std::fixed << std::setprecision(9) << tau;
    return out.str();
}

// 일단은 original reduced cost form을 사용 (누적 cost + terminal dependent cost로 계산 x)
double evaluate_sparse_column_reduced_cost(
    const CGDualSolution& dual_solution,
    CGPhase phase,
    double total_cost,
    const std::vector<std::pair<NodeId, int>>& visit_coefficients,
    const std::vector<std::pair<NodeStateKey, int>>& state_coefficients,
    const std::vector<std::pair<NodeStateKey, double>>& time_coefficients,
    const std::vector<int>& edge_incidence
) {
    double reduced_cost = (phase == CGPhase::PhaseII) ? total_cost : 0.0;

    for (const auto& entry : visit_coefficients) {
        reduced_cost -= static_cast<double>(entry.second) *
                        lookup_dual(dual_solution.visit_duals, entry.first);
    }
    for (const auto& entry : state_coefficients) {
        reduced_cost -= static_cast<double>(entry.second) *
                        lookup_dual(dual_solution.state_duals, entry.first);
    }
    for (const auto& entry : time_coefficients) {
        reduced_cost -= entry.second *
                        lookup_dual(dual_solution.time_duals, entry.first);
    }
    for (int edge_id : edge_incidence) {
        if (edge_id >= 0 &&
            static_cast<std::size_t>(edge_id) < dual_solution.edge_duals.size()) {
            reduced_cost -= dual_solution.edge_duals[static_cast<std::size_t>(edge_id)];
        }
    }

    return reduced_cost;
}

double compute_edge_score(
    const MultiDiGraph& graph,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    HeuristicStartScoreMode mode,
    int graph_edge_id
) {
    switch (mode) {
        case HeuristicStartScoreMode::OneStepMin: {
            const EdgeRecord& graph_edge =
                graph.edges()[static_cast<std::size_t>(graph_edge_id)];
            const double alpha_head =
                graph.is_physical_service_node(graph_edge.v)
                    ? lookup_dual(dual_solution.visit_duals, graph_edge.v)
                    : 0.0;
            const double delta =
                (graph_edge_id >= 0 &&
                 static_cast<std::size_t>(graph_edge_id) < dual_solution.edge_duals.size())
                    ? dual_solution.edge_duals[static_cast<std::size_t>(graph_edge_id)]
                    : 0.0;
            return phase_edge_objective_cost(phase, graph_edge.data.cost) -
                   delta - 2.0 * alpha_head;
        }
    }
    return std::numeric_limits<double>::infinity();
}

CGColumn build_shallow_column_from_graph_edges(
    const MultiDiGraph& graph,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    double time_limit,
    const std::vector<int>& edge_ids,
    double reduced_cost_tolerance,
    double tau
) {
    CGColumn column;
    if (edge_ids.empty()) {
        return column;
    }

    column.q = static_cast<int>(edge_ids.size());
    column.tau = tau;

    const EdgeRecord& first_edge = graph.edges()[static_cast<std::size_t>(edge_ids.front())];
    column.start_node_id = first_edge.u;
    column.start_state = canonicalize_state(first_edge.data.start_state);
    column.node_sequence.push_back(first_edge.u);
    column.state_sequence.push_back(column.start_state);

    for (int edge_id : edge_ids) {
        const EdgeRecord& edge = graph.edges()[static_cast<std::size_t>(edge_id)];
        column.edge_ids.push_back(edge_id);
        column.edge_incidence.push_back(edge_id);
        column.total_time += edge.data.time;
        column.total_cost += edge.data.cost;
        column.node_sequence.push_back(edge.v);
        column.state_sequence.push_back(canonicalize_state(edge.data.end_state));
    }

    column.last_node_id = column.node_sequence.back();
    column.last_state = canonicalize_state(column.state_sequence.back());

    if (!std::isfinite(
            choose_best_tau(
                graph,
                dual_solution,
                column.start_node_id,
                column.start_state,
                column.last_node_id,
                column.last_state,
                column.total_time,
                time_limit
            )
        )) {
        return CGColumn{};
    }

    for (std::size_t pos = 0; pos < column.node_sequence.size(); ++pos) {
        const NodeId node_id = column.node_sequence[pos];
        if (!graph.is_physical_service_node(node_id)) {
            continue;
        }
        const int coefficient =
            (pos == 0U || pos + 1U == column.node_sequence.size()) ? 1 : 2;
        column.visit_coefficients.push_back({node_id, coefficient});
    }

    if (graph.is_physical_service_node(column.start_node_id)) {
        column.state_coefficients.push_back(
            {NodeStateKey{column.start_node_id, canonicalize_state(column.start_state)}, 1}
        );
        column.time_coefficients.push_back(
            {NodeStateKey{column.start_node_id, canonicalize_state(column.start_state)}, tau}
        );
    }
    if (graph.is_physical_service_node(column.last_node_id)) {
        column.state_coefficients.push_back(
            {NodeStateKey{column.last_node_id, canonicalize_state(column.last_state)}, -1}
        );
        column.time_coefficients.push_back(
            {NodeStateKey{column.last_node_id, canonicalize_state(column.last_state)},
             -(tau + column.total_time)}
        );
    }

    column.reduced_cost = evaluate_sparse_column_reduced_cost(
        dual_solution,
        phase,
        column.total_cost,
        column.visit_coefficients,
        column.state_coefficients,
        column.time_coefficients,
        column.edge_incidence
    );
    if (!(column.reduced_cost < reduced_cost_tolerance)) {
        return CGColumn{};
    }
    return column;
}

void append_shallow_column_variants(
    const MultiDiGraph& graph,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    double time_limit,
    const std::vector<int>& edge_ids,
    double reduced_cost_tolerance,
    std::unordered_map<std::string, std::size_t>& seen_columns,
    std::vector<CGColumn>& columns
) {
    if (edge_ids.empty()) {
        return;
    }

    const EdgeRecord& first_edge = graph.edges()[static_cast<std::size_t>(edge_ids.front())];
    const EdgeRecord& last_edge = graph.edges()[static_cast<std::size_t>(edge_ids.back())];

    double total_time = 0.0;
    for (int edge_id : edge_ids) {
        total_time += graph.edges()[static_cast<std::size_t>(edge_id)].data.time;
    }

    const double latest_start = time_limit - total_time;
    if (latest_start < -kTolerance) {
        return;
    }

    std::vector<double> tau_candidates;
    if (first_edge.u == 0) {
        tau_candidates.push_back(0.0);
    } else if (last_edge.v == graph.end_node_id()) {
        tau_candidates.push_back(std::max(0.0, latest_start));
    } else {
        tau_candidates.push_back(0.0);
        const double latest = std::max(0.0, latest_start);
        if (!double_equal(latest, 0.0)) {
            tau_candidates.push_back(latest);
        }
    }

    for (double tau : tau_candidates) {
        const std::string key = build_column_key(edge_ids, tau);
        if (seen_columns.find(key) != seen_columns.end()) {
            continue;
        }
        CGColumn column = build_shallow_column_from_graph_edges(
            graph,
            dual_solution,
            phase,
            time_limit,
            edge_ids,
            reduced_cost_tolerance,
            tau
        );
        if (column.edge_ids.empty()) {
            continue;
        }
        columns.push_back(std::move(column));
        seen_columns.emplace(key, columns.size() - 1U);
    }
}

bool passes_forward_extension_filters(
    const ForwardPricingContext& context,
    const ForwardPricingContext::NodeStateInfo& next_info,
    const ShallowSearchState& current_state,
    const ForwardPricingOptions& options
) {
    if (next_info.is_physical_service_node &&
        mask_contains(current_state.visited_physical_mask, next_info.physical_bit_index)) {
        return false;
    }

    if (options.prune_pickup_symmetry_43 &&
        next_info.is_pickup &&
        next_info.pickup_request_class_index >= 0) {
        const int current_max =
            current_state.pickup_max_request_idx_by_class
                [static_cast<std::size_t>(next_info.pickup_request_class_index)];
        if (next_info.request_index <= current_max) {
            return false;
        }
    }

    if (options.prune_delivery_symmetry_43 &&
        next_info.is_delivery &&
        next_info.delivery_request_class_index >= 0) {
        const int current_max =
            current_state.delivery_max_request_idx_by_class
                [static_cast<std::size_t>(next_info.delivery_request_class_index)];
        if (next_info.request_index <= current_max) {
            return false;
        }
    }

    return true;
}

ShallowSearchState extend_shallow_state(
    const ForwardPricingContext& context,
    const ForwardPricingContext::EdgeInfo& pricing_edge,
    const ShallowSearchState& current_state,
    const ForwardPricingOptions& options
) {
    ShallowSearchState next_state = current_state;
    next_state.node_state_index = pricing_edge.to_node_state_index;
    next_state.depth = current_state.depth + 1;
    next_state.total_time = current_state.total_time + pricing_edge.time;
    next_state.edge_ids.push_back(pricing_edge.graph_edge_id);

    const auto& next_info =
        context.node_states[static_cast<std::size_t>(pricing_edge.to_node_state_index)];
    if (next_info.is_physical_service_node) {
        mask_add(next_state.visited_physical_mask, next_info.physical_bit_index);
    }
    if (options.prune_pickup_symmetry_43 &&
        next_info.is_pickup &&
        next_info.pickup_request_class_index >= 0) {
        next_state.pickup_max_request_idx_by_class
            [static_cast<std::size_t>(next_info.pickup_request_class_index)] =
                next_info.request_index;
    }
    if (options.prune_delivery_symmetry_43 &&
        next_info.is_delivery &&
        next_info.delivery_request_class_index >= 0) {
        next_state.delivery_max_request_idx_by_class
            [static_cast<std::size_t>(next_info.delivery_request_class_index)] =
                next_info.request_index;
    }
    return next_state;
}

std::vector<int> collect_feasible_shallow_edges(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    const ForwardPricingOptions& options,
    const ShallowSearchState& current_state,
    std::size_t limit
) {
    std::vector<std::pair<double, int>> scored_edges;
    for (int pricing_edge_index :
         context.outgoing_edges_by_node_state[static_cast<std::size_t>(current_state.node_state_index)]) {
        const auto& pricing_edge = context.edges[static_cast<std::size_t>(pricing_edge_index)];
        if (!double_less_or_equal(
                current_state.total_time + pricing_edge.time,
                context.time_limit
            )) {
            continue;
        }
        const auto& next_info =
            context.node_states[static_cast<std::size_t>(pricing_edge.to_node_state_index)];
        if (!passes_forward_extension_filters(context, next_info, current_state, options)) {
            continue;
        }

        scored_edges.push_back({
            compute_edge_score(
                graph,
                dual_solution,
                phase,
                options.heuristic_start_score_mode,
                pricing_edge.graph_edge_id
            ),
            pricing_edge_index,
        });
    }

    if (limit > 0U && scored_edges.size() > limit) {
        std::partial_sort(
            scored_edges.begin(),
            scored_edges.begin() + static_cast<std::ptrdiff_t>(limit),
            scored_edges.end(),
            [](const std::pair<double, int>& lhs, const std::pair<double, int>& rhs) {
                if (double_equal(lhs.first, rhs.first)) {
                    return lhs.second < rhs.second;
                }
                return lhs.first < rhs.first;
            }
        );
        scored_edges.resize(limit);
    } else {
        std::sort(
            scored_edges.begin(),
            scored_edges.end(),
            [](const std::pair<double, int>& lhs, const std::pair<double, int>& rhs) {
                if (double_equal(lhs.first, rhs.first)) {
                    return lhs.second < rhs.second;
                }
                return lhs.first < rhs.first;
            }
        );
    }

    std::vector<int> feasible_edges;
    feasible_edges.reserve(scored_edges.size());
    for (const auto& entry : scored_edges) {
        feasible_edges.push_back(entry.second);
    }
    return feasible_edges;
}

std::vector<int> collect_top_k_forward_labeling_edges(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    const ForwardPricingOptions& options,
    const Label& current_label,
    std::size_t top_k_next
) {
    std::vector<std::pair<double, int>> scored_edges;
    scored_edges.reserve(
        context.outgoing_edges_by_node_state[static_cast<std::size_t>(current_label.node_state_index)]
            .size()
    );

    for (int pricing_edge_index :
         context.outgoing_edges_by_node_state[static_cast<std::size_t>(current_label.node_state_index)]) {
        const ForwardPricingContext::EdgeInfo& pricing_edge =
            context.edges[static_cast<std::size_t>(pricing_edge_index)];
        const ForwardPricingContext::NodeStateInfo& next_info =
            context.node_states[static_cast<std::size_t>(pricing_edge.to_node_state_index)];

        if (!double_less_or_equal(
                current_label.total_time + pricing_edge.time,
                context.time_limit
            )) {
            continue;
        }

        if (next_info.is_physical_service_node &&
            mask_contains(current_label.visited_physical_mask, next_info.physical_bit_index)) {
            continue;
        }

        if (options.prune_pickup_symmetry_43 &&
            next_info.is_pickup &&
            next_info.pickup_request_class_index >= 0) {
            const int current_max =
                current_label.pickup_max_request_idx_by_class
                    [static_cast<std::size_t>(next_info.pickup_request_class_index)];
            if (next_info.request_index <= current_max) {
                continue;
            }
        }

        if (options.prune_delivery_symmetry_43 &&
            next_info.is_delivery &&
            next_info.delivery_request_class_index >= 0) {
            const int current_max =
                current_label.delivery_max_request_idx_by_class
                    [static_cast<std::size_t>(next_info.delivery_request_class_index)];
            if (next_info.request_index <= current_max) {
                continue;
            }
        }

        scored_edges.push_back({
            compute_edge_score(
                graph,
                dual_solution,
                phase,
                options.heuristic_start_score_mode,
                pricing_edge.graph_edge_id
            ),
            pricing_edge_index,
        });
    }

    if (top_k_next > 0U && scored_edges.size() > top_k_next) {
        std::partial_sort(
            scored_edges.begin(),
            scored_edges.begin() + static_cast<std::ptrdiff_t>(top_k_next),
            scored_edges.end(),
            [](const std::pair<double, int>& lhs, const std::pair<double, int>& rhs) {
                if (double_equal(lhs.first, rhs.first)) {
                    return lhs.second < rhs.second;
                }
                return lhs.first < rhs.first;
            }
        );
        scored_edges.resize(top_k_next);
    } else {
        std::sort(
            scored_edges.begin(),
            scored_edges.end(),
            [](const std::pair<double, int>& lhs, const std::pair<double, int>& rhs) {
                if (double_equal(lhs.first, rhs.first)) {
                    return lhs.second < rhs.second;
                }
                return lhs.first < rhs.first;
            }
        );
    }

    std::vector<int> pricing_edge_indices;
    pricing_edge_indices.reserve(scored_edges.size());
    for (const auto& entry : scored_edges) {
        pricing_edge_indices.push_back(entry.second);
    }
    return pricing_edge_indices;
}

bool evaluate_shallow_complete_state(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    const ForwardPricingOptions& options,
    const ShallowSearchState& complete_state,
    std::unordered_map<std::string, std::size_t>& seen_columns_for_start,
    std::vector<CGColumn>& generated_columns_for_start,
    std::vector<CGColumn>& accepted_columns_for_start,
    std::size_t search_max_columns_per_start,
    ForwardPricingResult& result
) {
    ++result.complete_label_count;
    const std::size_t previous_generated_count = generated_columns_for_start.size();
    append_shallow_column_variants(
        graph,
        dual_solution,
        phase,
        context.time_limit,
        complete_state.edge_ids,
        options.reduced_cost_tolerance,
        seen_columns_for_start,
        generated_columns_for_start
    );

    for (std::size_t idx = previous_generated_count;
         idx < generated_columns_for_start.size();
         ++idx) {
        const CGColumn& column = generated_columns_for_start[idx];
        result.best_reduced_cost = std::min(result.best_reduced_cost, column.reduced_cost);
        accepted_columns_for_start.push_back(column);
        if (accepted_columns_for_start.size() >= search_max_columns_per_start) {
            return true;
        }
    }

    return false;
}

CGColumn build_column_from_complete_label(
    const MultiDiGraph& graph,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    double time_limit,
    const ForwardPricingContext& context,
    const std::vector<Label>& labels,
    int label_index,
    double reduced_cost_tolerance
) {
    const Label& complete_label = labels[static_cast<std::size_t>(label_index)];

    std::vector<int> reversed_edge_ids;
    std::vector<int> reversed_node_state_indices;
    int cursor = label_index;
    while (cursor >= 0) {
        const Label& label = labels[static_cast<std::size_t>(cursor)];
        reversed_node_state_indices.push_back(label.node_state_index);
        if (label.incoming_pricing_edge_index >= 0) {
            reversed_edge_ids.push_back(label.incoming_pricing_edge_index);
        }
        cursor = label.parent_label_index;
    }

    std::reverse(reversed_node_state_indices.begin(), reversed_node_state_indices.end());
    std::reverse(reversed_edge_ids.begin(), reversed_edge_ids.end());

    CGColumn column;
    column.q = complete_label.depth;
    column.total_time = complete_label.total_time;

    for (int node_state_index : reversed_node_state_indices) {
        const ForwardPricingContext::NodeStateInfo& info =
            context.node_states[static_cast<std::size_t>(node_state_index)];
        column.node_sequence.push_back(info.node_id);
        column.state_sequence.push_back(canonicalize_state(info.state));
    }

    for (int pricing_edge_index : reversed_edge_ids) {
        const ForwardPricingContext::EdgeInfo& pricing_edge =
            context.edges[static_cast<std::size_t>(pricing_edge_index)];
        const EdgeRecord& edge = graph.edges()[static_cast<std::size_t>(pricing_edge.graph_edge_id)];
        column.edge_ids.push_back(pricing_edge.graph_edge_id);
        column.edge_incidence.push_back(pricing_edge.graph_edge_id);
        column.total_cost += edge.data.cost;
    }

    column.start_node_id = column.node_sequence.front();
    column.last_node_id = column.node_sequence.back();
    column.start_state = canonicalize_state(column.state_sequence.front());
    column.last_state = canonicalize_state(column.state_sequence.back());

    // dynamic master coefficients independent of tau.
    for (std::size_t pos = 0; pos < column.node_sequence.size(); ++pos) {
        const NodeId node_id = column.node_sequence[pos];
        if (!graph.is_physical_service_node(node_id)) {
            continue;
        }
        const int coefficient =
            (pos == 0U || pos + 1U == column.node_sequence.size()) ? 1 : 2;
        column.visit_coefficients.push_back({node_id, coefficient});
    }

    if (graph.is_physical_service_node(column.start_node_id)) {
        column.state_coefficients.push_back(
            {NodeStateKey{column.start_node_id, canonicalize_state(column.start_state)}, 1}
        );
    }
    if (graph.is_physical_service_node(column.last_node_id)) {
        column.state_coefficients.push_back(
            {NodeStateKey{column.last_node_id, canonicalize_state(column.last_state)}, -1}
        );
    }

    column.tau = choose_best_tau(
        graph,
        dual_solution,
        column.start_node_id,
        column.start_state,
        column.last_node_id,
        column.last_state,
        column.total_time,
        time_limit
    );
    if (!std::isfinite(column.tau)) {
        column.edge_ids.clear();
        column.node_sequence.clear();
        column.state_sequence.clear();
        column.visit_coefficients.clear();
        column.state_coefficients.clear();
        column.time_coefficients.clear();
        column.edge_incidence.clear();
        return column;
    }

    if (graph.is_physical_service_node(column.start_node_id)) {
        column.time_coefficients.push_back(
            {NodeStateKey{column.start_node_id, canonicalize_state(column.start_state)}, column.tau}
        );
    }
    if (graph.is_physical_service_node(column.last_node_id)) {
        column.time_coefficients.push_back(
            {NodeStateKey{column.last_node_id, canonicalize_state(column.last_state)},
             -(column.tau + column.total_time)}
        );
    }

    column.reduced_cost = evaluate_sparse_column_reduced_cost(
        dual_solution,
        phase,
        column.total_cost,
        column.visit_coefficients,
        column.state_coefficients,
        column.time_coefficients,
        column.edge_incidence
    );

    if (!(column.reduced_cost < reduced_cost_tolerance)) {
        column.edge_ids.clear();
        column.node_sequence.clear();
        column.state_sequence.clear();
        column.visit_coefficients.clear();
        column.state_coefficients.clear();
        column.time_coefficients.clear();
        column.edge_incidence.clear();
    }

    return column;
}

bool dominates_label(
    const Label& lhs,
    const Label& rhs
) {
    return double_less_or_equal(lhs.base_cost, rhs.base_cost) &&
           double_less_or_equal(lhs.total_time, rhs.total_time) &&
           mask_subset_of(lhs.visited_physical_mask, rhs.visited_physical_mask) &&
           symmetry_state_less_or_equal(
               lhs.pickup_max_request_idx_by_class,
               rhs.pickup_max_request_idx_by_class
           ) &&
           symmetry_state_less_or_equal(
               lhs.delivery_max_request_idx_by_class,
               rhs.delivery_max_request_idx_by_class
           );
}

std::size_t count_surviving_labels(
    const std::vector<std::vector<std::vector<int>>>& bucket_label_indices
) {
    std::size_t count = 0;
    for (const auto& depth_buckets : bucket_label_indices) {
        for (const auto& bucket : depth_buckets) {
            count += bucket.size();
        }
    }
    return count;
}

std::vector<State> collect_sigma_values_for_node(
    const MultiDiGraph& graph,
    NodeId node_id
) {
    std::set<State> sigma_set;

    for (const EdgeRecord& edge : graph.edges()) {
        if (edge.u == node_id && graph.is_physical_service_node(edge.u)) {
            sigma_set.insert(canonicalize_state(edge.data.start_state));
        }
        if (edge.v == node_id && graph.is_physical_service_node(edge.v)) {
            sigma_set.insert(canonicalize_state(edge.data.end_state));
        }
    }

    return std::vector<State>(sigma_set.begin(), sigma_set.end());
}

}  // namespace

ForwardPricingContext build_forward_pricing_context(
    const MultiDiGraph& graph,
    int p,
    double time_limit
) {
    if (p < 1) {
        throw std::runtime_error("Forward pricing requires p >= 1.");
    }
    if (time_limit < 0.0) {
        throw std::runtime_error("Forward pricing requires a nonnegative time limit.");
    }

    ForwardPricingContext context;
    context.p = p;
    context.time_limit = time_limit;
    context.end_node_id = graph.end_node_id();

    std::map<NodeId, int> physical_bit_by_node;
    for (NodeId node_id = 1; node_id < context.end_node_id; ++node_id) {
        if (!graph.is_physical_service_node(node_id)) {
            continue;
        }
        const int next_bit = static_cast<int>(physical_bit_by_node.size());
        physical_bit_by_node.emplace(node_id, next_bit);
        context.sigma_by_node[node_id] = collect_sigma_values_for_node(graph, node_id);
    }
    context.physical_node_count = physical_bit_by_node.size();

    std::map<PricingRequestClassKey, int> pickup_request_class_index_by_key;
    std::map<PricingDeliveryClassKey, int> delivery_request_class_index_by_key;
    for (NodeId node_id = 1; node_id < context.end_node_id; ++node_id) {
        const std::optional<PricingRequestClassKey> pickup_key =
            pickup_request_class_of_node(graph, node_id);
        if (pickup_key.has_value() &&
            pickup_request_class_index_by_key.find(pickup_key.value()) ==
                pickup_request_class_index_by_key.end()) {
            const int next_index = static_cast<int>(pickup_request_class_index_by_key.size());
            pickup_request_class_index_by_key.emplace(pickup_key.value(), next_index);
        }

        const std::optional<PricingDeliveryClassKey> delivery_key =
            delivery_request_class_of_node(graph, node_id);
        if (delivery_key.has_value() &&
            delivery_request_class_index_by_key.find(delivery_key.value()) ==
                delivery_request_class_index_by_key.end()) {
            const int next_index = static_cast<int>(delivery_request_class_index_by_key.size());
            delivery_request_class_index_by_key.emplace(delivery_key.value(), next_index);
        }
    }
    context.pickup_request_class_count = pickup_request_class_index_by_key.size();
    context.delivery_request_class_count = delivery_request_class_index_by_key.size();

    for (const EdgeRecord& edge : graph.edges()) {
        const NodeStateKey start_key{edge.u, canonicalize_state(edge.data.start_state)};
        const NodeStateKey end_key{edge.v, canonicalize_state(edge.data.end_state)};

        if (context.node_state_index_by_key.find(start_key) == context.node_state_index_by_key.end()) {
            ForwardPricingContext::NodeStateInfo info;
            info.node_id = start_key.node_id;
            info.state = start_key.state;
            info.is_physical_service_node = graph.is_physical_service_node(start_key.node_id);
            const auto physical_found = physical_bit_by_node.find(start_key.node_id);
            if (physical_found != physical_bit_by_node.end()) {
                info.physical_bit_index = physical_found->second;
            }

            const NodeSpec& node = graph.node(start_key.node_id);
            if (node.kind == NodeSpec::Kind::Pickup) {
                info.is_pickup = true;
                info.request_index = node.request_idx.value();
                const std::optional<PricingRequestClassKey> class_key =
                    pickup_request_class_of_node(graph, start_key.node_id);
                if (class_key.has_value()) {
                    info.pickup_request_class_index =
                        pickup_request_class_index_by_key[class_key.value()];
                }
            } else if (node.kind == NodeSpec::Kind::Delivery) {
                info.is_delivery = true;
                info.request_index = node.request_idx.value();
                const std::optional<PricingDeliveryClassKey> class_key =
                    delivery_request_class_of_node(graph, start_key.node_id);
                if (class_key.has_value()) {
                    info.delivery_request_class_index =
                        delivery_request_class_index_by_key[class_key.value()];
                }
            }

            const int new_index = static_cast<int>(context.node_states.size());
            context.node_state_index_by_key.emplace(start_key, new_index);
            context.node_states.push_back(info);
        }

        if (context.node_state_index_by_key.find(end_key) == context.node_state_index_by_key.end()) {
            ForwardPricingContext::NodeStateInfo info;
            info.node_id = end_key.node_id;
            info.state = end_key.state;
            info.is_physical_service_node = graph.is_physical_service_node(end_key.node_id);
            const auto physical_found = physical_bit_by_node.find(end_key.node_id);
            if (physical_found != physical_bit_by_node.end()) {
                info.physical_bit_index = physical_found->second;
            }

            const NodeSpec& node = graph.node(end_key.node_id);
            if (node.kind == NodeSpec::Kind::Pickup) {
                info.is_pickup = true;
                info.request_index = node.request_idx.value();
                const std::optional<PricingRequestClassKey> class_key =
                    pickup_request_class_of_node(graph, end_key.node_id);
                if (class_key.has_value()) {
                    info.pickup_request_class_index =
                        pickup_request_class_index_by_key[class_key.value()];
                }
            } else if (node.kind == NodeSpec::Kind::Delivery) {
                info.is_delivery = true;
                info.request_index = node.request_idx.value();
                const std::optional<PricingDeliveryClassKey> class_key =
                    delivery_request_class_of_node(graph, end_key.node_id);
                if (class_key.has_value()) {
                    info.delivery_request_class_index =
                        delivery_request_class_index_by_key[class_key.value()];
                }
            }

            const int new_index = static_cast<int>(context.node_states.size());
            context.node_state_index_by_key.emplace(end_key, new_index);
            context.node_states.push_back(info);
        }
    }

    context.outgoing_edges_by_node_state.assign(context.node_states.size(), {});
    context.incoming_edges_by_node_state.assign(context.node_states.size(), {});

    for (std::size_t edge_idx = 0; edge_idx < graph.edges().size(); ++edge_idx) {
        const EdgeRecord& edge = graph.edges()[edge_idx];
        const NodeStateKey start_key{edge.u, canonicalize_state(edge.data.start_state)};
        const NodeStateKey end_key{edge.v, canonicalize_state(edge.data.end_state)};
        const int from_index = context.node_state_index_by_key.at(start_key);
        const int to_index = context.node_state_index_by_key.at(end_key);

        ForwardPricingContext::EdgeInfo info;
        info.graph_edge_id = static_cast<int>(edge_idx);
        info.from_node_state_index = from_index;
        info.to_node_state_index = to_index;
        info.time = edge.data.time;
        info.cost = edge.data.cost;

        const int pricing_edge_index = static_cast<int>(context.edges.size());
        context.edges.push_back(info);
        context.outgoing_edges_by_node_state[static_cast<std::size_t>(from_index)].push_back(
            pricing_edge_index
        );
        context.incoming_edges_by_node_state[static_cast<std::size_t>(to_index)].push_back(
            pricing_edge_index
        );
    }

    for (std::size_t node_state_index = 0; node_state_index < context.node_states.size();
         ++node_state_index) {
        const NodeId node_id = context.node_states[node_state_index].node_id;
        if (node_id == context.end_node_id) {
        } else if (!context.outgoing_edges_by_node_state[node_state_index].empty()) {
            context.start_node_state_indices.push_back(static_cast<int>(node_state_index));
        }

        if (node_id == 0) {
        } else if (!context.incoming_edges_by_node_state[node_state_index].empty()) {
            context.end_node_state_indices.push_back(static_cast<int>(node_state_index));
        }
    }

    return context;
}

ForwardPricingResult run_forward_pricing(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    const ForwardPricingOptions& options
) {
    ForwardPricingResult result;
    result.best_reduced_cost = std::numeric_limits<double>::infinity();
    const std::size_t search_max_columns_per_start =
        options.heuristic_pricing
            ? compute_heuristic_search_limit(
                  options.max_columns_per_start,
                  options.heuristic_search_column_ratio
              )
            : options.max_columns_per_start;
    const std::size_t search_max_total_columns =
        options.heuristic_pricing
            ? compute_heuristic_search_limit(
                  options.max_total_columns,
                  options.heuristic_search_column_ratio
              )
            : options.max_total_columns;

    const std::size_t mask_word_count = bit_word_count(context.physical_node_count); // physical node 방문 여부를 나타내는 set을 비트마스크로 표현할 때 필요한 64bit cell 수.
    std::vector<CGColumn> candidate_columns;
    std::vector<CGColumn> deferred_candidate_columns;
    std::size_t total_negative_column_count = 0U;
    const std::vector<int> ordered_start_node_state_indices =
        build_start_order(graph, context, dual_solution, phase, options);

    for (int start_node_state_index : ordered_start_node_state_indices) {
        if (options.heuristic_pricing &&
            total_negative_column_count >= search_max_total_columns) {
            break;
        }

        const ForwardPricingContext::NodeStateInfo& start_info =
            context.node_states[static_cast<std::size_t>(start_node_state_index)];

        std::vector<Label> labels;
        labels.reserve(128);

        Label root_label;
        root_label.node_state_index = start_node_state_index;
        root_label.depth = 0;
        root_label.base_cost = 0.0;
        root_label.total_time = 0.0;
        root_label.visited_physical_mask.assign(mask_word_count, 0U);
        root_label.pickup_max_request_idx_by_class.assign(
            context.pickup_request_class_count,
            -1
        );
        root_label.delivery_max_request_idx_by_class.assign(
            context.delivery_request_class_count,
            -1
        );
        if (start_info.is_physical_service_node) {
            mask_add(root_label.visited_physical_mask, start_info.physical_bit_index);
        }
        if (options.prune_pickup_symmetry_43 &&
            start_info.is_pickup &&
            start_info.pickup_request_class_index >= 0) {
            root_label.pickup_max_request_idx_by_class
                [static_cast<std::size_t>(start_info.pickup_request_class_index)] =
                    start_info.request_index;
        }
        if (options.prune_delivery_symmetry_43 &&
            start_info.is_delivery &&
            start_info.delivery_request_class_index >= 0) {
            root_label.delivery_max_request_idx_by_class
                [static_cast<std::size_t>(start_info.delivery_request_class_index)] =
                    start_info.request_index;
        }
        labels.push_back(root_label);
        ++result.start_label_count;

        // same (x, depth) bucket에서 dominance를 본다.
        std::vector<std::vector<std::vector<int>>> bucket_label_indices(
            context.p + 1,
            std::vector<std::vector<int>>(context.node_states.size())
        );
        bucket_label_indices[0][static_cast<std::size_t>(start_node_state_index)].push_back(0);

        std::vector<CGColumn> best_columns_for_start;
        std::vector<int> current_layer_indices{0};
        bool stop_current_start = false;

        for (int depth = 0; depth <= context.p && !current_layer_indices.empty(); ++depth) {
            for (int label_index : current_layer_indices) {
                const Label current_label = labels[static_cast<std::size_t>(label_index)];

                if (is_complete_label(start_info.node_id, current_label.depth, context.p)) {
                    ++result.complete_label_count;
                    CGColumn column = build_column_from_complete_label(
                        graph,
                        dual_solution,
                        phase,
                        context.time_limit,
                        context,
                        labels,
                        label_index,
                        options.reduced_cost_tolerance
                    );
                    if (!column.edge_ids.empty()) {
                        result.best_reduced_cost = std::min(result.best_reduced_cost, column.reduced_cost);
                        best_columns_for_start.push_back(std::move(column));
                        ++total_negative_column_count;
                        if (options.heuristic_pricing &&
                            (best_columns_for_start.size() >= search_max_columns_per_start ||
                             total_negative_column_count >= search_max_total_columns)) {
                            stop_current_start = true;
                            break;
                        }
                    }
                }

                // depot에서는 complete label이 depth 1 이상인 경우이므로 complete 평가 후에도 확장 가능할 수 있다.
                if (current_label.depth >= context.p) {
                    continue;
                }

                const bool use_top_k_next =
                    options.heuristic_pricing && options.labeling_top_k_next > 0U;
                std::vector<int> top_k_pricing_edge_indices;
                const std::vector<int>* pricing_edge_indices =
                    &context.outgoing_edges_by_node_state
                        [static_cast<std::size_t>(current_label.node_state_index)];
                if (use_top_k_next) {
                    top_k_pricing_edge_indices = collect_top_k_forward_labeling_edges(
                        graph,
                        context,
                        dual_solution,
                        phase,
                        options,
                        current_label,
                        options.labeling_top_k_next
                    );
                    pricing_edge_indices = &top_k_pricing_edge_indices;
                }

                for (int pricing_edge_index : *pricing_edge_indices) {
                    const ForwardPricingContext::EdgeInfo& pricing_edge =
                        context.edges[static_cast<std::size_t>(pricing_edge_index)];
                    const ForwardPricingContext::NodeStateInfo& next_info =
                        context.node_states[static_cast<std::size_t>(pricing_edge.to_node_state_index)];

                    if (!use_top_k_next) {
                        if (!double_less_or_equal(
                                current_label.total_time + pricing_edge.time,
                                context.time_limit
                            )) {
                            continue;
                        }

                        if (next_info.is_physical_service_node &&
                            mask_contains(
                                current_label.visited_physical_mask,
                                next_info.physical_bit_index
                            )) {
                            continue;
                        }

                        if (options.prune_pickup_symmetry_43 &&
                            next_info.is_pickup &&
                            next_info.pickup_request_class_index >= 0) {
                            const int current_max =
                                current_label.pickup_max_request_idx_by_class
                                    [static_cast<std::size_t>(next_info.pickup_request_class_index)];
                            if (next_info.request_index <= current_max) {
                                continue;
                            }
                        }

                        // 같은 (delivery location, type) class의 delivery들은
                        // request_idx 오름차순으로만 등장하도록 제한한다.
                        if (options.prune_delivery_symmetry_43 &&
                            next_info.is_delivery &&
                            next_info.delivery_request_class_index >= 0) {
                            const int current_max =
                                current_label.delivery_max_request_idx_by_class
                                    [static_cast<std::size_t>(next_info.delivery_request_class_index)];
                            if (next_info.request_index <= current_max) {
                                continue;
                            }
                        }
                    }

                    const EdgeRecord& graph_edge =
                        graph.edges()[static_cast<std::size_t>(pricing_edge.graph_edge_id)];
                    double alpha_head = 0.0;
                    if (graph.is_physical_service_node(graph_edge.v)) {
                        alpha_head = lookup_dual(dual_solution.visit_duals, graph_edge.v);
                    }
                    const double edge_objective_cost =
                        phase_edge_objective_cost(phase, graph_edge.data.cost);
                    const double delta =
                        dual_solution.edge_duals[static_cast<std::size_t>(pricing_edge.graph_edge_id)];

                    Label next_label = current_label;
                    next_label.node_state_index = pricing_edge.to_node_state_index;
                    next_label.depth = current_label.depth + 1;
                    next_label.base_cost =
                        current_label.base_cost + edge_objective_cost - delta - 2.0 * alpha_head;
                    next_label.total_time = current_label.total_time + pricing_edge.time;
                    next_label.parent_label_index = label_index;
                    next_label.incoming_pricing_edge_index = pricing_edge_index;

                    if (next_info.is_physical_service_node) {
                        mask_add(next_label.visited_physical_mask, next_info.physical_bit_index);
                    }
                    if (options.prune_pickup_symmetry_43 &&
                        next_info.is_pickup &&
                        next_info.pickup_request_class_index >= 0) {
                        next_label.pickup_max_request_idx_by_class
                            [static_cast<std::size_t>(next_info.pickup_request_class_index)] =
                                next_info.request_index;
                    }
                    if (options.prune_delivery_symmetry_43 &&
                        next_info.is_delivery &&
                        next_info.delivery_request_class_index >= 0) {
                        next_label.delivery_max_request_idx_by_class
                            [static_cast<std::size_t>(next_info.delivery_request_class_index)] =
                                next_info.request_index;
                    }

                    ++result.generated_label_count;

                    std::vector<int>& bucket =
                        bucket_label_indices[static_cast<std::size_t>(next_label.depth)]
                                          [static_cast<std::size_t>(next_label.node_state_index)];

                    bool dominated = false;
                    std::vector<std::size_t> dominated_bucket_positions;
                    for (std::size_t bucket_pos = 0; bucket_pos < bucket.size(); ++bucket_pos) {
                        const int existing_label_index = bucket[bucket_pos];
                        const Label& existing_label =
                            labels[static_cast<std::size_t>(existing_label_index)];
                        if (dominates_label(existing_label, next_label)) {
                            dominated = true;
                            break;
                        }
                        if (dominates_label(next_label, existing_label)) {
                            dominated_bucket_positions.push_back(bucket_pos);
                        }
                    }

                    if (dominated) {
                        ++result.dominated_label_count;
                        continue;
                    }

                    if (!dominated_bucket_positions.empty()) {
                        // Bucket order is irrelevant; remove dominated entries by swap-pop
                        // to avoid rebuilding large buckets when only a few labels are dominated.
                        for (auto pos_it = dominated_bucket_positions.rbegin();
                             pos_it != dominated_bucket_positions.rend();
                             ++pos_it) {
                            const std::size_t bucket_pos = *pos_it;
                            bucket[bucket_pos] = bucket.back();
                            bucket.pop_back();
                            ++result.dominated_label_count;
                        }
                    }

                    const int new_label_index = static_cast<int>(labels.size());
                    labels.push_back(std::move(next_label));
                    bucket.push_back(new_label_index);
                }
            }

            if (stop_current_start) {
                break;
            }

            current_layer_indices.clear();
            if (depth < context.p) {
                const auto& next_layer_buckets =
                    bucket_label_indices[static_cast<std::size_t>(depth + 1)];
                for (const std::vector<int>& bucket : next_layer_buckets) {
                    current_layer_indices.insert(
                        current_layer_indices.end(),
                        bucket.begin(),
                        bucket.end()
                    );
                }
            }
        }

        result.surviving_label_count += count_surviving_labels(bucket_label_indices);

        if (!best_columns_for_start.empty()) {
            std::sort(
                best_columns_for_start.begin(),
                best_columns_for_start.end(),
                [](const CGColumn& lhs, const CGColumn& rhs) {
                    return lhs.reduced_cost < rhs.reduced_cost;
                }
            );
            if (best_columns_for_start.size() > options.max_columns_per_start) {
                deferred_candidate_columns.insert(
                    deferred_candidate_columns.end(),
                    best_columns_for_start.begin() +
                        static_cast<std::ptrdiff_t>(options.max_columns_per_start),
                    best_columns_for_start.end()
                );
                best_columns_for_start.resize(options.max_columns_per_start);
            }
            candidate_columns.insert(
                candidate_columns.end(),
                best_columns_for_start.begin(),
                best_columns_for_start.end()
            );
        }
    }

    std::sort(
        candidate_columns.begin(),
        candidate_columns.end(),
        [](const CGColumn& lhs, const CGColumn& rhs) {
            return lhs.reduced_cost < rhs.reduced_cost;
        }
    );
    std::sort(
        deferred_candidate_columns.begin(),
        deferred_candidate_columns.end(),
        [](const CGColumn& lhs, const CGColumn& rhs) {
            return lhs.reduced_cost < rhs.reduced_cost;
        }
    );

    std::map<std::string, int> seen_column_keys;
    for (const CGColumn& column : candidate_columns) {
        const std::string column_key = build_column_key(column.edge_ids, column.tau);
        if (seen_column_keys.find(column_key) != seen_column_keys.end()) {
            continue;
        }
        seen_column_keys.emplace(column_key, 1);
        if (result.columns.size() < options.max_total_columns) {
            result.columns.push_back(column);
        } else {
            result.deferred_columns.push_back(column);
        }
    }
    for (const CGColumn& column : deferred_candidate_columns) {
        const std::string column_key = build_column_key(column.edge_ids, column.tau);
        if (seen_column_keys.find(column_key) != seen_column_keys.end()) {
            continue;
        }
        seen_column_keys.emplace(column_key, 1);
        result.deferred_columns.push_back(column);
    }

    if (!std::isfinite(result.best_reduced_cost)) {
        result.best_reduced_cost = 0.0;
    }
    if (result.complete_label_count == 0U) {
        result.status = ForwardPricingStatus::NoCompleteLabel;
    } else if (result.columns.empty()) {
        result.status = ForwardPricingStatus::NoNegativeColumn;
    } else {
        result.status = ForwardPricingStatus::ColumnsFound;
    }

    return result;
}

ForwardPricingResult run_phase_two_shallow_search(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    const ForwardPricingOptions& options
) {
    if (context.p != 3) {
        throw std::runtime_error("3-step shallow search currently requires p = 3.");
    }

    ForwardPricingResult result;
    result.best_reduced_cost = std::numeric_limits<double>::infinity();
    const std::size_t search_max_columns_per_start =
        options.heuristic_pricing
            ? compute_heuristic_search_limit(
                  options.max_columns_per_start,
                  options.heuristic_search_column_ratio
              )
            : options.max_columns_per_start;
    const std::size_t search_max_total_columns =
        options.heuristic_pricing
            ? compute_heuristic_search_limit(
                  options.max_total_columns,
                  options.heuristic_search_column_ratio
              )
            : options.max_total_columns;

    const std::vector<int> ordered_starts =
        options.explicit_start_node_state_indices.empty()
            ? build_heuristic_start_order(
                  graph,
                  context,
                  dual_solution,
                  phase,
                  options.heuristic_start_score_mode
              )
            : options.explicit_start_node_state_indices;

    const std::size_t mask_word_count = bit_word_count(context.physical_node_count);
    std::vector<CGColumn> candidate_columns;
    std::vector<CGColumn> deferred_candidate_columns;
    candidate_columns.reserve(search_max_total_columns);
    std::size_t total_negative_column_count = 0U;

    for (int start_node_state_index : ordered_starts) {
        if (total_negative_column_count >= search_max_total_columns) {
            break;
        }

        ++result.start_label_count;
        const ForwardPricingContext::NodeStateInfo& start_info =
            context.node_states[static_cast<std::size_t>(start_node_state_index)];
        const std::size_t generated_before_start = result.generated_label_count;

        ShallowSearchState root_state;
        root_state.node_state_index = start_node_state_index;
        root_state.visited_physical_mask.assign(mask_word_count, 0U);
        root_state.pickup_max_request_idx_by_class.assign(
            context.pickup_request_class_count,
            -1
        );
        root_state.delivery_max_request_idx_by_class.assign(
            context.delivery_request_class_count,
            -1
        );
        if (start_info.is_physical_service_node) {
            mask_add(root_state.visited_physical_mask, start_info.physical_bit_index);
        }
        if (options.prune_pickup_symmetry_43 &&
            start_info.is_pickup &&
            start_info.pickup_request_class_index >= 0) {
            root_state.pickup_max_request_idx_by_class
                [static_cast<std::size_t>(start_info.pickup_request_class_index)] =
                    start_info.request_index;
        }
        if (options.prune_delivery_symmetry_43 &&
            start_info.is_delivery &&
            start_info.delivery_request_class_index >= 0) {
            root_state.delivery_max_request_idx_by_class
                [static_cast<std::size_t>(start_info.delivery_request_class_index)] =
                    start_info.request_index;
        }

        std::unordered_map<std::string, std::size_t> seen_columns_for_start;
        std::vector<CGColumn> generated_columns_for_start;
        std::vector<CGColumn> accepted_columns_for_start;
        bool stop_current_start = false;

        const std::vector<int> first_edges = collect_feasible_shallow_edges(
            graph,
            context,
            dual_solution,
            phase,
            options,
            root_state,
            options.shallow_k1
        );
        for (int first_pricing_edge_index : first_edges) {
            const ForwardPricingContext::EdgeInfo& first_edge =
                context.edges[static_cast<std::size_t>(first_pricing_edge_index)];
            ShallowSearchState first_state = extend_shallow_state(
                context,
                first_edge,
                root_state,
                options
            );
            ++result.generated_label_count;

            if (is_complete_label(start_info.node_id, first_state.depth, context.p) &&
                evaluate_shallow_complete_state(
                    graph,
                    context,
                    dual_solution,
                    phase,
                    options,
                    first_state,
                    seen_columns_for_start,
                    generated_columns_for_start,
                    accepted_columns_for_start,
                    search_max_columns_per_start,
                    result
                )) {
                stop_current_start = true;
            }
            if (stop_current_start ||
                total_negative_column_count + accepted_columns_for_start.size() >=
                    search_max_total_columns) {
                break;
            }

            const std::vector<int> second_edges = collect_feasible_shallow_edges(
                graph,
                context,
                dual_solution,
                phase,
                options,
                first_state,
                options.shallow_k2
            );
            for (int second_pricing_edge_index : second_edges) {
                const ForwardPricingContext::EdgeInfo& second_edge =
                    context.edges[static_cast<std::size_t>(second_pricing_edge_index)];
                ShallowSearchState second_state = extend_shallow_state(
                    context,
                    second_edge,
                    first_state,
                    options
                );
                ++result.generated_label_count;

                if (is_complete_label(start_info.node_id, second_state.depth, context.p) &&
                    evaluate_shallow_complete_state(
                        graph,
                        context,
                        dual_solution,
                        phase,
                        options,
                        second_state,
                        seen_columns_for_start,
                        generated_columns_for_start,
                        accepted_columns_for_start,
                        search_max_columns_per_start,
                        result
                    )) {
                    stop_current_start = true;
                }
                if (stop_current_start ||
                    total_negative_column_count + accepted_columns_for_start.size() >=
                        search_max_total_columns) {
                    break;
                }

                const std::vector<int> third_edges = collect_feasible_shallow_edges(
                    graph,
                    context,
                    dual_solution,
                    phase,
                    options,
                    second_state,
                    0U
                );
                for (int third_pricing_edge_index : third_edges) {
                    const ForwardPricingContext::EdgeInfo& third_edge =
                        context.edges[static_cast<std::size_t>(third_pricing_edge_index)];
                    ShallowSearchState third_state = extend_shallow_state(
                        context,
                        third_edge,
                        second_state,
                        options
                    );
                    ++result.generated_label_count;

                    if (is_complete_label(start_info.node_id, third_state.depth, context.p) &&
                        evaluate_shallow_complete_state(
                            graph,
                            context,
                            dual_solution,
                            phase,
                            options,
                            third_state,
                            seen_columns_for_start,
                            generated_columns_for_start,
                            accepted_columns_for_start,
                            search_max_columns_per_start,
                            result
                        )) {
                        stop_current_start = true;
                    }
                    if (stop_current_start ||
                        total_negative_column_count + accepted_columns_for_start.size() >=
                            search_max_total_columns) {
                        break;
                    }
                }

                if (stop_current_start ||
                    total_negative_column_count + accepted_columns_for_start.size() >=
                        search_max_total_columns) {
                    break;
                }
            }

            if (stop_current_start ||
                total_negative_column_count + accepted_columns_for_start.size() >=
                    search_max_total_columns) {
                break;
            }
        }

        result.surviving_label_count +=
            1U + (result.generated_label_count - generated_before_start);

        if (!accepted_columns_for_start.empty()) {
            total_negative_column_count += accepted_columns_for_start.size();
            std::sort(
                accepted_columns_for_start.begin(),
                accepted_columns_for_start.end(),
                [](const CGColumn& lhs, const CGColumn& rhs) {
                    return lhs.reduced_cost < rhs.reduced_cost;
                }
            );
            if (accepted_columns_for_start.size() > options.max_columns_per_start) {
                deferred_candidate_columns.insert(
                    deferred_candidate_columns.end(),
                    accepted_columns_for_start.begin() +
                        static_cast<std::ptrdiff_t>(options.max_columns_per_start),
                    accepted_columns_for_start.end()
                );
                accepted_columns_for_start.resize(options.max_columns_per_start);
            }
            candidate_columns.insert(
                candidate_columns.end(),
                accepted_columns_for_start.begin(),
                accepted_columns_for_start.end()
            );
        }
    }

    std::sort(
        candidate_columns.begin(),
        candidate_columns.end(),
        [](const CGColumn& lhs, const CGColumn& rhs) {
            return lhs.reduced_cost < rhs.reduced_cost;
        }
    );
    std::sort(
        deferred_candidate_columns.begin(),
        deferred_candidate_columns.end(),
        [](const CGColumn& lhs, const CGColumn& rhs) {
            return lhs.reduced_cost < rhs.reduced_cost;
        }
    );

    std::map<std::string, int> seen_column_keys;
    for (const CGColumn& column : candidate_columns) {
        const std::string column_key = build_column_key(column.edge_ids, column.tau);
        if (seen_column_keys.find(column_key) != seen_column_keys.end()) {
            continue;
        }
        seen_column_keys.emplace(column_key, 1);
        if (result.columns.size() < options.max_total_columns) {
            result.columns.push_back(column);
        } else {
            result.deferred_columns.push_back(column);
        }
    }
    for (const CGColumn& column : deferred_candidate_columns) {
        const std::string column_key = build_column_key(column.edge_ids, column.tau);
        if (seen_column_keys.find(column_key) != seen_column_keys.end()) {
            continue;
        }
        seen_column_keys.emplace(column_key, 1);
        result.deferred_columns.push_back(column);
    }

    if (!std::isfinite(result.best_reduced_cost)) {
        result.best_reduced_cost = 0.0;
    }
    if (result.complete_label_count == 0U) {
        result.status = ForwardPricingStatus::NoCompleteLabel;
    } else if (result.columns.empty()) {
        result.status = ForwardPricingStatus::NoNegativeColumn;
    } else {
        result.status = ForwardPricingStatus::ColumnsFound;
    }

    return result;
}

void write_generated_columns(
    std::ostream& out,
    const std::vector<CGColumn>& columns
) {
    out << "[pricing] Generated column count: " << columns.size() << '\n';
    for (const CGColumn& column : columns) {
        out << "  rc=" << format_double(column.reduced_cost)
            << " q=" << column.q
            << " tau=" << format_double(column.tau)
            << " time=" << format_double(column.total_time)
            << " cost=" << format_double(column.total_cost)
            << " edges=[";
        for (std::size_t idx = 0; idx < column.edge_ids.size(); ++idx) {
            if (idx > 0U) {
                out << ", ";
            }
            out << column.edge_ids[idx];
        }
        out << "]";
        out << " visit={";
        for (std::size_t idx = 0; idx < column.visit_coefficients.size(); ++idx) {
            if (idx > 0U) {
                out << ", ";
            }
            out << column.visit_coefficients[idx].first << ":" << column.visit_coefficients[idx].second;
        }
        out << "} state={";
        for (std::size_t idx = 0; idx < column.state_coefficients.size(); ++idx) {
            if (idx > 0U) {
                out << ", ";
            }
            out << "(" << column.state_coefficients[idx].first.node_id << ","
                << state_to_str(column.state_coefficients[idx].first.state) << "):"
                << column.state_coefficients[idx].second;
        }
        out << "} time={";
        for (std::size_t idx = 0; idx < column.time_coefficients.size(); ++idx) {
            if (idx > 0U) {
                out << ", ";
            }
            out << "(" << column.time_coefficients[idx].first.node_id << ","
                << state_to_str(column.time_coefficients[idx].first.state) << "):"
                << format_double(column.time_coefficients[idx].second);
        }
        out << "}\n";
    }
}

std::vector<int> build_heuristic_start_order(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    HeuristicStartScoreMode mode
) {
    return build_heuristic_start_order_impl(
        graph,
        context,
        dual_solution,
        phase,
        mode
    );
}

double evaluate_column_reduced_cost(
    const CGDualSolution& dual_solution,
    CGPhase phase,
    const CGColumn& column
) {
    return evaluate_sparse_column_reduced_cost(
        dual_solution,
        phase,
        column.total_cost,
        column.visit_coefficients,
        column.state_coefficients,
        column.time_coefficients,
        column.edge_incidence
    );
}

const char* to_string(ForwardPricingStatus status) {
    switch (status) {
        case ForwardPricingStatus::ColumnsFound:
            return "columns-found";
        case ForwardPricingStatus::NoNegativeColumn:
            return "no-negative-column";
        case ForwardPricingStatus::NoCompleteLabel:
            return "no-complete-label";
    }
    return "unknown";
}

const char* to_string(HeuristicStartScoreMode mode) {
    switch (mode) {
        case HeuristicStartScoreMode::OneStepMin:
            return "one-step-min";
    }
    return "unknown";
}

}  // namespace spdp
