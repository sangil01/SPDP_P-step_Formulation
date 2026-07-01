#include "PstepPricing.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <map>
#include <optional>
#include <ostream>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef SPDP_HAVE_ONEMKL
#include <mkl.h>
#include <mkl_spblas.h>
#endif

namespace spdp {
namespace {

constexpr double kTolerance = 1e-9;

bool double_equal(double lhs, double rhs) {
    return std::fabs(lhs - rhs) <= kTolerance;
}

bool is_forbidden_pricing_edge(
    const ForwardPricingOptions& options,
    int graph_edge_id
) {
    return graph_edge_id >= 0 &&
           static_cast<std::size_t>(graph_edge_id) < options.forbidden_edge_mask.size() &&
           options.forbidden_edge_mask[static_cast<std::size_t>(graph_edge_id)] != 0U;
}

bool violates_required_pricing_edge_rule(
    const ForwardPricingOptions& options,
    const MultiDiGraph& graph,
    int graph_edge_id
) {
    if (graph_edge_id < 0 ||
        static_cast<std::size_t>(graph_edge_id) >= graph.edges().size()) {
        return false;
    }

    const EdgeRecord& edge = graph.edges()[static_cast<std::size_t>(graph_edge_id)];
    if (edge.u >= 0 &&
        static_cast<std::size_t>(edge.u) < options.required_outgoing_edge_by_node_id.size()) {
        const int required_edge_id =
            options.required_outgoing_edge_by_node_id[static_cast<std::size_t>(edge.u)];
        if (required_edge_id >= 0 && required_edge_id != graph_edge_id) {
            return true;
        }
    }
    if (edge.v >= 0 &&
        static_cast<std::size_t>(edge.v) < options.required_incoming_edge_by_node_id.size()) {
        const int required_edge_id =
            options.required_incoming_edge_by_node_id[static_cast<std::size_t>(edge.v)];
        if (required_edge_id >= 0 && required_edge_id != graph_edge_id) {
            return true;
        }
    }
    return false;
}

bool is_disallowed_pricing_edge(
    const ForwardPricingOptions& options,
    const MultiDiGraph& graph,
    int graph_edge_id
) {
    return is_forbidden_pricing_edge(options, graph_edge_id) ||
           violates_required_pricing_edge_rule(options, graph, graph_edge_id);
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

struct TopKForwardEdgeSelection {
    std::vector<int> pricing_edge_indices;
    std::size_t feasible_edge_count_before_top_k = 0;
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
        if (is_disallowed_pricing_edge(options, graph, pricing_edge.graph_edge_id)) {
            continue;
        }
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

TopKForwardEdgeSelection collect_top_k_forward_labeling_edges(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    const ForwardPricingOptions& options,
    const Label& current_label,
    std::size_t top_k_next
) {
    TopKForwardEdgeSelection selection;
    std::vector<std::pair<double, int>> scored_edges;
    scored_edges.reserve(
        context.outgoing_edges_by_node_state[static_cast<std::size_t>(current_label.node_state_index)]
            .size()
    );

    for (int pricing_edge_index :
         context.outgoing_edges_by_node_state[static_cast<std::size_t>(current_label.node_state_index)]) {
        const ForwardPricingContext::EdgeInfo& pricing_edge =
            context.edges[static_cast<std::size_t>(pricing_edge_index)];
        if (is_disallowed_pricing_edge(options, graph, pricing_edge.graph_edge_id)) {
            continue;
        }
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

    selection.feasible_edge_count_before_top_k = scored_edges.size();

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

    selection.pricing_edge_indices.reserve(scored_edges.size());
    for (const auto& entry : scored_edges) {
        selection.pricing_edge_indices.push_back(entry.second);
    }
    return selection;
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
    const std::vector<int> ordered_start_node_state_indices =
        build_start_order(graph, context, dual_solution, phase, options);

    for (int start_node_state_index : ordered_start_node_state_indices) {
        if (options.heuristic_pricing &&
            result.total_negative_column_count >= search_max_total_columns) {
            result.hit_global_search_cap = true;
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
        bool hit_start_search_cap = false;

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
                        ++result.total_negative_column_count;
                        if (options.heuristic_pricing &&
                            best_columns_for_start.size() >= search_max_columns_per_start) {
                            hit_start_search_cap = true;
                            stop_current_start = true;
                            break;
                        }
                        if (options.heuristic_pricing &&
                            result.total_negative_column_count >= search_max_total_columns) {
                            result.hit_global_search_cap = true;
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
                    const TopKForwardEdgeSelection top_k_selection =
                        collect_top_k_forward_labeling_edges(
                        graph,
                        context,
                        dual_solution,
                        phase,
                        options,
                        current_label,
                        options.labeling_top_k_next
                    );
                    ++result.top_k_next_applied_label_count;
                    result.top_k_next_feasible_edges_before +=
                        top_k_selection.feasible_edge_count_before_top_k;
                    result.top_k_next_feasible_edges_after +=
                        top_k_selection.pricing_edge_indices.size();
                    top_k_pricing_edge_indices = top_k_selection.pricing_edge_indices;
                    pricing_edge_indices = &top_k_pricing_edge_indices;
                }

                for (int pricing_edge_index : *pricing_edge_indices) {
                    const ForwardPricingContext::EdgeInfo& pricing_edge =
                        context.edges[static_cast<std::size_t>(pricing_edge_index)];
                    if (is_disallowed_pricing_edge(options, graph, pricing_edge.graph_edge_id)) {
                        continue;
                    }
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

        if (hit_start_search_cap) {
            ++result.per_start_search_cap_hit_count;
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

    for (int start_node_state_index : ordered_starts) {
        if (result.total_negative_column_count >= search_max_total_columns) {
            result.hit_global_search_cap = true;
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
                result.total_negative_column_count + accepted_columns_for_start.size() >=
                    search_max_total_columns) {
                if (result.total_negative_column_count + accepted_columns_for_start.size() >=
                    search_max_total_columns) {
                    result.hit_global_search_cap = true;
                }
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
                    result.total_negative_column_count + accepted_columns_for_start.size() >=
                        search_max_total_columns) {
                    if (result.total_negative_column_count + accepted_columns_for_start.size() >=
                        search_max_total_columns) {
                        result.hit_global_search_cap = true;
                    }
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
                        result.total_negative_column_count + accepted_columns_for_start.size() >=
                            search_max_total_columns) {
                        if (result.total_negative_column_count + accepted_columns_for_start.size() >=
                            search_max_total_columns) {
                            result.hit_global_search_cap = true;
                        }
                        break;
                    }
                }

                if (stop_current_start ||
                    result.total_negative_column_count + accepted_columns_for_start.size() >=
                        search_max_total_columns) {
                    if (result.total_negative_column_count + accepted_columns_for_start.size() >=
                        search_max_total_columns) {
                        result.hit_global_search_cap = true;
                    }
                    break;
                }
            }

            if (stop_current_start ||
                result.total_negative_column_count + accepted_columns_for_start.size() >=
                    search_max_total_columns) {
                if (result.total_negative_column_count + accepted_columns_for_start.size() >=
                    search_max_total_columns) {
                    result.hit_global_search_cap = true;
                }
                break;
            }
        }

        result.surviving_label_count +=
            1U + (result.generated_label_count - generated_before_start);

        if (!accepted_columns_for_start.empty()) {
            result.total_negative_column_count += accepted_columns_for_start.size();
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
        if (stop_current_start) {
            ++result.per_start_search_cap_hit_count;
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

namespace {

struct FullEnumerationEntryMatrixBuildData {
    std::vector<std::pair<int, int>> visit_terms;
    std::vector<std::pair<int, int>> state_terms;
    std::vector<int> edge_terms;
};

int ensure_visit_index(
    FullEnumerationStaticPool& pool,
    NodeId node_id
) {
    const auto found = pool.visit_index_by_node.find(node_id);
    if (found != pool.visit_index_by_node.end()) {
        return found->second;
    }
    const int index = static_cast<int>(pool.visit_node_ids.size());
    pool.visit_index_by_node.emplace(node_id, index);
    pool.visit_node_ids.push_back(node_id);
    return index;
}

int ensure_node_state_index(
    FullEnumerationStaticPool& pool,
    const NodeStateKey& key
) {
    const NodeStateKey canonical_key{key.node_id, canonicalize_state(key.state)};
    const auto found = pool.node_state_index_by_key.find(canonical_key);
    if (found != pool.node_state_index_by_key.end()) {
        return found->second;
    }
    const int index = static_cast<int>(pool.node_state_keys.size());
    pool.node_state_index_by_key.emplace(canonical_key, index);
    pool.node_state_keys.push_back(canonical_key);
    return index;
}

std::size_t full_enumeration_entry_count(
    const FullEnumerationStaticPool& pool
) {
    return pool.entry_total_cost.size();
}

std::size_t full_enumeration_variant_count(
    const FullEnumerationStaticPool& pool
) {
    return pool.variant_tau.size();
}

std::size_t full_enumeration_variant_entry_index(
    std::size_t variant_index
) {
    return variant_index / 2U;
}

std::size_t full_enumeration_variant_tau_index(
    std::size_t variant_index
) {
    return variant_index % 2U;
}

std::size_t full_enumeration_variant_index(
    std::size_t entry_index,
    std::size_t tau_index
) {
    return 2U * entry_index + tau_index;
}

bool full_enumeration_variant_is_valid(
    const FullEnumerationStaticPool& pool,
    std::size_t variant_index
) {
    return pool.variant_valid[variant_index] != 0U;
}

bool full_enumeration_variant_is_active(
    const FullEnumerationPoolNodeState& node_state,
    std::size_t variant_index
) {
    return node_state.variant_active[variant_index] != 0U;
}

bool full_enumeration_entry_is_forbidden(
    const FullEnumerationPoolNodeState& node_state,
    std::size_t entry_index
) {
    return !node_state.forbidden_entry.empty() &&
           node_state.forbidden_entry[entry_index] != 0U;
}

bool full_enumeration_variant_is_available(
    const FullEnumerationStaticPool& static_pool,
    const FullEnumerationPoolNodeState& node_state,
    std::size_t variant_index
) {
    if (!full_enumeration_variant_is_valid(static_pool, variant_index) ||
        full_enumeration_variant_is_active(node_state, variant_index)) {
        return false;
    }
    return !full_enumeration_entry_is_forbidden(
        node_state,
        full_enumeration_variant_entry_index(variant_index)
    );
}

double full_enumeration_entry_latest_start(
    const FullEnumerationStaticPool& pool,
    std::size_t entry_index
) {
    return std::max(0.0, pool.time_limit - pool.entry_total_time[entry_index]);
}

std::size_t full_enumeration_selected_variant_index(
    const FullEnumerationStaticPool& pool,
    const FullEnumerationPoolNodeState& node_state,
    const std::vector<double>& time_duals,
    std::size_t entry_index
) {
    if (full_enumeration_entry_is_forbidden(node_state, entry_index)) {
        return std::numeric_limits<std::size_t>::max();
    }
    const std::size_t slot0 = full_enumeration_variant_index(entry_index, 0U);
    const std::size_t slot1 = full_enumeration_variant_index(entry_index, 1U);
    const bool valid0 = full_enumeration_variant_is_available(pool, node_state, slot0);
    const bool valid1 = full_enumeration_variant_is_available(pool, node_state, slot1);

    if (!valid0 && !valid1) {
        return std::numeric_limits<std::size_t>::max();
    }
    if (valid0 && !valid1) {
        return slot0;
    }
    if (!valid0 && valid1) {
        return slot1;
    }

    const int start_row_index = pool.entry_start_time_row_index[entry_index];
    const int last_row_index = pool.entry_last_time_row_index[entry_index];
    const double gamma_start =
        start_row_index >= 0
            ? time_duals[static_cast<std::size_t>(start_row_index)]
            : 0.0;
    const double gamma_end =
        last_row_index >= 0
            ? time_duals[static_cast<std::size_t>(last_row_index)]
            : 0.0;

    return (gamma_end - gamma_start < 0.0) ? slot1 : slot0;
}

void append_full_enumeration_entry(
    FullEnumerationStaticPool& pool,
    std::vector<FullEnumerationEntryMatrixBuildData>& entry_matrix_build_data,
    const MultiDiGraph& graph,
    const CompactPStep& pstep,
    const CompactPStepCoefficients& coefficients
) {
    const State canonical_start_state = canonicalize_state(pstep.start_state);
    const State canonical_last_state = canonicalize_state(pstep.last_state);
    int start_time_row_index = -1;
    int last_time_row_index = -1;
    if (graph.is_physical_service_node(pstep.start_node_id)) {
        start_time_row_index = ensure_node_state_index(
            pool,
            NodeStateKey{pstep.start_node_id, canonical_start_state}
        );
    }
    if (graph.is_physical_service_node(pstep.last_node_id)) {
        last_time_row_index = ensure_node_state_index(
            pool,
            NodeStateKey{pstep.last_node_id, canonical_last_state}
        );
    }

    pool.entry_q.push_back(pstep.q);
    pool.entry_start_node_id.push_back(pstep.start_node_id);
    pool.entry_last_node_id.push_back(pstep.last_node_id);
    pool.entry_start_state.push_back(canonical_start_state);
    pool.entry_last_state.push_back(canonical_last_state);
    pool.entry_total_time.push_back(pstep.total_time);
    pool.entry_total_cost.push_back(pstep.total_cost);
    pool.entry_start_time_row_index.push_back(start_time_row_index);
    pool.entry_last_time_row_index.push_back(last_time_row_index);

    pool.entry_edge_ids.insert(
        pool.entry_edge_ids.end(),
        pstep.edge_ids.begin(),
        pstep.edge_ids.end()
    );
    pool.entry_edge_ids_row_ptr.push_back(pool.entry_edge_ids.size());

    pool.entry_node_sequence.insert(
        pool.entry_node_sequence.end(),
        pstep.node_sequence.begin(),
        pstep.node_sequence.end()
    );
    pool.entry_node_sequence_row_ptr.push_back(pool.entry_node_sequence.size());

    pool.entry_state_sequence.insert(
        pool.entry_state_sequence.end(),
        pstep.state_sequence.begin(),
        pstep.state_sequence.end()
    );
    pool.entry_state_sequence_row_ptr.push_back(pool.entry_state_sequence.size());

    pool.variant_tau.push_back(0.0);
    pool.variant_tau.push_back(0.0);
    pool.variant_valid.push_back(0U);
    pool.variant_valid.push_back(0U);

    FullEnumerationEntryMatrixBuildData entry_terms;
    entry_terms.edge_terms =
        coefficients.edge_incidence_by_pstep.at(static_cast<std::size_t>(pstep.id));

    const auto& visit_coefficients =
        coefficients.visit_coefficients_by_pstep.at(static_cast<std::size_t>(pstep.id));
    for (const auto& term : visit_coefficients) {
        entry_terms.visit_terms.push_back({ensure_visit_index(pool, term.first), term.second});
    }

    const auto& state_coefficients =
        coefficients.state_coefficients_by_pstep.at(static_cast<std::size_t>(pstep.id));
    for (const auto& term : state_coefficients) {
        const NodeStateKey key{term.first.node_id, canonicalize_state(term.first.state)};
        entry_terms.state_terms.push_back({ensure_node_state_index(pool, key), term.second});
    }

    entry_matrix_build_data.push_back(std::move(entry_terms));
}

void finalize_full_enumeration_stage1_matrix(
    FullEnumerationStaticPool& pool,
    const std::vector<FullEnumerationEntryMatrixBuildData>& entry_matrix_build_data,
    std::size_t edge_column_count
) {
    pool.stage1_visit_column_count = pool.visit_node_ids.size();
    pool.stage1_state_column_count = pool.node_state_keys.size();
    pool.stage1_edge_column_count = edge_column_count;
    pool.stage1_matrix.column_count =
        pool.stage1_visit_column_count +
        pool.stage1_state_column_count +
        pool.stage1_edge_column_count;
    pool.stage1_matrix.row_ptr.clear();
    pool.stage1_matrix.column_index.clear();
    pool.stage1_matrix.value.clear();
    pool.stage1_matrix.row_ptr.reserve(entry_matrix_build_data.size() + 1U);
    pool.stage1_matrix.row_ptr.push_back(0U);

    const int state_column_offset = static_cast<int>(pool.stage1_visit_column_count);
    const int edge_column_offset = static_cast<int>(
        pool.stage1_visit_column_count + pool.stage1_state_column_count
    );
    for (const FullEnumerationEntryMatrixBuildData& entry_terms : entry_matrix_build_data) {
        for (const auto& term : entry_terms.visit_terms) {
            pool.stage1_matrix.column_index.push_back(term.first);
            pool.stage1_matrix.value.push_back(static_cast<double>(term.second));
        }
        for (const auto& term : entry_terms.state_terms) {
            pool.stage1_matrix.column_index.push_back(state_column_offset + term.first);
            pool.stage1_matrix.value.push_back(static_cast<double>(term.second));
        }
        for (int edge_id : entry_terms.edge_terms) {
            pool.stage1_matrix.column_index.push_back(edge_column_offset + edge_id);
            pool.stage1_matrix.value.push_back(1.0);
        }
        pool.stage1_matrix.row_ptr.push_back(pool.stage1_matrix.column_index.size());
    }
}

void finalize_full_enumeration_stage2_matrix(FullEnumerationStaticPool& pool) {
    const std::size_t entry_count = full_enumeration_entry_count(pool);
    pool.stage2_matrix.column_count = entry_count + pool.node_state_keys.size();
    pool.stage2_matrix.row_ptr.clear();
    pool.stage2_matrix.column_index.clear();
    pool.stage2_matrix.value.clear();
    pool.stage2_matrix.row_ptr.reserve(full_enumeration_variant_count(pool) + 1U);
    pool.stage2_matrix.row_ptr.push_back(0U);

    for (std::size_t variant_index = 0; variant_index < full_enumeration_variant_count(pool);
         ++variant_index) {
        if (!full_enumeration_variant_is_valid(pool, variant_index)) {
            pool.stage2_matrix.row_ptr.push_back(pool.stage2_matrix.column_index.size());
            continue;
        }

        const std::size_t entry_index = full_enumeration_variant_entry_index(variant_index);
        pool.stage2_matrix.column_index.push_back(static_cast<int>(entry_index));
        pool.stage2_matrix.value.push_back(1.0);

        const int start_time_row_index = pool.entry_start_time_row_index[entry_index];
        if (start_time_row_index >= 0) {
            pool.stage2_matrix.column_index.push_back(
                static_cast<int>(entry_count + static_cast<std::size_t>(start_time_row_index))
            );
            pool.stage2_matrix.value.push_back(-pool.variant_tau[variant_index]);
        }

        const int last_time_row_index = pool.entry_last_time_row_index[entry_index];
        if (last_time_row_index >= 0) {
            pool.stage2_matrix.column_index.push_back(
                static_cast<int>(entry_count + static_cast<std::size_t>(last_time_row_index))
            );
            pool.stage2_matrix.value.push_back(
                pool.variant_tau[variant_index] + pool.entry_total_time[entry_index]
            );
        }

        pool.stage2_matrix.row_ptr.push_back(pool.stage2_matrix.column_index.size());
    }
}

void finalize_full_enumeration_direct_visit_lookup(FullEnumerationStaticPool& pool) {
    NodeId max_node_id = -1;
    for (NodeId node_id : pool.visit_node_ids) {
        max_node_id = std::max(max_node_id, node_id);
    }

    if (max_node_id < 0) {
        pool.dense_visit_row_index_by_node_id.clear();
        return;
    }

    pool.dense_visit_row_index_by_node_id.assign(
        static_cast<std::size_t>(max_node_id) + 1U,
        -1
    );
    for (std::size_t visit_index = 0; visit_index < pool.visit_node_ids.size(); ++visit_index) {
        const NodeId node_id = pool.visit_node_ids[visit_index];
        if (node_id >= 0) {
            pool.dense_visit_row_index_by_node_id[static_cast<std::size_t>(node_id)] =
                static_cast<int>(visit_index);
        }
    }
}

struct FullPoolSelection {
    double reduced_cost = 0.0;
    std::size_t variant_index = 0;
    std::size_t entry_index = 0;
    std::size_t tau_index = 0;
};

struct FullPoolWorstReducedCostFirst {
    bool operator()(const FullPoolSelection& lhs, const FullPoolSelection& rhs) const {
        return lhs.reduced_cost < rhs.reduced_cost;
    }
};

using FullPoolTopKQueue =
    std::priority_queue<
        FullPoolSelection,
        std::vector<FullPoolSelection>,
        FullPoolWorstReducedCostFirst>;

CGColumn build_column_from_full_pool_selection(
    const FullEnumerationStaticPool& pool,
    const FullPoolSelection& selection
) {
    const std::size_t entry_index = selection.entry_index;
    const std::size_t variant_index = selection.variant_index;

    CGColumn column;
    column.q = pool.entry_q[entry_index];
    column.start_node_id = pool.entry_start_node_id[entry_index];
    column.last_node_id = pool.entry_last_node_id[entry_index];
    column.start_state = pool.entry_start_state[entry_index];
    column.last_state = pool.entry_last_state[entry_index];
    column.total_time = pool.entry_total_time[entry_index];
    column.total_cost = pool.entry_total_cost[entry_index];
    column.tau = pool.variant_tau[variant_index];
    column.reduced_cost = selection.reduced_cost;
    column.pool_variant_index = static_cast<int>(variant_index);

    for (std::size_t k = pool.entry_edge_ids_row_ptr[entry_index];
         k < pool.entry_edge_ids_row_ptr[entry_index + 1U];
         ++k) {
        column.edge_ids.push_back(pool.entry_edge_ids[k]);
    }
    for (std::size_t k = pool.entry_node_sequence_row_ptr[entry_index];
         k < pool.entry_node_sequence_row_ptr[entry_index + 1U];
         ++k) {
        column.node_sequence.push_back(pool.entry_node_sequence[k]);
    }
    for (std::size_t k = pool.entry_state_sequence_row_ptr[entry_index];
         k < pool.entry_state_sequence_row_ptr[entry_index + 1U];
         ++k) {
        column.state_sequence.push_back(pool.entry_state_sequence[k]);
    }

    const std::size_t visit_limit = pool.stage1_visit_column_count;
    const std::size_t state_limit =
        pool.stage1_visit_column_count + pool.stage1_state_column_count;
    for (std::size_t k = pool.stage1_matrix.row_ptr[entry_index];
         k < pool.stage1_matrix.row_ptr[entry_index + 1U];
         ++k) {
        const std::size_t column_index =
            static_cast<std::size_t>(pool.stage1_matrix.column_index[k]);
        const int coefficient = static_cast<int>(std::llround(pool.stage1_matrix.value[k]));
        if (column_index < visit_limit) {
            column.visit_coefficients.push_back({
                pool.visit_node_ids[column_index],
                coefficient,
            });
        } else if (column_index < state_limit) {
            const std::size_t state_index = column_index - visit_limit;
            column.state_coefficients.push_back({
                pool.node_state_keys[state_index],
                coefficient,
            });
        } else {
            column.edge_incidence.push_back(static_cast<int>(column_index - state_limit));
        }
    }

    const int start_time_row_index = pool.entry_start_time_row_index[entry_index];
    if (start_time_row_index >= 0) {
        column.time_coefficients.push_back({
            pool.node_state_keys[static_cast<std::size_t>(start_time_row_index)],
            pool.variant_tau[variant_index],
        });
    }
    const int last_time_row_index = pool.entry_last_time_row_index[entry_index];
    if (last_time_row_index >= 0) {
        column.time_coefficients.push_back(
            {
                pool.node_state_keys[static_cast<std::size_t>(last_time_row_index)],
                -(pool.variant_tau[variant_index] + pool.entry_total_time[entry_index]),
            }
        );
    }
    return column;
}

struct FullEnumerationDenseDualVectors {
    std::vector<double> stage1_input;
    std::vector<double> time;
};

FullEnumerationDenseDualVectors build_full_enumeration_dense_duals(
    const FullEnumerationStaticPool& pool,
    const CGDualSolution& dual_solution
) {
    FullEnumerationDenseDualVectors dense_duals;
    dense_duals.stage1_input.assign(pool.stage1_matrix.column_count, 0.0);
    for (std::size_t idx = 0; idx < pool.visit_node_ids.size(); ++idx) {
        dense_duals.stage1_input[idx] =
            lookup_dual(dual_solution.visit_duals, pool.visit_node_ids[idx]);
    }
    dense_duals.time.assign(pool.node_state_keys.size(), 0.0);
    for (std::size_t idx = 0; idx < pool.node_state_keys.size(); ++idx) {
        dense_duals.stage1_input[pool.stage1_visit_column_count + idx] =
            lookup_dual(dual_solution.state_duals, pool.node_state_keys[idx]);
        dense_duals.time[idx] = lookup_dual(dual_solution.time_duals, pool.node_state_keys[idx]);
    }
    const std::size_t edge_offset =
        pool.stage1_visit_column_count + pool.stage1_state_column_count;
    const std::size_t edge_count = std::min(
        pool.stage1_edge_column_count,
        dual_solution.edge_duals.size()
    );
    for (std::size_t idx = 0; idx < edge_count; ++idx) {
        dense_duals.stage1_input[edge_offset + idx] = dual_solution.edge_duals[idx];
    }
    return dense_duals;
}

struct FullEnumerationBaseRCBreakdown {
    double entry_cost_term = 0.0;
    double visit_contribution = 0.0;
    double state_contribution = 0.0;
    double edge_contribution = 0.0;
    double base_reduced_cost = 0.0;
};

double compute_full_enumeration_entry_base_reduced_cost_direct(
    const FullEnumerationStaticPool& pool,
    const FullEnumerationDenseDualVectors& dense_duals,
    CGPhase phase,
    std::size_t entry_index
) {
    const double entry_cost_term =
        (phase == CGPhase::PhaseII) ? pool.entry_total_cost[entry_index] : 0.0;
    double row_sum = 0.0;

    const std::size_t node_begin = pool.entry_node_sequence_row_ptr[entry_index];
    const std::size_t node_end = pool.entry_node_sequence_row_ptr[entry_index + 1U];
    for (std::size_t pos = node_begin; pos < node_end; ++pos) {
        const NodeId node_id = pool.entry_node_sequence[pos];
        if (node_id < 0 ||
            static_cast<std::size_t>(node_id) >= pool.dense_visit_row_index_by_node_id.size()) {
            continue;
        }
        const int visit_row_index =
            pool.dense_visit_row_index_by_node_id[static_cast<std::size_t>(node_id)];
        if (visit_row_index < 0) {
            continue;
        }
        const double coefficient =
            (pos == node_begin || pos + 1U == node_end) ? 1.0 : 2.0;
        row_sum += coefficient * dense_duals.stage1_input[static_cast<std::size_t>(visit_row_index)];
    }

    const std::size_t state_offset = pool.stage1_visit_column_count;
    const int start_state_row_index = pool.entry_start_time_row_index[entry_index];
    if (start_state_row_index >= 0) {
        row_sum += dense_duals.stage1_input[
            state_offset + static_cast<std::size_t>(start_state_row_index)
        ];
    }
    const int last_state_row_index = pool.entry_last_time_row_index[entry_index];
    if (last_state_row_index >= 0) {
        row_sum -= dense_duals.stage1_input[
            state_offset + static_cast<std::size_t>(last_state_row_index)
        ];
    }

    const std::size_t edge_offset =
        pool.stage1_visit_column_count + pool.stage1_state_column_count;
    for (std::size_t k = pool.entry_edge_ids_row_ptr[entry_index];
         k < pool.entry_edge_ids_row_ptr[entry_index + 1U];
         ++k) {
        const int edge_id = pool.entry_edge_ids[k];
        if (edge_id < 0 ||
            static_cast<std::size_t>(edge_id) >= pool.stage1_edge_column_count) {
            continue;
        }
        row_sum += dense_duals.stage1_input[edge_offset + static_cast<std::size_t>(edge_id)];
    }

    return entry_cost_term - row_sum;
}

FullEnumerationBaseRCBreakdown compute_full_enumeration_entry_base_reduced_cost_direct_breakdown(
    const FullEnumerationStaticPool& pool,
    const FullEnumerationDenseDualVectors& dense_duals,
    CGPhase phase,
    std::size_t entry_index
) {
    FullEnumerationBaseRCBreakdown breakdown;
    breakdown.entry_cost_term =
        (phase == CGPhase::PhaseII) ? pool.entry_total_cost[entry_index] : 0.0;

    const std::size_t node_begin = pool.entry_node_sequence_row_ptr[entry_index];
    const std::size_t node_end = pool.entry_node_sequence_row_ptr[entry_index + 1U];
    for (std::size_t pos = node_begin; pos < node_end; ++pos) {
        const NodeId node_id = pool.entry_node_sequence[pos];
        if (node_id < 0 ||
            static_cast<std::size_t>(node_id) >= pool.dense_visit_row_index_by_node_id.size()) {
            continue;
        }
        const int visit_row_index =
            pool.dense_visit_row_index_by_node_id[static_cast<std::size_t>(node_id)];
        if (visit_row_index < 0) {
            continue;
        }
        const double coefficient =
            (pos == node_begin || pos + 1U == node_end) ? 1.0 : 2.0;
        breakdown.visit_contribution +=
            coefficient * dense_duals.stage1_input[static_cast<std::size_t>(visit_row_index)];
    }

    const std::size_t state_offset = pool.stage1_visit_column_count;
    const int start_state_row_index = pool.entry_start_time_row_index[entry_index];
    if (start_state_row_index >= 0) {
        breakdown.state_contribution += dense_duals.stage1_input[
            state_offset + static_cast<std::size_t>(start_state_row_index)
        ];
    }
    const int last_state_row_index = pool.entry_last_time_row_index[entry_index];
    if (last_state_row_index >= 0) {
        breakdown.state_contribution -= dense_duals.stage1_input[
            state_offset + static_cast<std::size_t>(last_state_row_index)
        ];
    }

    const std::size_t edge_offset =
        pool.stage1_visit_column_count + pool.stage1_state_column_count;
    for (std::size_t k = pool.entry_edge_ids_row_ptr[entry_index];
         k < pool.entry_edge_ids_row_ptr[entry_index + 1U];
         ++k) {
        const int edge_id = pool.entry_edge_ids[k];
        if (edge_id < 0 ||
            static_cast<std::size_t>(edge_id) >= pool.stage1_edge_column_count) {
            continue;
        }
        breakdown.edge_contribution +=
            dense_duals.stage1_input[edge_offset + static_cast<std::size_t>(edge_id)];
    }

    breakdown.base_reduced_cost =
        breakdown.entry_cost_term - breakdown.visit_contribution -
        breakdown.state_contribution - breakdown.edge_contribution;
    return breakdown;
}

FullEnumerationBaseRCBreakdown compute_full_enumeration_entry_base_reduced_cost(
    const FullEnumerationStaticPool& pool,
    const FullEnumerationDenseDualVectors& dense_duals,
    CGPhase phase,
    std::size_t entry_index
) {
    FullEnumerationBaseRCBreakdown breakdown;
    breakdown.entry_cost_term =
        (phase == CGPhase::PhaseII) ? pool.entry_total_cost[entry_index] : 0.0;

    const std::size_t visit_limit = pool.stage1_visit_column_count;
    const std::size_t state_limit =
        pool.stage1_visit_column_count + pool.stage1_state_column_count;
    for (std::size_t k = pool.stage1_matrix.row_ptr[entry_index];
         k < pool.stage1_matrix.row_ptr[entry_index + 1U];
         ++k) {
        const std::size_t column_index =
            static_cast<std::size_t>(pool.stage1_matrix.column_index[k]);
        const double contribution =
            pool.stage1_matrix.value[k] * dense_duals.stage1_input[column_index];
        if (column_index < visit_limit) {
            breakdown.visit_contribution += contribution;
        } else if (column_index < state_limit) {
            breakdown.state_contribution += contribution;
        } else {
            breakdown.edge_contribution += contribution;
        }
    }

    breakdown.base_reduced_cost =
        breakdown.entry_cost_term - breakdown.visit_contribution -
        breakdown.state_contribution - breakdown.edge_contribution;
    return breakdown;
}

struct FullEnumerationVariantRCBreakdown {
    double base_reduced_cost = 0.0;
    double start_time_contribution = 0.0;
    double end_time_contribution = 0.0;
    double reduced_cost = 0.0;
};

FullEnumerationVariantRCBreakdown compute_full_enumeration_variant_reduced_cost(
    const FullEnumerationStaticPool& pool,
    const FullEnumerationDenseDualVectors& dense_duals,
    const std::vector<double>& base_reduced_costs,
    std::size_t variant_index
) {
    const std::size_t entry_index = full_enumeration_variant_entry_index(variant_index);
    const double tau = pool.variant_tau[variant_index];

    FullEnumerationVariantRCBreakdown breakdown;
    breakdown.base_reduced_cost = base_reduced_costs[entry_index];

    const int start_row_index = pool.entry_start_time_row_index[entry_index];
    if (start_row_index >= 0) {
        breakdown.start_time_contribution =
            -tau * dense_duals.time[static_cast<std::size_t>(start_row_index)];
    }

    const int last_row_index = pool.entry_last_time_row_index[entry_index];
    if (last_row_index >= 0) {
        breakdown.end_time_contribution =
            (tau + pool.entry_total_time[entry_index]) *
            dense_duals.time[static_cast<std::size_t>(last_row_index)];
    }

    breakdown.reduced_cost =
        breakdown.base_reduced_cost + breakdown.start_time_contribution +
        breakdown.end_time_contribution;
    return breakdown;
}

std::size_t resolve_full_enumeration_thread_count(
    std::size_t requested_thread_count,
    std::size_t item_count
) {
    if (item_count == 0U) {
        return 1U;
    }

    const std::size_t hardware_thread_count = std::max<std::size_t>(
        1U,
        static_cast<std::size_t>(std::thread::hardware_concurrency())
    );
    const std::size_t target =
        requested_thread_count == 0U ? hardware_thread_count : requested_thread_count;
    return std::max<std::size_t>(1U, std::min(target, item_count));
}

std::pair<std::size_t, std::size_t> balanced_work_range(
    std::size_t item_count,
    std::size_t worker_count,
    std::size_t worker_index
) {
    const std::size_t begin = (item_count * worker_index) / worker_count;
    const std::size_t end = (item_count * (worker_index + 1U)) / worker_count;
    return {begin, end};
}

bool full_pool_selection_less(
    const FullPoolSelection& lhs,
    const FullPoolSelection& rhs
) {
    if (!double_equal(lhs.reduced_cost, rhs.reduced_cost)) {
        return lhs.reduced_cost < rhs.reduced_cost;
    }
    if (lhs.entry_index != rhs.entry_index) {
        return lhs.entry_index < rhs.entry_index;
    }
    return lhs.tau_index < rhs.tau_index;
}

void maybe_push_full_pool_top_k(
    FullPoolTopKQueue& top_k_queue,
    const FullPoolSelection& selection,
    std::size_t max_total_columns
) {
    if (max_total_columns == 0U) {
        return;
    }
    if (top_k_queue.size() < max_total_columns) {
        top_k_queue.push(selection);
        return;
    }
    if (full_pool_selection_less(selection, top_k_queue.top())) {
        top_k_queue.pop();
        top_k_queue.push(selection);
    }
}

std::vector<FullPoolSelection> extract_full_pool_top_k(FullPoolTopKQueue& top_k_queue) {
    std::vector<FullPoolSelection> selected;
    selected.reserve(top_k_queue.size());
    while (!top_k_queue.empty()) {
        selected.push_back(top_k_queue.top());
        top_k_queue.pop();
    }
    return selected;
}

void log_full_enumeration_pool_selection_summary(
    std::ostream& out,
    const char* mode_name,
    std::size_t negative_column_count,
    std::size_t selected_column_count,
    double best_reduced_cost
) {
    out << "[full-enum-rc] mode=" << mode_name
        << " negative_columns=" << negative_column_count
        << " selected_columns=" << selected_column_count
        << " best_rc=" << format_double(best_reduced_cost)
        << '\n';
}

struct FullEnumerationParallelStage2ScanResult {
    std::vector<FullPoolSelection> negative_selections;
    std::size_t total_negative_column_count = 0;
    double best_reduced_cost = std::numeric_limits<double>::infinity();
};

void compute_full_enumeration_stage1_custom(
    const FullEnumerationStaticPool& pool,
    const FullEnumerationPoolNodeState& node_state,
    const FullEnumerationDenseDualVectors& dense_duals,
    CGPhase phase,
    std::size_t requested_thread_count,
    std::vector<double>& base_reduced_costs,
    bool rc_detail_log_enabled,
    std::ostream* log_stream
) {
    const std::size_t entry_thread_count =
        resolve_full_enumeration_thread_count(
            requested_thread_count,
            full_enumeration_entry_count(pool)
        );

    std::vector<std::thread> entry_threads;
    std::vector<std::ostringstream> entry_thread_logs(
        rc_detail_log_enabled ? entry_thread_count : 0U
    );
    entry_threads.reserve(entry_thread_count);
    for (std::size_t thread_index = 0; thread_index < entry_thread_count; ++thread_index) {
        entry_threads.emplace_back(
            [&, thread_index]() {
                const auto [begin, end] = balanced_work_range(
                    full_enumeration_entry_count(pool),
                    entry_thread_count,
                    thread_index
                );
                std::ostringstream* thread_log =
                    rc_detail_log_enabled ? &entry_thread_logs[thread_index] : nullptr;
                if (thread_log != nullptr) {
                    *thread_log
                        << "[full-enum-rc] mode=parallel stage=entry backend=custom thread="
                        << thread_index
                        << " assigned_entries=[" << begin << "," << end << ")\n";
                }

                for (std::size_t entry_index = begin; entry_index < end; ++entry_index) {
                    if (full_enumeration_entry_is_forbidden(node_state, entry_index) ||
                        node_state.entry_available_variant_count[entry_index] == 0U) {
                        base_reduced_costs[entry_index] = 0.0;
                        continue;
                    }
                    if (thread_log == nullptr) {
                        base_reduced_costs[entry_index] =
                            compute_full_enumeration_entry_base_reduced_cost_direct(
                                pool,
                                dense_duals,
                                phase,
                                entry_index
                            );
                        continue;
                    }

                    const FullEnumerationBaseRCBreakdown breakdown =
                        compute_full_enumeration_entry_base_reduced_cost_direct_breakdown(
                            pool,
                            dense_duals,
                            phase,
                            entry_index
                        );
                    base_reduced_costs[entry_index] = breakdown.base_reduced_cost;
                    if (node_state.entry_available_variant_count[entry_index] > 0U) {
                        *thread_log
                            << "[full-enum-rc] mode=parallel stage=entry backend=custom thread="
                            << thread_index
                            << " entry=" << entry_index
                            << " cost_term=" << format_double(breakdown.entry_cost_term)
                            << " visit_part=" << format_double(breakdown.visit_contribution)
                            << " state_part=" << format_double(breakdown.state_contribution)
                            << " edge_part=" << format_double(breakdown.edge_contribution)
                            << " base_rc=" << format_double(breakdown.base_reduced_cost)
                            << '\n';
                    }
                }
            }
        );
    }
    for (std::thread& worker : entry_threads) {
        worker.join();
    }
    if (rc_detail_log_enabled && log_stream != nullptr) {
        for (const std::ostringstream& thread_log : entry_thread_logs) {
            *log_stream << thread_log.str();
        }
    }
}

FullEnumerationParallelStage2ScanResult compute_full_enumeration_stage2_custom(
    const FullEnumerationStaticPool& pool,
    const FullEnumerationPoolNodeState& node_state,
    const FullEnumerationDenseDualVectors& dense_duals,
    const std::vector<double>& base_reduced_costs,
    double reduced_cost_tolerance,
    std::size_t max_total_columns,
    std::size_t requested_thread_count,
    bool rc_detail_log_enabled,
    std::ostream* log_stream
) {
    FullEnumerationParallelStage2ScanResult result;
    const std::size_t variant_thread_count =
        resolve_full_enumeration_thread_count(
            requested_thread_count,
            full_enumeration_variant_count(pool)
        );

    std::vector<std::thread> variant_threads;
    std::vector<std::ostringstream> variant_thread_logs(
        rc_detail_log_enabled ? variant_thread_count : 0U
    );
    std::vector<std::vector<FullPoolSelection>> thread_top_negative_selections(
        variant_thread_count
    );
    std::vector<std::size_t> thread_negative_counts(variant_thread_count, 0U);
    std::vector<double> thread_best_reduced_costs(
        variant_thread_count,
        std::numeric_limits<double>::infinity()
    );
    variant_threads.reserve(variant_thread_count);
    for (std::size_t thread_index = 0; thread_index < variant_thread_count; ++thread_index) {
        variant_threads.emplace_back(
            [&, thread_index]() {
                const auto [begin, end] = balanced_work_range(
                    full_enumeration_variant_count(pool),
                    variant_thread_count,
                    thread_index
                );
                std::ostringstream* thread_log =
                    rc_detail_log_enabled ? &variant_thread_logs[thread_index] : nullptr;
                if (thread_log != nullptr) {
                    *thread_log
                        << "[full-enum-rc] mode=parallel stage=variant backend=custom thread="
                        << thread_index
                        << " assigned_variants=[" << begin << "," << end << ")\n";
                }

                FullPoolTopKQueue local_top_negative_columns;
                for (std::size_t variant_index = begin; variant_index < end; ++variant_index) {
                    if (!full_enumeration_variant_is_available(pool, node_state, variant_index)) {
                        continue;
                    }
                    const std::size_t entry_index =
                        full_enumeration_variant_entry_index(variant_index);
                    const FullEnumerationVariantRCBreakdown breakdown =
                        compute_full_enumeration_variant_reduced_cost(
                            pool,
                            dense_duals,
                            base_reduced_costs,
                            variant_index
                        );
                    thread_best_reduced_costs[thread_index] = std::min(
                        thread_best_reduced_costs[thread_index],
                        breakdown.reduced_cost
                    );
                    if (thread_log != nullptr) {
                        *thread_log
                            << "[full-enum-rc] mode=parallel stage=variant backend=custom thread="
                            << thread_index
                            << " entry=" << entry_index
                            << " tau_index="
                            << full_enumeration_variant_tau_index(variant_index)
                            << " tau=" << format_double(pool.variant_tau[variant_index])
                            << " base_rc=" << format_double(breakdown.base_reduced_cost)
                            << " start_time_part="
                            << format_double(breakdown.start_time_contribution)
                            << " end_time_part="
                            << format_double(breakdown.end_time_contribution)
                            << " rc=" << format_double(breakdown.reduced_cost)
                            << " negative="
                            << (breakdown.reduced_cost < reduced_cost_tolerance ? 1 : 0)
                            << '\n';
                    }
                    if (breakdown.reduced_cost < reduced_cost_tolerance) {
                        ++thread_negative_counts[thread_index];
                        maybe_push_full_pool_top_k(
                            local_top_negative_columns,
                            {
                            breakdown.reduced_cost,
                            variant_index,
                            entry_index,
                            full_enumeration_variant_tau_index(variant_index),
                            },
                            max_total_columns
                        );
                    }
                }
                thread_top_negative_selections[thread_index] =
                    extract_full_pool_top_k(local_top_negative_columns);
            }
        );
    }
    for (std::thread& worker : variant_threads) {
        worker.join();
    }
    if (rc_detail_log_enabled && log_stream != nullptr) {
        for (const std::ostringstream& thread_log : variant_thread_logs) {
            *log_stream << thread_log.str();
        }
    }

    FullPoolTopKQueue global_top_negative_columns;
    for (const std::vector<FullPoolSelection>& local_top_negative_selections :
         thread_top_negative_selections) {
        for (const FullPoolSelection& selection : local_top_negative_selections) {
            maybe_push_full_pool_top_k(global_top_negative_columns, selection, max_total_columns);
        }
    }
    result.negative_selections = extract_full_pool_top_k(global_top_negative_columns);
    for (std::size_t thread_negative_count : thread_negative_counts) {
        result.total_negative_column_count += thread_negative_count;
    }
    for (double thread_best_reduced_cost : thread_best_reduced_costs) {
        result.best_reduced_cost = std::min(result.best_reduced_cost, thread_best_reduced_cost);
    }
    return result;
}

#ifdef SPDP_HAVE_ONEMKL

MKL_INT to_mkl_int(
    std::size_t value,
    const char* context
) {
    const std::size_t max_value =
        static_cast<std::size_t>(std::numeric_limits<MKL_INT>::max());
    if (value > max_value) {
        throw std::runtime_error(
            std::string("oneMKL index overflow in ") + context
        );
    }
    return static_cast<MKL_INT>(value);
}

void check_mkl_status(
    sparse_status_t status,
    const char* context
) {
    if (status != SPARSE_STATUS_SUCCESS) {
        throw std::runtime_error(
            std::string("oneMKL sparse operation failed in ") + context +
            " (status=" + std::to_string(static_cast<int>(status)) + ")"
        );
    }
}

class ScopedMklThreads {
public:
    explicit ScopedMklThreads(std::size_t thread_count)
        : previous_thread_count_(
              mkl_set_num_threads_local(static_cast<int>(std::max<std::size_t>(1U, thread_count)))
          ) {}

    ~ScopedMklThreads() {
        mkl_set_num_threads_local(previous_thread_count_);
    }

private:
    int previous_thread_count_ = 1;
};

class ScopedMklSparseMatrix {
public:
    ~ScopedMklSparseMatrix() {
        if (handle_ != nullptr) {
            mkl_sparse_destroy(handle_);
        }
    }

    sparse_matrix_t* out() {
        return &handle_;
    }

    sparse_matrix_t get() const {
        return handle_;
    }

private:
    sparse_matrix_t handle_ = nullptr;
};

struct FullEnumerationMklCsrMatrix {
    MKL_INT row_count = 0;
    MKL_INT column_count = 0;
    std::vector<MKL_INT> row_start;
    std::vector<MKL_INT> row_end;
    std::vector<MKL_INT> column_index;
    std::vector<double> values;
};

FullEnumerationMklCsrMatrix build_full_enumeration_mkl_matrix(
    const FullEnumerationStaticPool::CSRMatrix& csr_matrix,
    const char* context
) {
    FullEnumerationMklCsrMatrix mkl_matrix;
    const std::size_t row_count =
        csr_matrix.row_ptr.empty() ? 0U : (csr_matrix.row_ptr.size() - 1U);
    mkl_matrix.row_count = to_mkl_int(row_count, context);
    mkl_matrix.column_count = to_mkl_int(csr_matrix.column_count, context);
    mkl_matrix.row_start.reserve(row_count);
    mkl_matrix.row_end.reserve(row_count);
    for (std::size_t row = 0; row < row_count; ++row) {
        mkl_matrix.row_start.push_back(to_mkl_int(csr_matrix.row_ptr[row], context));
        mkl_matrix.row_end.push_back(to_mkl_int(csr_matrix.row_ptr[row + 1U], context));
    }
    mkl_matrix.column_index.reserve(csr_matrix.column_index.size());
    mkl_matrix.values.reserve(csr_matrix.value.size());
    for (std::size_t idx = 0; idx < csr_matrix.column_index.size(); ++idx) {
        mkl_matrix.column_index.push_back(
            to_mkl_int(static_cast<std::size_t>(csr_matrix.column_index[idx]), context)
        );
        mkl_matrix.values.push_back(csr_matrix.value[idx]);
    }
    return mkl_matrix;
}

std::vector<double> build_full_enumeration_stage2_mkl_input(
    const FullEnumerationStaticPool& pool,
    const FullEnumerationDenseDualVectors& dense_duals,
    const std::vector<double>& base_reduced_costs
) {
    std::vector<double> input;
    input.reserve(full_enumeration_entry_count(pool) + dense_duals.time.size());
    input.insert(input.end(), base_reduced_costs.begin(), base_reduced_costs.end());
    input.insert(input.end(), dense_duals.time.begin(), dense_duals.time.end());
    return input;
}

std::vector<double> run_full_enumeration_mkl_spmv(
    const FullEnumerationMklCsrMatrix& matrix,
    const std::vector<double>& input,
    std::size_t thread_count,
    const char* context
) {
    if (matrix.row_count == 0) {
        return {};
    }

    ScopedMklSparseMatrix handle;
    check_mkl_status(
        mkl_sparse_d_create_csr(
            handle.out(),
            SPARSE_INDEX_BASE_ZERO,
            matrix.row_count,
            matrix.column_count,
            const_cast<MKL_INT*>(matrix.row_start.data()),
            const_cast<MKL_INT*>(matrix.row_end.data()),
            const_cast<MKL_INT*>(matrix.column_index.data()),
            const_cast<double*>(matrix.values.data())
        ),
        context
    );

    matrix_descr descriptor;
    descriptor.type = SPARSE_MATRIX_TYPE_GENERAL;
    descriptor.mode = SPARSE_FILL_MODE_FULL;
    descriptor.diag = SPARSE_DIAG_NON_UNIT;

    check_mkl_status(mkl_sparse_optimize(handle.get()), context);

    std::vector<double> output(static_cast<std::size_t>(matrix.row_count), 0.0);
    ScopedMklThreads thread_guard(thread_count);
    check_mkl_status(
        mkl_sparse_d_mv(
            SPARSE_OPERATION_NON_TRANSPOSE,
            1.0,
            handle.get(),
            descriptor,
            input.data(),
            0.0,
            output.data()
        ),
        context
    );
    return output;
}

void compute_full_enumeration_stage1_onemkl(
    const FullEnumerationStaticPool& pool,
    const FullEnumerationPoolNodeState& node_state,
    const FullEnumerationDenseDualVectors& dense_duals,
    CGPhase phase,
    std::size_t requested_thread_count,
    std::vector<double>& base_reduced_costs,
    bool rc_detail_log_enabled,
    std::ostream* log_stream
) {
    const std::size_t thread_count =
        resolve_full_enumeration_thread_count(
            requested_thread_count,
            full_enumeration_entry_count(pool)
        );
    const FullEnumerationMklCsrMatrix matrix =
        build_full_enumeration_mkl_matrix(pool.stage1_matrix, "stage1 onemkl matrix");
    const std::vector<double> output =
        run_full_enumeration_mkl_spmv(
            matrix,
            dense_duals.stage1_input,
            thread_count,
            "stage1 onemkl SpMV"
        );

    if (rc_detail_log_enabled && log_stream != nullptr) {
        *log_stream << "[full-enum-rc] mode=parallel stage=entry backend=onemkl rows="
                    << full_enumeration_entry_count(pool)
                    << " cols=" << dense_duals.stage1_input.size()
                    << " threads=" << thread_count << '\n';
    }

    for (std::size_t entry_index = 0; entry_index < full_enumeration_entry_count(pool);
         ++entry_index) {
        const double entry_cost_term =
            (phase == CGPhase::PhaseII) ? pool.entry_total_cost[entry_index] : 0.0;
        base_reduced_costs[entry_index] = entry_cost_term - output[entry_index];

        if (rc_detail_log_enabled &&
            log_stream != nullptr &&
            !full_enumeration_entry_is_forbidden(node_state, entry_index) &&
            node_state.entry_available_variant_count[entry_index] > 0U) {
            const FullEnumerationBaseRCBreakdown breakdown =
                compute_full_enumeration_entry_base_reduced_cost(
                    pool,
                    dense_duals,
                    phase,
                    entry_index
                );
            *log_stream << "[full-enum-rc] mode=parallel stage=entry backend=onemkl entry="
                        << entry_index
                        << " cost_term=" << format_double(breakdown.entry_cost_term)
                        << " visit_part=" << format_double(breakdown.visit_contribution)
                        << " state_part=" << format_double(breakdown.state_contribution)
                        << " edge_part=" << format_double(breakdown.edge_contribution)
                        << " base_rc=" << format_double(base_reduced_costs[entry_index])
                        << '\n';
        }
    }
}

FullEnumerationParallelStage2ScanResult compute_full_enumeration_stage2_onemkl(
    const FullEnumerationStaticPool& pool,
    const FullEnumerationPoolNodeState& node_state,
    const FullEnumerationDenseDualVectors& dense_duals,
    const std::vector<double>& base_reduced_costs,
    double reduced_cost_tolerance,
    std::size_t max_total_columns,
    std::size_t requested_thread_count,
    bool rc_detail_log_enabled,
    std::ostream* log_stream
) {
    FullEnumerationParallelStage2ScanResult result;
    const std::size_t thread_count =
        resolve_full_enumeration_thread_count(
            requested_thread_count,
            full_enumeration_variant_count(pool)
        );
    const FullEnumerationMklCsrMatrix matrix =
        build_full_enumeration_mkl_matrix(pool.stage2_matrix, "stage2 onemkl matrix");
    const std::vector<double> input =
        build_full_enumeration_stage2_mkl_input(pool, dense_duals, base_reduced_costs);
    const std::vector<double> output =
        run_full_enumeration_mkl_spmv(matrix, input, thread_count, "stage2 onemkl SpMV");

    if (rc_detail_log_enabled && log_stream != nullptr) {
        *log_stream << "[full-enum-rc] mode=parallel stage=variant backend=onemkl rows="
                    << full_enumeration_variant_count(pool)
                    << " cols=" << input.size()
                    << " threads=" << thread_count << '\n';
    }

    std::vector<std::thread> variant_threads;
    std::vector<std::ostringstream> variant_thread_logs(
        rc_detail_log_enabled ? thread_count : 0U
    );
    std::vector<std::vector<FullPoolSelection>> thread_top_negative_selections(thread_count);
    std::vector<std::size_t> thread_negative_counts(thread_count, 0U);
    std::vector<double> thread_best_reduced_costs(
        thread_count,
        std::numeric_limits<double>::infinity()
    );
    variant_threads.reserve(thread_count);
    for (std::size_t thread_index = 0; thread_index < thread_count; ++thread_index) {
        variant_threads.emplace_back(
            [&, thread_index]() {
                const auto [begin, end] = balanced_work_range(
                    full_enumeration_variant_count(pool),
                    thread_count,
                    thread_index
                );
                std::ostringstream* thread_log =
                    rc_detail_log_enabled ? &variant_thread_logs[thread_index] : nullptr;
                if (thread_log != nullptr) {
                    *thread_log
                        << "[full-enum-rc] mode=parallel stage=variant backend=onemkl thread="
                        << thread_index
                        << " assigned_variants=[" << begin << "," << end << ")\n";
                }

                FullPoolTopKQueue local_top_negative_columns;
                for (std::size_t variant_index = begin; variant_index < end; ++variant_index) {
                    if (!full_enumeration_variant_is_available(pool, node_state, variant_index)) {
                        continue;
                    }
                    const std::size_t entry_index =
                        full_enumeration_variant_entry_index(variant_index);
                    const int start_row_index = pool.entry_start_time_row_index[entry_index];
                    const int last_row_index = pool.entry_last_time_row_index[entry_index];
                    const double start_time_contribution =
                        start_row_index >= 0
                            ? -pool.variant_tau[variant_index] *
                                  dense_duals.time[static_cast<std::size_t>(start_row_index)]
                            : 0.0;
                    const double end_time_contribution =
                        last_row_index >= 0
                            ? (pool.variant_tau[variant_index] +
                               pool.entry_total_time[entry_index]) *
                                  dense_duals.time[static_cast<std::size_t>(last_row_index)]
                            : 0.0;
                    const double reduced_cost = output[variant_index];
                    thread_best_reduced_costs[thread_index] = std::min(
                        thread_best_reduced_costs[thread_index],
                        reduced_cost
                    );

                    if (thread_log != nullptr) {
                        *thread_log
                            << "[full-enum-rc] mode=parallel stage=variant backend=onemkl thread="
                            << thread_index
                            << " entry=" << entry_index
                            << " tau_index="
                            << full_enumeration_variant_tau_index(variant_index)
                            << " tau=" << format_double(pool.variant_tau[variant_index])
                            << " base_rc=" << format_double(base_reduced_costs[entry_index])
                            << " start_time_part=" << format_double(start_time_contribution)
                            << " end_time_part=" << format_double(end_time_contribution)
                            << " rc=" << format_double(reduced_cost)
                            << " negative=" << (reduced_cost < reduced_cost_tolerance ? 1 : 0)
                            << '\n';
                    }

                    if (reduced_cost < reduced_cost_tolerance) {
                        ++thread_negative_counts[thread_index];
                        maybe_push_full_pool_top_k(
                            local_top_negative_columns,
                            {
                                reduced_cost,
                                variant_index,
                                entry_index,
                                full_enumeration_variant_tau_index(variant_index),
                            },
                            max_total_columns
                        );
                    }
                }
                thread_top_negative_selections[thread_index] =
                    extract_full_pool_top_k(local_top_negative_columns);
            }
        );
    }
    for (std::thread& worker : variant_threads) {
        worker.join();
    }
    if (rc_detail_log_enabled && log_stream != nullptr) {
        for (const std::ostringstream& thread_log : variant_thread_logs) {
            *log_stream << thread_log.str();
        }
    }

    FullPoolTopKQueue global_top_negative_columns;
    for (const std::vector<FullPoolSelection>& local_top_negative_selections :
         thread_top_negative_selections) {
        for (const FullPoolSelection& selection : local_top_negative_selections) {
            maybe_push_full_pool_top_k(global_top_negative_columns, selection, max_total_columns);
        }
    }
    result.negative_selections = extract_full_pool_top_k(global_top_negative_columns);
    for (std::size_t thread_negative_count : thread_negative_counts) {
        result.total_negative_column_count += thread_negative_count;
    }
    for (double thread_best_reduced_cost : thread_best_reduced_costs) {
        result.best_reduced_cost = std::min(result.best_reduced_cost, thread_best_reduced_cost);
    }
    return result;
}

#else

void compute_full_enumeration_stage1_onemkl(
    const FullEnumerationStaticPool&,
    const FullEnumerationPoolNodeState&,
    const FullEnumerationDenseDualVectors&,
    CGPhase,
    std::size_t,
    std::vector<double>&,
    bool,
    std::ostream*
) {
    throw std::runtime_error(
        "parallel stage1 backend=onemkl requested, but this build does not include oneMKL."
    );
}

FullEnumerationParallelStage2ScanResult compute_full_enumeration_stage2_onemkl(
    const FullEnumerationStaticPool&,
    const FullEnumerationPoolNodeState&,
    const FullEnumerationDenseDualVectors&,
    const std::vector<double>&,
    double,
    std::size_t,
    std::size_t,
    bool,
    std::ostream*
) {
    throw std::runtime_error(
        "parallel stage2 backend=onemkl requested, but this build does not include oneMKL."
    );
}

#endif

}  // namespace

FullEnumerationStaticPool build_full_enumeration_static_pool(
    const MultiDiGraph& graph,
    const FullEnumerationPoolBuildOptions& options,
    FullEnumerationPoolBuildStats* stats
) {
    const auto start_time = std::chrono::steady_clock::now();

    CompactPStepOptions compact_options;
    compact_options.p = options.p;
    compact_options.time_limit = options.time_limit;
    compact_options.validate = false;
    compact_options.prune_pickup_symmetry_43 = options.prune_pickup_symmetry_43;
    compact_options.prune_delivery_symmetry_43 = options.prune_delivery_symmetry_43;

    const CompactPStepArtifacts artifacts =
        build_compact_pstep_artifacts(graph, compact_options, nullptr);

    FullEnumerationStaticPool pool;
    pool.time_limit = options.time_limit;
    pool.graph_edge_source_node_id.reserve(graph.number_of_edges());
    pool.graph_edge_target_node_id.reserve(graph.number_of_edges());
    for (const EdgeRecord& edge : graph.edges()) {
        pool.graph_edge_source_node_id.push_back(edge.u);
        pool.graph_edge_target_node_id.push_back(edge.v);
    }
    std::map<int, std::size_t> entry_index_by_raw_path_id;
    std::size_t skipped_duplicate_columns = 0;
    std::vector<FullEnumerationEntryMatrixBuildData> entry_matrix_build_data;
    entry_matrix_build_data.reserve(artifacts.raw_paths.size());
    pool.entry_edge_ids_row_ptr.push_back(0U);
    pool.entry_node_sequence_row_ptr.push_back(0U);
    pool.entry_state_sequence_row_ptr.push_back(0U);

    for (const CompactPStep& pstep : artifacts.compact_psteps) {
        const std::string key = build_column_key(pstep.edge_ids, pstep.tau);
        if (pool.variant_index_by_key.find(key) != pool.variant_index_by_key.end()) {
            ++skipped_duplicate_columns;
            continue;
        }

        std::size_t entry_index = 0;
        const auto found_entry = entry_index_by_raw_path_id.find(pstep.raw_path_id);
        if (found_entry == entry_index_by_raw_path_id.end()) {
            entry_index = full_enumeration_entry_count(pool);
            append_full_enumeration_entry(
                pool,
                entry_matrix_build_data,
                graph,
                pstep,
                artifacts.coefficients
            );
            entry_index_by_raw_path_id.emplace(pstep.raw_path_id, entry_index);
        } else {
            entry_index = found_entry->second;
        }

        const double latest_start = full_enumeration_entry_latest_start(pool, entry_index);
        std::size_t tau_index = 0U;
        if (double_equal(pstep.tau, 0.0)) {
            tau_index = 0U;
        } else if (double_equal(pstep.tau, latest_start)) {
            tau_index = 1U;
        } else {
            throw std::runtime_error(
                "Full-enumeration pool encountered a tau value that is neither 0 nor T-t(entry)."
            );
        }

        const std::size_t variant_index =
            full_enumeration_variant_index(entry_index, tau_index);
        if (pool.variant_valid[variant_index] != 0U) {
            throw std::runtime_error(
                "Full-enumeration pool encountered duplicate tau endpoint for one entry."
            );
        }
        pool.variant_tau[variant_index] = pstep.tau;
        pool.variant_valid[variant_index] = 1U;
        pool.variant_index_by_key.emplace(key, variant_index);
    }

    finalize_full_enumeration_stage1_matrix(
        pool,
        entry_matrix_build_data,
        graph.number_of_edges()
    );
    finalize_full_enumeration_stage2_matrix(pool);
    finalize_full_enumeration_direct_visit_lookup(pool);

    if (stats != nullptr) {
        stats->raw_path_count = artifacts.raw_paths.size();
        stats->compact_pstep_count = artifacts.compact_psteps.size();
        stats->inactive_column_count = 0U;
        stats->inactive_path_count = 0U;
        stats->skipped_master_column_count = 0U;
        stats->skipped_duplicate_column_count = skipped_duplicate_columns;
        stats->runtime_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();
    }

    return pool;
}

FullEnumerationPoolNodeState build_full_enumeration_pool_node_state(
    const FullEnumerationStaticPool& static_pool,
    const std::map<std::string, int>& active_column_id_by_key,
    FullEnumerationPoolBuildStats* stats
) {
    FullEnumerationPoolNodeState node_state;
    node_state.variant_active.assign(full_enumeration_variant_count(static_pool), 0U);
    node_state.forbidden_entry.assign(full_enumeration_entry_count(static_pool), 0U);
    node_state.entry_available_variant_count.assign(
        full_enumeration_entry_count(static_pool),
        0U
    );

    for (std::size_t variant_index = 0; variant_index < full_enumeration_variant_count(static_pool);
         ++variant_index) {
        if (!full_enumeration_variant_is_valid(static_pool, variant_index)) {
            continue;
        }
        const std::size_t entry_index = full_enumeration_variant_entry_index(variant_index);
        if (node_state.entry_available_variant_count[entry_index] == 0U) {
            ++node_state.available_entry_count;
        }
        ++node_state.entry_available_variant_count[entry_index];
        ++node_state.available_variant_count;
    }

    std::size_t skipped_master_columns = 0U;
    for (const auto& [key, _] : active_column_id_by_key) {
        const auto found = static_pool.variant_index_by_key.find(key);
        if (found == static_pool.variant_index_by_key.end()) {
            continue;
        }
        const std::size_t variant_index = found->second;
        if (!full_enumeration_variant_is_valid(static_pool, variant_index) ||
            node_state.variant_active[variant_index] != 0U) {
            continue;
        }
        node_state.variant_active[variant_index] = 1U;
        if (node_state.available_variant_count > 0U) {
            --node_state.available_variant_count;
        }
        const std::size_t entry_index = full_enumeration_variant_entry_index(variant_index);
        if (node_state.entry_available_variant_count[entry_index] > 0U) {
            --node_state.entry_available_variant_count[entry_index];
            if (node_state.entry_available_variant_count[entry_index] == 0U &&
                node_state.available_entry_count > 0U) {
                --node_state.available_entry_count;
            }
        }
        ++skipped_master_columns;
    }

    if (stats != nullptr) {
        stats->inactive_column_count = node_state.available_variant_count;
        stats->inactive_path_count = node_state.available_entry_count;
        stats->skipped_master_column_count = skipped_master_columns;
    }

    return node_state;
}

namespace {

ForwardPricingResult run_full_enumeration_pool_pricing_sequential(
    const FullEnumerationStaticPool& pool,
    const FullEnumerationPoolNodeState& node_state,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    double reduced_cost_tolerance,
    std::size_t max_total_columns,
    bool rc_detail_log_enabled,
    std::ostream* log_stream
) {
    ForwardPricingResult result;
    result.status = ForwardPricingStatus::NoNegativeColumn;
    result.best_reduced_cost = std::numeric_limits<double>::infinity();

    const FullEnumerationDenseDualVectors dense_duals =
        build_full_enumeration_dense_duals(pool, dual_solution);
    std::vector<double> base_reduced_costs(full_enumeration_entry_count(pool), 0.0);

    if (log_stream != nullptr) {
        *log_stream << "[full-enum-rc] mode=sequential inactive_entries="
                    << node_state.available_entry_count
                    << " inactive_variants=" << node_state.available_variant_count
                    << " max_total_columns=" << max_total_columns << '\n';
    }

    std::priority_queue<
        FullPoolSelection,
        std::vector<FullPoolSelection>,
        FullPoolWorstReducedCostFirst>
        top_negative_columns;

    for (std::size_t entry_index = 0; entry_index < full_enumeration_entry_count(pool);
         ++entry_index) {
        if (full_enumeration_entry_is_forbidden(node_state, entry_index) ||
            node_state.entry_available_variant_count[entry_index] == 0U) {
            continue;
        }
        const FullEnumerationBaseRCBreakdown breakdown =
            compute_full_enumeration_entry_base_reduced_cost(
                pool,
                dense_duals,
                phase,
                entry_index
            );
        base_reduced_costs[entry_index] = breakdown.base_reduced_cost;
        if (rc_detail_log_enabled &&
            log_stream != nullptr &&
            node_state.entry_available_variant_count[entry_index] > 0U) {
            *log_stream << "[full-enum-rc] mode=sequential stage=entry entry=" << entry_index
                        << " cost_term=" << format_double(breakdown.entry_cost_term)
                        << " visit_part=" << format_double(breakdown.visit_contribution)
                        << " state_part=" << format_double(breakdown.state_contribution)
                        << " edge_part=" << format_double(breakdown.edge_contribution)
                        << " base_rc=" << format_double(breakdown.base_reduced_cost)
                        << '\n';
        }
    }

    for (std::size_t variant_index = 0; variant_index < full_enumeration_variant_count(pool);
         ++variant_index) {
        if (!full_enumeration_variant_is_available(pool, node_state, variant_index)) {
            continue;
        }
        const std::size_t entry_index = full_enumeration_variant_entry_index(variant_index);
        const FullEnumerationVariantRCBreakdown breakdown =
            compute_full_enumeration_variant_reduced_cost(
                pool,
                dense_duals,
                base_reduced_costs,
                variant_index
            );
        result.best_reduced_cost = std::min(result.best_reduced_cost, breakdown.reduced_cost);
        if (rc_detail_log_enabled && log_stream != nullptr) {
            *log_stream << "[full-enum-rc] mode=sequential stage=variant entry="
                        << entry_index
                        << " tau_index=" << full_enumeration_variant_tau_index(variant_index)
                        << " tau=" << format_double(pool.variant_tau[variant_index])
                        << " base_rc=" << format_double(breakdown.base_reduced_cost)
                        << " start_time_part="
                        << format_double(breakdown.start_time_contribution)
                        << " end_time_part="
                        << format_double(breakdown.end_time_contribution)
                        << " rc=" << format_double(breakdown.reduced_cost)
                        << " negative="
                        << (breakdown.reduced_cost < reduced_cost_tolerance ? 1 : 0)
                        << '\n';
        }

        if (breakdown.reduced_cost < reduced_cost_tolerance) {
            ++result.total_negative_column_count;
            FullPoolSelection selection{
                breakdown.reduced_cost,
                variant_index,
                entry_index,
                full_enumeration_variant_tau_index(variant_index),
            };
            maybe_push_full_pool_top_k(top_negative_columns, selection, max_total_columns);
        }
    }

    std::vector<FullPoolSelection> selected = extract_full_pool_top_k(top_negative_columns);
    std::sort(selected.begin(), selected.end(), full_pool_selection_less);

    for (const FullPoolSelection& selection : selected) {
        result.columns.push_back(build_column_from_full_pool_selection(pool, selection));
    }

    if (!std::isfinite(result.best_reduced_cost)) {
        result.best_reduced_cost = 0.0;
    }
    result.status = result.columns.empty()
                        ? ForwardPricingStatus::NoNegativeColumn
                        : ForwardPricingStatus::ColumnsFound;

    if (log_stream != nullptr) {
        log_full_enumeration_pool_selection_summary(
            *log_stream,
            "sequential",
            result.total_negative_column_count,
            result.columns.size(),
            result.best_reduced_cost
        );
    }
    return result;
}

ForwardPricingResult run_full_enumeration_pool_pricing_parallel(
    const FullEnumerationStaticPool& pool,
    const FullEnumerationPoolNodeState& node_state,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    double reduced_cost_tolerance,
    std::size_t max_total_columns,
    FullEnumerationRCUpdateBackend parallel_stage1_backend,
    FullEnumerationRCUpdateBackend parallel_stage2_backend,
    std::size_t requested_thread_count,
    bool rc_detail_log_enabled,
    std::ostream* log_stream
) {
    ForwardPricingResult result;
    result.status = ForwardPricingStatus::NoNegativeColumn;
    result.best_reduced_cost = std::numeric_limits<double>::infinity();

    const FullEnumerationDenseDualVectors dense_duals =
        build_full_enumeration_dense_duals(pool, dual_solution);
    std::vector<double> base_reduced_costs(full_enumeration_entry_count(pool), 0.0);

    const std::size_t entry_thread_count =
        resolve_full_enumeration_thread_count(
            requested_thread_count,
            full_enumeration_entry_count(pool)
        );
    const std::size_t variant_thread_count =
        resolve_full_enumeration_thread_count(
            requested_thread_count,
            full_enumeration_variant_count(pool)
        );

    if (log_stream != nullptr) {
        *log_stream << "[full-enum-rc] mode=parallel inactive_entries="
                    << node_state.available_entry_count
                    << " inactive_variants=" << node_state.available_variant_count
                    << " requested_threads=" << requested_thread_count
                    << " entry_threads=" << entry_thread_count
                    << " variant_threads=" << variant_thread_count
                    << " max_total_columns=" << max_total_columns
                    << " stage1_backend=" << to_string(parallel_stage1_backend)
                    << " stage2_backend=" << to_string(parallel_stage2_backend)
                    << " selection=top-k\n";
    }

    switch (parallel_stage1_backend) {
        case FullEnumerationRCUpdateBackend::Custom:
            compute_full_enumeration_stage1_custom(
                pool,
                node_state,
                dense_duals,
                phase,
                requested_thread_count,
                base_reduced_costs,
                rc_detail_log_enabled,
                log_stream
            );
            break;
        case FullEnumerationRCUpdateBackend::OneMKL:
            compute_full_enumeration_stage1_onemkl(
                pool,
                node_state,
                dense_duals,
                phase,
                requested_thread_count,
                base_reduced_costs,
                rc_detail_log_enabled,
                log_stream
            );
            break;
    }

    FullEnumerationParallelStage2ScanResult stage2_result;
    switch (parallel_stage2_backend) {
        case FullEnumerationRCUpdateBackend::Custom:
            stage2_result = compute_full_enumeration_stage2_custom(
                pool,
                node_state,
                dense_duals,
                base_reduced_costs,
                reduced_cost_tolerance,
                max_total_columns,
                requested_thread_count,
                rc_detail_log_enabled,
                log_stream
            );
            break;
        case FullEnumerationRCUpdateBackend::OneMKL:
            stage2_result = compute_full_enumeration_stage2_onemkl(
                pool,
                node_state,
                dense_duals,
                base_reduced_costs,
                reduced_cost_tolerance,
                max_total_columns,
                requested_thread_count,
                rc_detail_log_enabled,
                log_stream
            );
            break;
    }

    std::vector<FullPoolSelection> selected = std::move(stage2_result.negative_selections);
    result.total_negative_column_count = stage2_result.total_negative_column_count;
    result.best_reduced_cost = stage2_result.best_reduced_cost;
    std::sort(selected.begin(), selected.end(), full_pool_selection_less);

    result.columns.reserve(selected.size());
    for (const FullPoolSelection& selection : selected) {
        result.columns.push_back(build_column_from_full_pool_selection(pool, selection));
    }

    if (!std::isfinite(result.best_reduced_cost)) {
        result.best_reduced_cost = 0.0;
    }
    result.status = result.columns.empty()
                        ? ForwardPricingStatus::NoNegativeColumn
                        : ForwardPricingStatus::ColumnsFound;

    if (log_stream != nullptr) {
        log_full_enumeration_pool_selection_summary(
            *log_stream,
            "parallel",
            result.total_negative_column_count,
            result.columns.size(),
            result.best_reduced_cost
        );
    }
    return result;
}

}  // namespace

ForwardPricingResult run_full_enumeration_pool_pricing(
    const FullEnumerationStaticPool& pool,
    const FullEnumerationPoolNodeState& node_state,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    double reduced_cost_tolerance,
    std::size_t max_total_columns,
    FullEnumerationRCUpdateMode rc_update_mode,
    FullEnumerationRCUpdateBackend parallel_stage1_backend,
    FullEnumerationRCUpdateBackend parallel_stage2_backend,
    std::size_t requested_thread_count,
    bool rc_detail_log_enabled,
    std::ostream* log_stream
) {
    switch (rc_update_mode) {
        case FullEnumerationRCUpdateMode::Sequential:
            return run_full_enumeration_pool_pricing_sequential(
                pool,
                node_state,
                dual_solution,
                phase,
                reduced_cost_tolerance,
                max_total_columns,
                rc_detail_log_enabled,
                log_stream
            );
        case FullEnumerationRCUpdateMode::Parallel:
            return run_full_enumeration_pool_pricing_parallel(
                pool,
                node_state,
                dual_solution,
                phase,
                reduced_cost_tolerance,
                max_total_columns,
                parallel_stage1_backend,
                parallel_stage2_backend,
                requested_thread_count,
                rc_detail_log_enabled,
                log_stream
            );
    }
    throw std::runtime_error("Unsupported full-enumeration RC update mode.");
}

void mark_full_enumeration_pool_columns_active(
    const FullEnumerationStaticPool& pool,
    FullEnumerationPoolNodeState& node_state,
    const std::vector<CGColumn>& columns
) {
    for (const CGColumn& column : columns) {
        std::size_t variant_index = std::numeric_limits<std::size_t>::max();
        if (column.pool_variant_index >= 0 &&
            static_cast<std::size_t>(column.pool_variant_index) < full_enumeration_variant_count(pool)) {
            variant_index = static_cast<std::size_t>(column.pool_variant_index);
        } else {
            const std::string key = build_column_key(column.edge_ids, column.tau);
            const auto found = pool.variant_index_by_key.find(key);
            if (found == pool.variant_index_by_key.end()) {
                continue;
            }
            variant_index = found->second;
        }

        if (!full_enumeration_variant_is_available(pool, node_state, variant_index)) {
            continue;
        }

        node_state.variant_active[variant_index] = 1U;
        if (node_state.available_variant_count > 0U) {
            --node_state.available_variant_count;
        }
        const std::size_t entry_index = full_enumeration_variant_entry_index(variant_index);
        if (node_state.entry_available_variant_count[entry_index] > 0U) {
            --node_state.entry_available_variant_count[entry_index];
            if (node_state.entry_available_variant_count[entry_index] == 0U &&
                node_state.available_entry_count > 0U) {
                --node_state.available_entry_count;
            }
        }
    }
}

void mark_full_enumeration_pool_entries_forbidden_by_edge_mask(
    const FullEnumerationStaticPool& static_pool,
    const std::vector<std::uint8_t>& forbidden_edge_mask,
    FullEnumerationPoolNodeState& node_state
) {
    if (forbidden_edge_mask.empty()) {
        return;
    }

    for (std::size_t entry_index = 0; entry_index < full_enumeration_entry_count(static_pool);
         ++entry_index) {
        if (full_enumeration_entry_is_forbidden(node_state, entry_index)) {
            continue;
        }

        bool contains_forbidden_edge = false;
        for (std::size_t k = static_pool.entry_edge_ids_row_ptr[entry_index];
             k < static_pool.entry_edge_ids_row_ptr[entry_index + 1U];
             ++k) {
            const int edge_id = static_pool.entry_edge_ids[k];
            if (edge_id >= 0 &&
                static_cast<std::size_t>(edge_id) < forbidden_edge_mask.size() &&
                forbidden_edge_mask[static_cast<std::size_t>(edge_id)] != 0U) {
                contains_forbidden_edge = true;
                break;
            }
        }

        if (!contains_forbidden_edge) {
            continue;
        }

        node_state.forbidden_entry[entry_index] = 1U;
        const std::size_t available_variant_count =
            node_state.entry_available_variant_count[entry_index];
        if (available_variant_count > 0U) {
            if (node_state.available_variant_count >= available_variant_count) {
                node_state.available_variant_count -= available_variant_count;
            } else {
                node_state.available_variant_count = 0U;
            }
            node_state.entry_available_variant_count[entry_index] = 0U;
            if (node_state.available_entry_count > 0U) {
                --node_state.available_entry_count;
            }
        }
    }
}

void mark_full_enumeration_pool_entries_forbidden_by_required_edges(
    const FullEnumerationStaticPool& static_pool,
    const std::vector<int>& required_outgoing_edge_by_node_id,
    const std::vector<int>& required_incoming_edge_by_node_id,
    FullEnumerationPoolNodeState& node_state
) {
    if (required_outgoing_edge_by_node_id.empty() &&
        required_incoming_edge_by_node_id.empty()) {
        return;
    }

    for (std::size_t entry_index = 0; entry_index < full_enumeration_entry_count(static_pool);
         ++entry_index) {
        if (full_enumeration_entry_is_forbidden(node_state, entry_index)) {
            continue;
        }

        bool violates_required_rule = false;
        for (std::size_t k = static_pool.entry_edge_ids_row_ptr[entry_index];
             k < static_pool.entry_edge_ids_row_ptr[entry_index + 1U];
             ++k) {
            const int edge_id = static_pool.entry_edge_ids[k];
            if (edge_id < 0 ||
                static_cast<std::size_t>(edge_id) >= static_pool.graph_edge_source_node_id.size() ||
                static_cast<std::size_t>(edge_id) >= static_pool.graph_edge_target_node_id.size()) {
                continue;
            }

            const NodeId source_node_id =
                static_pool.graph_edge_source_node_id[static_cast<std::size_t>(edge_id)];
            if (source_node_id >= 0 &&
                static_cast<std::size_t>(source_node_id) <
                    required_outgoing_edge_by_node_id.size()) {
                const int required_edge_id =
                    required_outgoing_edge_by_node_id[static_cast<std::size_t>(source_node_id)];
                if (required_edge_id >= 0 && required_edge_id != edge_id) {
                    violates_required_rule = true;
                    break;
                }
            }

            const NodeId target_node_id =
                static_pool.graph_edge_target_node_id[static_cast<std::size_t>(edge_id)];
            if (target_node_id >= 0 &&
                static_cast<std::size_t>(target_node_id) <
                    required_incoming_edge_by_node_id.size()) {
                const int required_edge_id =
                    required_incoming_edge_by_node_id[static_cast<std::size_t>(target_node_id)];
                if (required_edge_id >= 0 && required_edge_id != edge_id) {
                    violates_required_rule = true;
                    break;
                }
            }
        }

        if (!violates_required_rule) {
            continue;
        }

        node_state.forbidden_entry[entry_index] = 1U;
        const std::size_t available_variant_count =
            node_state.entry_available_variant_count[entry_index];
        if (available_variant_count > 0U) {
            if (node_state.available_variant_count >= available_variant_count) {
                node_state.available_variant_count -= available_variant_count;
            } else {
                node_state.available_variant_count = 0U;
            }
            node_state.entry_available_variant_count[entry_index] = 0U;
            if (node_state.available_entry_count > 0U) {
                --node_state.available_entry_count;
            }
        }
    }
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

const char* to_string(FullEnumerationRCUpdateMode mode) {
    switch (mode) {
        case FullEnumerationRCUpdateMode::Sequential:
            return "sequential";
        case FullEnumerationRCUpdateMode::Parallel:
            return "parallel";
    }
    return "unknown";
}

const char* to_string(FullEnumerationRCUpdateBackend backend) {
    switch (backend) {
        case FullEnumerationRCUpdateBackend::Custom:
            return "custom";
        case FullEnumerationRCUpdateBackend::OneMKL:
            return "onemkl";
    }
    return "unknown";
}

}  // namespace spdp
