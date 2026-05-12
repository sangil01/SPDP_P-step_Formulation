#include "PstepBnP.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <map>
#include <optional>
#include <ostream>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace spdp {
namespace {

constexpr double kTolerance = 1e-9;

bool double_less_or_equal(double lhs, double rhs) {
    return lhs <= rhs + kTolerance;
}

bool double_equal(double lhs, double rhs) {
    return std::fabs(lhs - rhs) <= kTolerance;
}

State canonicalize_state(State state) {
    if (state[1] < state[0]) {
        std::swap(state[0], state[1]);
    }
    return state;
}

std::string format_double(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(6) << value;
    return out.str();
}

bool is_lp_status_usable(int status) {
    return status == GRB_OPTIMAL || status == GRB_SUBOPTIMAL;
}

double remaining_time_seconds(
    const std::chrono::steady_clock::time_point start_time,
    double total_time_limit_seconds
) {
    if (total_time_limit_seconds < 0.0) {
        return total_time_limit_seconds;
    }
    const double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();
    return std::max(0.0, total_time_limit_seconds - elapsed);
}

enum class PiRequirement {
    Empty,
    NonEmpty,
};

using NodeStateTemplate = std::vector<std::optional<State>>;
using NodeStateTemplateList = std::vector<NodeStateTemplate>;

std::uint64_t node_pair_key(NodeId u, NodeId v) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(u)) << 32U) |
           static_cast<std::uint32_t>(v);
}

std::string build_column_key(
    const std::vector<int>& edge_ids,
    double tau
) {
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

bool matches_pi_requirement(
    const std::vector<int>& sequence_pi,
    PiRequirement requirement
) {
    switch (requirement) {
        case PiRequirement::Empty:
            return sequence_pi.empty();
        case PiRequirement::NonEmpty:
            return !sequence_pi.empty();
    }
    return false;
}

StateToken make_n_token() {
    return StateToken{'N', -1, -1};
}

StateToken make_e_token(int container_type) {
    return StateToken{'E', container_type, -1};
}

StateToken make_f_token(const Request& request) {
    return StateToken{'F', request.container_type, request.to_id};
}

StateToken empty_of_full_token(const StateToken& full_token) {
    if (full_token.kind != 'F') {
        throw std::runtime_error("Expected a full token when converting to empty token.");
    }
    return make_e_token(full_token.container_type);
}

State make_state(const StateToken& first, const StateToken& second) {
    State state{first, second};
    return canonicalize_state(state);
}

int request_index_for_node(const MultiDiGraph& graph, NodeId node_id) {
    const NodeSpec& node = graph.node(node_id);
    if (!node.request_idx.has_value()) {
        throw std::runtime_error("Expected request-indexed service node in heuristic-cg.");
    }
    return node.request_idx.value();
}

const Request& request_for_node(
    const SPDPData& data,
    const MultiDiGraph& graph,
    NodeId node_id
) {
    return data.requests.at(static_cast<std::size_t>(request_index_for_node(graph, node_id)));
}

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

bool violates_pickup_symmetry_43(
    const MultiDiGraph& graph,
    const std::vector<NodeId>& node_sequence
) {
    std::map<PricingRequestClassKey, int> max_request_idx_by_class;
    for (NodeId node_id : node_sequence) {
        const std::optional<PricingRequestClassKey> class_key =
            pickup_request_class_of_node(graph, node_id);
        if (!class_key.has_value()) {
            continue;
        }
        const int request_idx = request_index_for_node(graph, node_id);
        const auto found = max_request_idx_by_class.find(class_key.value());
        if (found != max_request_idx_by_class.end() && request_idx <= found->second) {
            return true;
        }
        max_request_idx_by_class[class_key.value()] = request_idx;
    }
    return false;
}

bool violates_delivery_symmetry_43(
    const MultiDiGraph& graph,
    const std::vector<NodeId>& node_sequence
) {
    std::map<PricingDeliveryClassKey, int> max_request_idx_by_class;
    for (NodeId node_id : node_sequence) {
        const std::optional<PricingDeliveryClassKey> class_key =
            delivery_request_class_of_node(graph, node_id);
        if (!class_key.has_value()) {
            continue;
        }
        const int request_idx = request_index_for_node(graph, node_id);
        const auto found = max_request_idx_by_class.find(class_key.value());
        if (found != max_request_idx_by_class.end() && request_idx <= found->second) {
            return true;
        }
        max_request_idx_by_class[class_key.value()] = request_idx;
    }
    return false;
}

std::vector<StateToken> collect_unique_full_tokens(const SPDPData& data) {
    std::set<StateToken> unique_tokens;
    for (const Request& request : data.requests) {
        unique_tokens.insert(make_f_token(request));
    }
    return {unique_tokens.begin(), unique_tokens.end()};
}

bool matches_state_at(
    const NodeStateTemplate& state_template,
    std::size_t position,
    const State& state
) {
    if (position >= state_template.size()) {
        return false;
    }
    if (!state_template[position].has_value()) {
        return true;
    }
    return canonicalize_state(state_template[position].value()) == canonicalize_state(state);
}

State fixed_pair_post_first_delivery_state(
    const Request& pickup_first_request,
    const Request& pickup_second_request,
    int pickup_first_request_index,
    int delivered_request_index
) {
    if (delivered_request_index == pickup_first_request_index) {
        return make_state(make_n_token(), make_e_token(pickup_second_request.container_type));
    }
    return make_state(make_e_token(pickup_first_request.container_type), make_n_token());
}

std::vector<double> heuristic_tau_candidates(
    const MultiDiGraph& graph,
    NodeId start_node_id,
    NodeId last_node_id,
    double total_time,
    double time_limit
) {
    const double latest_start = time_limit - total_time;
    if (latest_start < -kTolerance) {
        return {};
    }

    if (start_node_id == 0) {
        return {0.0};
    }
    if (last_node_id == graph.end_node_id()) {
        return {std::max(0.0, latest_start)};
    }

    const double latest = std::max(0.0, latest_start);
    if (double_equal(latest, 0.0)) {
        return {0.0};
    }
    return {0.0, latest};
}

CGColumn build_heuristic_column_from_graph_edges(
    const MultiDiGraph& graph,
    const std::vector<int>& edge_ids,
    double tau
) {
    if (edge_ids.empty()) {
        throw std::runtime_error("Heuristic column requires at least one graph edge.");
    }

    CGColumn column;
    column.q = static_cast<int>(edge_ids.size());
    column.tau = tau;
    column.reduced_cost = 0.0;

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

    return column;
}

void append_heuristic_column_variants(
    const MultiDiGraph& graph,
    double time_limit,
    const std::vector<int>& edge_ids,
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

    const std::vector<double> tau_candidates = heuristic_tau_candidates(
        graph,
        first_edge.u,
        last_edge.v,
        total_time,
        time_limit
    );

    for (double tau : tau_candidates) {
        const std::string key = build_column_key(edge_ids, tau);
        if (seen_columns.find(key) != seen_columns.end()) {
            continue;
        }
        columns.push_back(build_heuristic_column_from_graph_edges(graph, edge_ids, tau));
        seen_columns.emplace(key, columns.size() - 1U);
    }
}

struct FixedPatternOrder {
    int pickup_first = -1;
    int pickup_second = -1;
    int delivery_first = -1;
    int delivery_second = -1;
};

using PricingEdgePairMap = std::unordered_map<std::uint64_t, std::vector<int>>;

PricingEdgePairMap build_pricing_edge_pair_map(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context
) {
    PricingEdgePairMap result;
    for (std::size_t idx = 0; idx < context.edges.size(); ++idx) {
        const auto& pricing_edge = context.edges[idx];
        const EdgeRecord& edge =
            graph.edges()[static_cast<std::size_t>(pricing_edge.graph_edge_id)];
        result[node_pair_key(edge.u, edge.v)].push_back(static_cast<int>(idx));
    }
    return result;
}

template <typename Callback>
void enumerate_pricing_paths_for_node_sequence(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    const std::vector<NodeId>& node_sequence,
    const std::vector<PiRequirement>& pi_requirements,
    const NodeStateTemplateList& state_templates,
    Callback&& callback
) {
    if (node_sequence.size() < 2U ||
        pi_requirements.size() + 1U != node_sequence.size()) {
        throw std::runtime_error("Node sequence / pi requirement size mismatch.");
    }
    if (state_templates.empty()) {
        throw std::runtime_error("Node-state template list must not be empty.");
    }
    for (const NodeStateTemplate& state_template : state_templates) {
        if (state_template.size() != node_sequence.size()) {
            throw std::runtime_error("Node sequence / state template size mismatch.");
        }
    }

    std::vector<int> current_graph_edge_ids;
    current_graph_edge_ids.reserve(pi_requirements.size());
    std::vector<int> initial_templates;
    initial_templates.reserve(state_templates.size());
    for (std::size_t idx = 0; idx < state_templates.size(); ++idx) {
        initial_templates.push_back(static_cast<int>(idx));
    }

    const auto dfs = [&](const auto& self,
                         std::size_t leg_index,
                         int expected_from_node_state_index,
                         const std::vector<int>& active_templates)
        -> void {
        if (leg_index == pi_requirements.size()) {
            callback(current_graph_edge_ids);
            return;
        }

        const NodeId u = node_sequence[leg_index];
        const NodeId v = node_sequence[leg_index + 1U];
        const auto found = pricing_edges_by_pair.find(node_pair_key(u, v));
        if (found == pricing_edges_by_pair.end()) {
            return;
        }

        for (int pricing_edge_index : found->second) {
            const auto& pricing_edge =
                context.edges[static_cast<std::size_t>(pricing_edge_index)];
            if (expected_from_node_state_index >= 0 &&
                pricing_edge.from_node_state_index != expected_from_node_state_index) {
                continue;
            }

            const EdgeRecord& edge =
                graph.edges()[static_cast<std::size_t>(pricing_edge.graph_edge_id)];
            if (!matches_pi_requirement(edge.data.sequence_pi, pi_requirements[leg_index])) {
                continue;
            }
            const State start_state = canonicalize_state(edge.data.start_state);
            const State end_state = canonicalize_state(edge.data.end_state);
            std::vector<int> next_templates;
            next_templates.reserve(active_templates.size());
            for (int template_index : active_templates) {
                const NodeStateTemplate& state_template =
                    state_templates[static_cast<std::size_t>(template_index)];
                if (!matches_state_at(state_template, leg_index, start_state)) {
                    continue;
                }
                if (!matches_state_at(state_template, leg_index + 1U, end_state)) {
                    continue;
                }
                next_templates.push_back(template_index);
            }
            if (next_templates.empty()) {
                continue;
            }

            current_graph_edge_ids.push_back(pricing_edge.graph_edge_id);
            self(self, leg_index + 1U, pricing_edge.to_node_state_index, next_templates);
            current_graph_edge_ids.pop_back();
        }
    };

    dfs(dfs, 0U, -1, initial_templates);
}

std::vector<CGColumn> build_phase_one_heuristic_columns(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const NodeCGOptions& options
) {
    if (options.p != 3) {
        throw std::runtime_error("heuristic-cg currently requires p = 3.");
    }

    const int request_count = static_cast<int>(data.requests.size());
    const NodeId end_node_id = graph.end_node_id();
    const PricingEdgePairMap pricing_edges_by_pair = build_pricing_edge_pair_map(graph, context);
    const std::vector<StateToken> unique_full_tokens = collect_unique_full_tokens(data);

    std::unordered_map<std::string, std::size_t> seen_columns;
    std::vector<CGColumn> columns;

    const auto pickup_node_id = [](int request_idx) -> NodeId {
        return request_idx + 1;
    };
    const auto delivery_node_id = [request_count](int request_idx) -> NodeId {
        return request_count + request_idx + 1;
    };

    auto append_columns_for_template =
        [&](const std::vector<NodeId>& node_sequence,
            const std::vector<PiRequirement>& pi_reqs,
            const NodeStateTemplateList& state_templates) {
            if (options.prune_pickup_symmetry_43 &&
                violates_pickup_symmetry_43(graph, node_sequence)) {
                return;
            }
            if (options.prune_delivery_symmetry_43 &&
                violates_delivery_symmetry_43(graph, node_sequence)) {
                return;
            }
            enumerate_pricing_paths_for_node_sequence(
                graph,
                context,
                pricing_edges_by_pair,
                node_sequence,
                pi_reqs,
                state_templates,
                [&](const std::vector<int>& edge_ids) {
                    append_heuristic_column_variants(
                        graph,
                        options.time_limit,
                        edge_ids,
                        seen_columns,
                        columns
                    );
                }
            );
        };

    // 1. depot boundary columns
    for (int request_idx = 0; request_idx < request_count; ++request_idx) {
        const Request& request = data.requests.at(static_cast<std::size_t>(request_idx));
        append_columns_for_template(
            {0, pickup_node_id(request_idx)},
            {PiRequirement::Empty},
            {{std::nullopt, make_state(make_n_token(), make_f_token(request))}}
        );
        append_columns_for_template(
            {delivery_node_id(request_idx), end_node_id},
            {PiRequirement::Empty},
            {{make_state(make_n_token(), make_n_token()), std::nullopt}}
        );
    }

    // 2-5. fixed-return families
    for (int i = 0; i < request_count; ++i) {
        for (int j = i + 1; j < request_count; ++j) {
            const Request& request_i = data.requests.at(static_cast<std::size_t>(i));
            const Request& request_j = data.requests.at(static_cast<std::size_t>(j));
            const NodeId pickup_i = pickup_node_id(i);
            const NodeId pickup_j = pickup_node_id(j);
            const NodeId delivery_i = delivery_node_id(i);
            const NodeId delivery_j = delivery_node_id(j);
            const std::vector<FixedPatternOrder> patterns = {
                {pickup_i, pickup_j, delivery_i, delivery_j},
                {pickup_i, pickup_j, delivery_j, delivery_i},
                {pickup_j, pickup_i, delivery_i, delivery_j},
                {pickup_j, pickup_i, delivery_j, delivery_i},
            };

            for (const FixedPatternOrder& pattern : patterns) {
                const Request& pickup_first_request =
                    request_for_node(data, graph, pattern.pickup_first);
                const Request& pickup_second_request =
                    request_for_node(data, graph, pattern.pickup_second);
                const int pickup_first_request_index =
                    request_index_for_node(graph, pattern.pickup_first);
                const int delivery_first_request_index =
                    request_index_for_node(graph, pattern.delivery_first);
                const State first_delivery_state = fixed_pair_post_first_delivery_state(
                    pickup_first_request,
                    pickup_second_request,
                    pickup_first_request_index,
                    delivery_first_request_index
                );
                append_columns_for_template(
                    {
                        pattern.pickup_first,
                        pattern.pickup_second,
                        pattern.delivery_first,
                        pattern.delivery_second,
                    },
                    {PiRequirement::Empty, PiRequirement::NonEmpty, PiRequirement::Empty},
                    {{
                        make_state(make_n_token(), make_f_token(pickup_first_request)),
                        make_state(
                            make_f_token(pickup_first_request),
                            make_f_token(pickup_second_request)
                        ),
                        first_delivery_state,
                        make_state(make_n_token(), make_n_token()),
                    }}
                );

                for (int k = 0; k < request_count; ++k) {
                    if (k == i || k == j) {
                        continue;
                    }
                    const NodeId pickup_k = pickup_node_id(k);
                    const NodeId delivery_k = delivery_node_id(k);
                    const Request& request_k = data.requests.at(static_cast<std::size_t>(k));

                    append_columns_for_template(
                        {
                            pattern.pickup_second,
                            pattern.delivery_first,
                            pattern.delivery_second,
                            pickup_k,
                        },
                        {PiRequirement::NonEmpty, PiRequirement::Empty, PiRequirement::Empty},
                        {{
                            make_state(
                                make_f_token(pickup_first_request),
                                make_f_token(pickup_second_request)
                            ),
                            first_delivery_state,
                            make_state(make_n_token(), make_n_token()),
                            make_state(make_n_token(), make_f_token(request_k)),
                        }}
                    );
                }

            }

            for (int k = 0; k < request_count; ++k) {
                if (k == i || k == j) {
                    continue;
                }
                const NodeId pickup_k = pickup_node_id(k);
                const NodeId delivery_k = delivery_node_id(k);
                const Request& request_k = data.requests.at(static_cast<std::size_t>(k));
                const std::vector<std::tuple<NodeId, NodeId, State>> delivery_orders = {
                    {
                        delivery_i,
                        delivery_j,
                        make_state(make_n_token(), make_e_token(request_j.container_type)),
                    },
                    {
                        delivery_j,
                        delivery_i,
                        make_state(make_n_token(), make_e_token(request_i.container_type)),
                    },
                };

                for (int m = 0; m < request_count; ++m) {
                    if (m == i || m == j || m == k) {
                        continue;
                    }
                    const NodeId pickup_m = pickup_node_id(m);
                    const NodeId delivery_m = delivery_node_id(m);
                    const Request& request_m = data.requests.at(static_cast<std::size_t>(m));

                    for (const auto& delivery_order : delivery_orders) {
                        const NodeId delivery_first = std::get<0>(delivery_order);
                        const NodeId delivery_second = std::get<1>(delivery_order);
                        const State first_delivery_state = std::get<2>(delivery_order);

                        append_columns_for_template(
                            {delivery_first, delivery_second, pickup_k, pickup_m},
                            {PiRequirement::Empty, PiRequirement::Empty, PiRequirement::Empty},
                            {{
                                first_delivery_state,
                                make_state(make_n_token(), make_n_token()),
                                make_state(make_n_token(), make_f_token(request_k)),
                                make_state(make_f_token(request_k), make_f_token(request_m)),
                            }}
                        );

                        append_columns_for_template(
                            {delivery_second, pickup_k, pickup_m, delivery_k},
                            {PiRequirement::Empty, PiRequirement::Empty, PiRequirement::NonEmpty},
                            {{
                                make_state(make_n_token(), make_n_token()),
                                make_state(make_n_token(), make_f_token(request_k)),
                                make_state(make_f_token(request_k), make_f_token(request_m)),
                                make_state(
                                    make_n_token(),
                                    make_e_token(request_m.container_type)
                                ),
                            }}
                        );

                        append_columns_for_template(
                            {delivery_second, pickup_k, pickup_m, delivery_m},
                            {PiRequirement::Empty, PiRequirement::Empty, PiRequirement::NonEmpty},
                            {{
                                make_state(make_n_token(), make_n_token()),
                                make_state(make_n_token(), make_f_token(request_k)),
                                make_state(make_f_token(request_k), make_f_token(request_m)),
                                make_state(
                                    make_e_token(request_k.container_type),
                                    make_n_token()
                                ),
                            }}
                        );
                    }
                }
            }
        }
    }

    // 6-8. flexible-return families
    for (int i = 0; i < request_count; ++i) {
        const NodeId pickup_i = pickup_node_id(i);
        const NodeId delivery_i = delivery_node_id(i);
        const Request& request_i = data.requests.at(static_cast<std::size_t>(i));

        for (int j = 0; j < request_count; ++j) {
            if (j == i) {
                continue;
            }
            const NodeId pickup_j = pickup_node_id(j);
            const NodeId delivery_j = delivery_node_id(j);
            const Request& request_j = data.requests.at(static_cast<std::size_t>(j));

            append_columns_for_template(
                {delivery_i, pickup_i, delivery_j, pickup_j},
                {PiRequirement::Empty, PiRequirement::Empty, PiRequirement::Empty},
                {{
                    make_state(make_n_token(), make_e_token(request_j.container_type)),
                    make_state(make_f_token(request_i), make_e_token(request_j.container_type)),
                    make_state(make_f_token(request_i), make_n_token()),
                    make_state(make_f_token(request_i), make_f_token(request_j)),
                }}
            );

            NodeStateTemplateList middle_templates_3;
            for (const StateToken& carried_full : unique_full_tokens) {
                middle_templates_3.push_back({
                    make_state(carried_full, make_n_token()),
                    make_state(carried_full, make_f_token(request_i)),
                    make_state(make_n_token(), make_e_token(request_i.container_type)),
                    make_state(make_f_token(request_j), make_e_token(request_i.container_type)),
                });
                middle_templates_3.push_back({
                    make_state(carried_full, make_n_token()),
                    make_state(carried_full, make_f_token(request_i)),
                    make_state(empty_of_full_token(carried_full), make_n_token()),
                    make_state(empty_of_full_token(carried_full), make_f_token(request_j)),
                });
            }
            append_columns_for_template(
                {delivery_i, pickup_i, delivery_j, pickup_j},
                {PiRequirement::Empty, PiRequirement::NonEmpty, PiRequirement::Empty},
                middle_templates_3
            );

            for (int k = 0; k < request_count; ++k) {
                if (k == i || k == j) {
                    continue;
                }
                const NodeId pickup_k = pickup_node_id(k);
                const NodeId delivery_k = delivery_node_id(k);
                const Request& request_k = data.requests.at(static_cast<std::size_t>(k));
                NodeStateTemplateList entry_templates;
                entry_templates.push_back({
                    make_state(make_n_token(), make_f_token(request_i)),
                    make_state(make_f_token(request_i), make_f_token(request_j)),
                    make_state(make_n_token(), make_e_token(request_j.container_type)),
                });
                entry_templates.push_back({
                    make_state(make_n_token(), make_f_token(request_i)),
                    make_state(make_f_token(request_i), make_f_token(request_j)),
                    make_state(make_e_token(request_i.container_type), make_n_token()),
                });

                append_columns_for_template(
                    {pickup_i, pickup_j, delivery_k},
                    {PiRequirement::Empty, PiRequirement::NonEmpty},
                    entry_templates
                );

                NodeStateTemplateList middle_templates_2;
                middle_templates_2.push_back({
                    make_state(make_f_token(request_i), make_e_token(request_j.container_type)),
                    make_state(make_f_token(request_i), make_n_token()),
                    make_state(make_f_token(request_i), make_f_token(request_j)),
                    make_state(make_n_token(), make_e_token(request_j.container_type)),
                });
                middle_templates_2.push_back({
                    make_state(make_f_token(request_i), make_e_token(request_j.container_type)),
                    make_state(make_f_token(request_i), make_n_token()),
                    make_state(make_f_token(request_i), make_f_token(request_j)),
                    make_state(make_e_token(request_i.container_type), make_n_token()),
                });
                append_columns_for_template(
                    {pickup_i, delivery_j, pickup_j, delivery_k},
                    {PiRequirement::Empty, PiRequirement::Empty, PiRequirement::NonEmpty},
                    middle_templates_2
                );
                NodeStateTemplateList middle_templates_4;
                for (const StateToken& carried_full : unique_full_tokens) {
                    middle_templates_4.push_back({
                        make_state(carried_full, make_f_token(request_i)),
                        make_state(make_n_token(), make_e_token(request_i.container_type)),
                        make_state(
                            make_f_token(request_j),
                            make_e_token(request_i.container_type)
                        ),
                        make_state(make_f_token(request_j), make_n_token()),
                    });
                    middle_templates_4.push_back({
                        make_state(carried_full, make_f_token(request_i)),
                        make_state(empty_of_full_token(carried_full), make_n_token()),
                        make_state(
                            empty_of_full_token(carried_full),
                            make_f_token(request_j)
                        ),
                        make_state(make_n_token(), make_f_token(request_j)),
                    });
                }
                append_columns_for_template(
                    {pickup_i, delivery_j, pickup_j, delivery_k},
                    {PiRequirement::NonEmpty, PiRequirement::Empty, PiRequirement::Empty},
                    middle_templates_4
                );

                NodeStateTemplateList exit_templates;
                for (const StateToken& carried_full : unique_full_tokens) {
                    exit_templates.push_back({
                        make_state(carried_full, make_f_token(request_i)),
                        make_state(make_n_token(), make_e_token(request_i.container_type)),
                        make_state(make_n_token(), make_n_token()),
                    });
                    exit_templates.push_back({
                        make_state(carried_full, make_f_token(request_i)),
                        make_state(empty_of_full_token(carried_full), make_n_token()),
                        make_state(make_n_token(), make_n_token()),
                    });
                }
                append_columns_for_template(
                    {pickup_i, delivery_j, delivery_k},
                    {PiRequirement::NonEmpty, PiRequirement::Empty},
                    exit_templates
                );
            }
        }
    }

    return columns;
}

struct PhaseTwoColumnPoolEntry {
    CGColumn column;
    double best_reduced_cost = std::numeric_limits<double>::infinity();
};

std::vector<int> build_ladder_level_starts(
    const std::vector<int>& ordered_starts,
    std::size_t max_start_cap,
    double start_ratio,
    std::size_t ladder_levels,
    std::size_t level_index,
    std::size_t& explored_prefix_size
) {
    if (explored_prefix_size >= ordered_starts.size()) {
        return {};
    }

    std::size_t current_prefix_size = ordered_starts.size();
    if (ladder_levels > 0U) {
        const std::size_t ratio_limit = static_cast<std::size_t>(
            std::ceil(start_ratio * static_cast<double>(ordered_starts.size()))
        );
        const std::size_t final_prefix_size = std::min(
            max_start_cap,
            std::min(ratio_limit, ordered_starts.size())
        );
        const std::size_t effective_final_prefix_size =
            std::max<std::size_t>(1U, final_prefix_size);
        current_prefix_size = static_cast<std::size_t>(
            std::ceil(
                static_cast<double>(level_index) *
                static_cast<double>(effective_final_prefix_size) /
                static_cast<double>(ladder_levels)
            )
        );
        current_prefix_size = std::min(current_prefix_size, effective_final_prefix_size);
    }

    current_prefix_size = std::min(current_prefix_size, ordered_starts.size());
    if (current_prefix_size <= explored_prefix_size) {
        return {};
    }

    std::vector<int> level_starts;
    level_starts.reserve(current_prefix_size - explored_prefix_size);
    for (std::size_t idx = explored_prefix_size; idx < current_prefix_size; ++idx) {
        level_starts.push_back(ordered_starts[idx]);
    }
    explored_prefix_size = current_prefix_size;
    return level_starts;
}

void rebuild_column_pool_index(
    const std::vector<PhaseTwoColumnPoolEntry>& column_pool,
    std::unordered_map<std::string, std::size_t>& pool_index_by_key
) {
    pool_index_by_key.clear();
    for (std::size_t idx = 0; idx < column_pool.size(); ++idx) {
        pool_index_by_key.emplace(
            build_column_key(column_pool[idx].column.edge_ids, column_pool[idx].column.tau),
            idx
        );
    }
}

void prune_column_pool(
    std::vector<PhaseTwoColumnPoolEntry>& column_pool,
    std::unordered_map<std::string, std::size_t>& pool_index_by_key,
    const CGMasterProblem& master_problem,
    std::size_t max_pool_size
) {
    column_pool.erase(
        std::remove_if(
            column_pool.begin(),
            column_pool.end(),
            [&](const PhaseTwoColumnPoolEntry& entry) {
                const std::string key =
                    build_column_key(entry.column.edge_ids, entry.column.tau);
                return master_problem.column_id_by_key.find(key) !=
                       master_problem.column_id_by_key.end();
            }
        ),
        column_pool.end()
    );

    if (column_pool.size() > max_pool_size) {
        std::sort(
            column_pool.begin(),
            column_pool.end(),
            [](const PhaseTwoColumnPoolEntry& lhs, const PhaseTwoColumnPoolEntry& rhs) {
                if (!double_equal(lhs.best_reduced_cost, rhs.best_reduced_cost)) {
                    return lhs.best_reduced_cost < rhs.best_reduced_cost;
                }
                return lhs.column.reduced_cost < rhs.column.reduced_cost;
            }
        );
        column_pool.resize(max_pool_size);
    }

    rebuild_column_pool_index(column_pool, pool_index_by_key);
}

void add_columns_to_phase_two_pool(
    const std::vector<CGColumn>& columns,
    const CGMasterProblem& master_problem,
    std::size_t max_pool_size,
    std::vector<PhaseTwoColumnPoolEntry>& column_pool,
    std::unordered_map<std::string, std::size_t>& pool_index_by_key
) {
    for (const CGColumn& column : columns) {
        const std::string key = build_column_key(column.edge_ids, column.tau);
        if (master_problem.column_id_by_key.find(key) != master_problem.column_id_by_key.end()) {
            continue;
        }
        const auto found = pool_index_by_key.find(key);
        if (found != pool_index_by_key.end()) {
            auto& entry = column_pool[found->second];
            entry.best_reduced_cost = std::min(entry.best_reduced_cost, column.reduced_cost);
            entry.column = column;
            continue;
        }

        PhaseTwoColumnPoolEntry entry;
        entry.column = column;
        entry.best_reduced_cost = column.reduced_cost;
        column_pool.push_back(std::move(entry));
    }

    prune_column_pool(column_pool, pool_index_by_key, master_problem, max_pool_size);
}

ForwardPricingResult reprice_phase_two_column_pool(
    const CGMasterProblem& master_problem,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    double reduced_cost_tolerance,
    std::size_t max_total_columns,
    std::size_t max_reprice,
    std::vector<PhaseTwoColumnPoolEntry>& column_pool,
    std::unordered_map<std::string, std::size_t>& pool_index_by_key
) {
    ForwardPricingResult result;
    result.best_reduced_cost = std::numeric_limits<double>::infinity();

    prune_column_pool(
        column_pool,
        pool_index_by_key,
        master_problem,
        column_pool.size()
    );

    std::vector<std::size_t> candidate_indices(column_pool.size(), 0U);
    for (std::size_t idx = 0; idx < column_pool.size(); ++idx) {
        candidate_indices[idx] = idx;
    }
    std::sort(
        candidate_indices.begin(),
        candidate_indices.end(),
        [&](std::size_t lhs, std::size_t rhs) {
            if (!double_equal(
                    column_pool[lhs].best_reduced_cost,
                    column_pool[rhs].best_reduced_cost
                )) {
                return column_pool[lhs].best_reduced_cost <
                       column_pool[rhs].best_reduced_cost;
            }
            return column_pool[lhs].column.reduced_cost <
                   column_pool[rhs].column.reduced_cost;
        }
    );

    const std::size_t reprice_count = std::min(max_reprice, candidate_indices.size());
    std::vector<CGColumn> negative_columns;
    negative_columns.reserve(reprice_count);
    for (std::size_t pos = 0; pos < reprice_count; ++pos) {
        auto& entry = column_pool[candidate_indices[pos]];
        entry.column.reduced_cost = evaluate_column_reduced_cost(
            dual_solution,
            phase,
            entry.column
        );
        entry.best_reduced_cost = std::min(entry.best_reduced_cost, entry.column.reduced_cost);
        result.best_reduced_cost = std::min(result.best_reduced_cost, entry.column.reduced_cost);
        if (entry.column.reduced_cost < reduced_cost_tolerance) {
            negative_columns.push_back(entry.column);
        }
    }

    std::sort(
        negative_columns.begin(),
        negative_columns.end(),
        [](const CGColumn& lhs, const CGColumn& rhs) {
            return lhs.reduced_cost < rhs.reduced_cost;
        }
    );
    for (const CGColumn& column : negative_columns) {
        if (result.columns.size() >= max_total_columns) {
            break;
        }
        result.columns.push_back(column);
    }

    if (!std::isfinite(result.best_reduced_cost)) {
        result.best_reduced_cost = 0.0;
    }
    result.status = result.columns.empty()
                        ? ForwardPricingStatus::NoNegativeColumn
                        : ForwardPricingStatus::ColumnsFound;
    return result;
}

void write_iteration_log(
    std::ostream& out,
    const NodeCGIterationLog& log
) {
    out << "[node-cg] phase="
        << (log.phase == CGPhase::PhaseI ? "I" : "II")
        << " iter=" << log.iteration_index
        << " lp_obj=" << format_double(log.lp_objective_value)
        << " artificial_sum=" << format_double(log.artificial_sum)
        << " best_rc=" << format_double(log.best_reduced_cost)
        << " pricing_status=" << to_string(log.pricing_status)
        << " pricing_seconds=" << format_double(log.pricing_runtime_seconds)
        << " added_columns=" << log.added_column_count
        << " complete_labels=" << log.complete_label_count
        << " start_labels=" << log.start_label_count
        << " generated_labels=" << log.generated_label_count
        << " surviving_labels=" << log.surviving_label_count
        << " dominated_labels=" << log.dominated_label_count
        << " heuristic_attempted=" << (log.heuristic_pricing_attempted ? 1 : 0)
        << " exact_fallback=" << (log.exact_pricing_fallback_used ? 1 : 0)
        << " heuristic_pricing_starts=" << log.heuristic_explored_start_count
        << " heuristic_pricing_cols=" << log.heuristic_found_column_count
        << '\n';
}

void capture_column_values_if_available(
    const CGMasterProblem& master_problem,
    NodeCGResult& result
) {
    result.columns = master_problem.columns;
    result.column_values.clear();

    if (master_problem.model->get(GRB_IntAttr_SolCount) <= 0) {
        return;
    }

    result.column_values.resize(master_problem.x_vars.size(), 0.0);
    for (std::size_t idx = 0; idx < master_problem.x_vars.size(); ++idx) {
        result.column_values[idx] = master_problem.x_vars[idx].get(GRB_DoubleAttr_X);
    }
}

bool run_node_pricing_phase(
    const MultiDiGraph& graph,
    const ForwardPricingContext& pricing_context,
    CGMasterProblem& master_problem,
    const NodeCGOptions& options,
    CGPhase phase,
    const std::chrono::steady_clock::time_point global_start_time,
    NodeCGResult& result,
    std::ostream* log_stream
) {
    ForwardPricingOptions exact_pricing_options;
    exact_pricing_options.max_columns_per_start =
        options.exact_pricing_max_columns_per_start;
    exact_pricing_options.max_total_columns =
        options.exact_pricing_max_total_columns_per_round;
    exact_pricing_options.reduced_cost_tolerance = options.reduced_cost_tolerance;
    exact_pricing_options.prune_pickup_symmetry_43 = options.prune_pickup_symmetry_43;
    exact_pricing_options.prune_delivery_symmetry_43 = options.prune_delivery_symmetry_43;

    ForwardPricingOptions heuristic_pricing_options;
    heuristic_pricing_options.max_columns_per_start =
        options.phase_two_heuristic_max_columns_per_start;
    heuristic_pricing_options.max_total_columns =
        options.phase_two_heuristic_max_total_columns;
    heuristic_pricing_options.reduced_cost_tolerance = options.reduced_cost_tolerance;
    heuristic_pricing_options.prune_pickup_symmetry_43 = options.prune_pickup_symmetry_43;
    heuristic_pricing_options.prune_delivery_symmetry_43 = options.prune_delivery_symmetry_43;
    heuristic_pricing_options.heuristic_pricing = true;
    heuristic_pricing_options.heuristic_max_starts = options.phase_two_heuristic_max_starts;
    heuristic_pricing_options.heuristic_start_ratio = options.phase_two_heuristic_start_ratio;
    heuristic_pricing_options.heuristic_search_column_ratio =
        options.phase_two_heuristic_search_column_ratio;
    heuristic_pricing_options.heuristic_start_score_mode =
        options.phase_two_heuristic_start_score_mode;
    heuristic_pricing_options.labeling_top_k_next =
        options.phase_two_labeling_top_k_next;
    heuristic_pricing_options.shallow_k1 = options.phase_two_shallow_k1;
    heuristic_pricing_options.shallow_k2 = options.phase_two_shallow_k2;

    CGLPSnapshot last_solved_snapshot;
    bool has_last_solved_snapshot = false;
    std::vector<PhaseTwoColumnPoolEntry> phase_two_column_pool;
    std::unordered_map<std::string, std::size_t> phase_two_pool_index_by_key;

    for (std::size_t iteration = 1; iteration <= options.max_iterations_per_phase; ++iteration) {
        const double remaining = remaining_time_seconds(global_start_time, options.solver_time_limit);
        if (remaining >= 0.0) {
            master_problem.model->set(GRB_DoubleParam_TimeLimit, remaining);
        }

        solve_cg_master_lp(master_problem);
        const CGLPSnapshot snapshot = capture_cg_master_snapshot(master_problem);
        last_solved_snapshot = snapshot;
        has_last_solved_snapshot = true;

        NodeCGIterationLog iteration_log;
        iteration_log.phase = phase;
        iteration_log.iteration_index = iteration;
        iteration_log.lp_objective_value = snapshot.objective_value;
        iteration_log.artificial_sum = snapshot.artificial_sum;

        if (log_stream != nullptr) {
            write_cg_master_snapshot(*log_stream, snapshot);
        }

        if (!is_lp_status_usable(snapshot.gurobi_status)) {
            result.final_snapshot = snapshot;
            result.iteration_logs.push_back(iteration_log);
            if (phase == CGPhase::PhaseI) {
                result.phase_one_objective = snapshot.objective_value;
                result.phase_one_iterations = iteration;
            } else {
                result.phase_two_objective = snapshot.objective_value;
                result.phase_two_iterations = iteration;
            }
            if (log_stream != nullptr) {
                *log_stream << "[node-cg] stopping because LP status is not usable for pricing.\n";
            }
            return false;
        }

        if (phase == CGPhase::PhaseI && snapshot.artificial_sum <= 1e-6) {
            result.phase_one_feasible = true;
            result.phase_one_objective = snapshot.objective_value;
            result.phase_one_iterations = iteration;
            result.final_snapshot = snapshot;
            result.iteration_logs.push_back(iteration_log);
            if (log_stream != nullptr) {
                *log_stream << "[node-cg] Phase I reached zero artificial sum.\n";
            }
            return true;
        }

        const CGDualSolution dual_solution = extract_cg_master_duals(master_problem);
        const auto pricing_start_time = std::chrono::steady_clock::now();
        ForwardPricingResult pricing_result;
        std::vector<CGColumn> columns_to_add_to_pool;
        if (phase == CGPhase::PhaseII &&
            options.phase_two_pricing_mode ==
                NodeCGPhaseTwoPricingMode::HeuristicPricingThenExact) {
            iteration_log.heuristic_pricing_attempted = true;
            ForwardPricingResult heuristic_stats;
            heuristic_stats.best_reduced_cost = std::numeric_limits<double>::infinity();
            ForwardPricingResult selected_heuristic_result;
            bool heuristic_found = false;

            if (options.phase_two_column_pool_enabled) {
                const ForwardPricingResult pool_result = reprice_phase_two_column_pool(
                    master_problem,
                    dual_solution,
                    phase,
                    options.reduced_cost_tolerance,
                    options.phase_two_heuristic_max_total_columns,
                    options.phase_two_column_pool_max_reprice,
                    phase_two_column_pool,
                    phase_two_pool_index_by_key
                );
                heuristic_stats.best_reduced_cost =
                    std::min(heuristic_stats.best_reduced_cost, pool_result.best_reduced_cost);
                if (!pool_result.columns.empty()) {
                    selected_heuristic_result = pool_result;
                    heuristic_found = true;
                }
            }

            if (!heuristic_found) {
                const std::vector<int> ordered_starts = build_heuristic_start_order(
                    graph,
                    pricing_context,
                    dual_solution,
                    phase,
                    options.phase_two_heuristic_start_score_mode
                );
                std::size_t explored_prefix_size = 0U;
                const std::size_t ladder_levels =
                    (options.phase_two_heuristic_ladder_levels == 0U)
                        ? 1U
                        : options.phase_two_heuristic_ladder_levels;

                for (std::size_t level = 1; level <= ladder_levels; ++level) {
                    const std::vector<int> level_starts = build_ladder_level_starts(
                        ordered_starts,
                        options.phase_two_heuristic_max_starts,
                        options.phase_two_heuristic_start_ratio,
                        options.phase_two_heuristic_ladder_levels,
                        level,
                        explored_prefix_size
                    );
                    if (level_starts.empty()) {
                        continue;
                    }

                    ForwardPricingOptions level_options = heuristic_pricing_options;
                    level_options.explicit_start_node_state_indices = level_starts;

                    ForwardPricingResult level_result;
                    if (options.phase_two_heuristic_engine ==
                        NodeCGPhaseTwoHeuristicEngine::ShallowSearch) {
                        level_result = run_phase_two_shallow_search(
                            graph,
                            pricing_context,
                            dual_solution,
                            phase,
                            level_options
                        );
                    } else {
                        level_result = run_forward_pricing(
                            graph,
                            pricing_context,
                            dual_solution,
                            phase,
                            level_options
                        );
                    }

                    heuristic_stats.best_reduced_cost =
                        std::min(heuristic_stats.best_reduced_cost, level_result.best_reduced_cost);
                    heuristic_stats.complete_label_count += level_result.complete_label_count;
                    heuristic_stats.start_label_count += level_result.start_label_count;
                    heuristic_stats.generated_label_count += level_result.generated_label_count;
                    heuristic_stats.surviving_label_count += level_result.surviving_label_count;
                    heuristic_stats.dominated_label_count += level_result.dominated_label_count;

                    if (!level_result.columns.empty()) {
                        selected_heuristic_result = level_result;
                        heuristic_found = true;
                        columns_to_add_to_pool = level_result.deferred_columns;
                        break;
                    }
                }
            }

            iteration_log.heuristic_explored_start_count = heuristic_stats.start_label_count;
            iteration_log.heuristic_found_column_count =
                heuristic_found ? selected_heuristic_result.columns.size() : 0U;

            if (heuristic_found) {
                pricing_result = selected_heuristic_result;
                pricing_result.best_reduced_cost = heuristic_stats.best_reduced_cost;
                pricing_result.complete_label_count = heuristic_stats.complete_label_count;
                pricing_result.start_label_count = heuristic_stats.start_label_count;
                pricing_result.generated_label_count = heuristic_stats.generated_label_count;
                pricing_result.surviving_label_count = heuristic_stats.surviving_label_count;
                pricing_result.dominated_label_count = heuristic_stats.dominated_label_count;
            } else {
                iteration_log.exact_pricing_fallback_used = true;
                const ForwardPricingResult exact_result = run_forward_pricing(
                    graph,
                    pricing_context,
                    dual_solution,
                    phase,
                    exact_pricing_options
                );
                pricing_result = exact_result;
                pricing_result.best_reduced_cost =
                    std::min(heuristic_stats.best_reduced_cost, exact_result.best_reduced_cost);
                pricing_result.complete_label_count += heuristic_stats.complete_label_count;
                pricing_result.start_label_count += heuristic_stats.start_label_count;
                pricing_result.generated_label_count += heuristic_stats.generated_label_count;
                pricing_result.surviving_label_count += heuristic_stats.surviving_label_count;
                pricing_result.dominated_label_count += heuristic_stats.dominated_label_count;
                columns_to_add_to_pool = exact_result.deferred_columns;
            }
        } else {
            pricing_result = run_forward_pricing(
                graph,
                pricing_context,
                dual_solution,
                phase,
                exact_pricing_options
            );
        }
        const auto pricing_end_time = std::chrono::steady_clock::now();
        iteration_log.best_reduced_cost = pricing_result.best_reduced_cost;
        iteration_log.pricing_status = pricing_result.status;
        iteration_log.pricing_runtime_seconds =
            std::chrono::duration<double>(pricing_end_time - pricing_start_time).count();
        iteration_log.complete_label_count = pricing_result.complete_label_count;
        iteration_log.start_label_count = pricing_result.start_label_count;
        iteration_log.generated_label_count = pricing_result.generated_label_count;
        iteration_log.surviving_label_count = pricing_result.surviving_label_count;
        iteration_log.dominated_label_count = pricing_result.dominated_label_count;

        const std::size_t added_columns =
            add_columns_to_cg_master(master_problem, pricing_result.columns);
        if (phase == CGPhase::PhaseII &&
            options.phase_two_column_pool_enabled &&
            !columns_to_add_to_pool.empty()) {
            add_columns_to_phase_two_pool(
                columns_to_add_to_pool,
                master_problem,
                options.phase_two_column_pool_max_size,
                phase_two_column_pool,
                phase_two_pool_index_by_key
            );
        }
        iteration_log.added_column_count = added_columns;
        result.iteration_logs.push_back(iteration_log);
        result.total_columns_added += added_columns;

        if (log_stream != nullptr) {
            write_iteration_log(*log_stream, iteration_log);
        }

        if (added_columns == 0U) {
            result.final_snapshot = snapshot;
            if (phase == CGPhase::PhaseI) {
                result.phase_one_feasible = false;
                result.phase_one_objective = snapshot.objective_value;
                result.phase_one_iterations = iteration;
            } else {
                result.phase_two_objective = snapshot.objective_value;
                result.phase_two_iterations = iteration;
                result.solved_to_completion = true;
            }
            return phase == CGPhase::PhaseII;
        }
    }

    result.final_snapshot = has_last_solved_snapshot
                                ? last_solved_snapshot
                                : capture_cg_master_snapshot(master_problem);
    if (phase == CGPhase::PhaseI) {
        result.phase_one_objective = result.final_snapshot.objective_value;
        result.phase_one_iterations = options.max_iterations_per_phase;
    } else {
        result.phase_two_objective = result.final_snapshot.objective_value;
        result.phase_two_iterations = options.max_iterations_per_phase;
    }
    return false;
}

bool run_node_phase_one_exact_cg(
    const MultiDiGraph& graph,
    const ForwardPricingContext& pricing_context,
    CGMasterProblem& master_problem,
    const NodeCGOptions& options,
    const std::chrono::steady_clock::time_point global_start_time,
    NodeCGResult& result,
    std::ostream* log_stream
) {
    if (log_stream != nullptr) {
        *log_stream << "[node-cg] Phase I mode=exact-cg\n";
    }
    return run_node_pricing_phase(
        graph,
        pricing_context,
        master_problem,
        options,
        CGPhase::PhaseI,
        global_start_time,
        result,
        log_stream
    );
}

bool run_node_phase_one_heuristic_cg(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& pricing_context,
    CGMasterProblem& master_problem,
    const NodeCGOptions& options,
    const std::chrono::steady_clock::time_point global_start_time,
    NodeCGResult& result,
    std::ostream* log_stream
) {
    if (log_stream != nullptr) {
        *log_stream << "[node-cg] Phase I mode=heuristic-cg\n";
    }

    const auto heuristic_start_time = std::chrono::steady_clock::now();
    const std::vector<CGColumn> heuristic_columns = build_phase_one_heuristic_columns(
        data,
        graph,
        pricing_context,
        options
    );
    const auto heuristic_end_time = std::chrono::steady_clock::now();
    const std::size_t added_seed_columns =
        add_columns_to_cg_master(master_problem, heuristic_columns);
    result.total_columns_added += added_seed_columns;

    if (log_stream != nullptr) {
        *log_stream << "[node-cg] heuristic seeds generated=" << heuristic_columns.size()
                    << " inserted=" << added_seed_columns
                    << " runtime_seconds="
                    << format_double(std::chrono::duration<double>(
                                         heuristic_end_time - heuristic_start_time)
                                         .count())
                    << '\n';
    }

    return run_node_pricing_phase(
        graph,
        pricing_context,
        master_problem,
        options,
        CGPhase::PhaseI,
        global_start_time,
        result,
        log_stream
    );
}

bool run_node_phase_one(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& pricing_context,
    CGMasterProblem& master_problem,
    const NodeCGOptions& options,
    const std::chrono::steady_clock::time_point global_start_time,
    NodeCGResult& result,
    std::ostream* log_stream
) {
    switch (options.phase_one_mode) {
        case NodeCGPhaseOneMode::ExactCG:
            return run_node_phase_one_exact_cg(
                graph,
                pricing_context,
                master_problem,
                options,
                global_start_time,
                result,
                log_stream
            );
        case NodeCGPhaseOneMode::HeuristicCG:
            return run_node_phase_one_heuristic_cg(
                data,
                graph,
                pricing_context,
                master_problem,
                options,
                global_start_time,
                result,
                log_stream
            );
    }
    throw std::runtime_error("Unsupported node CG phase-one mode.");
}

bool run_node_phase_two(
    const MultiDiGraph& graph,
    const ForwardPricingContext& pricing_context,
    CGMasterProblem& master_problem,
    const NodeCGOptions& options,
    const std::chrono::steady_clock::time_point global_start_time,
    NodeCGResult& result,
    std::ostream* log_stream
) {
    switch_cg_master_to_phase_two(master_problem);
    result.phase_two_reached = true;
    if (log_stream != nullptr) {
        *log_stream << "[node-cg] Entering Phase II mode="
                    << to_string(options.phase_two_pricing_mode) << '\n';
    }
    return run_node_pricing_phase(
        graph,
        pricing_context,
        master_problem,
        options,
        CGPhase::PhaseII,
        global_start_time,
        result,
        log_stream
    );
}

}  // namespace

NodeCGResult solve_node_column_generation(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const NodeCGOptions& options,
    std::ostream* log_stream
) {
    if (options.p < 1) {
        throw std::runtime_error("Node CG requires p >= 1.");
    }
    if (!double_less_or_equal(0.0, options.time_limit)) {
        throw std::runtime_error("Node CG requires a nonnegative time limit.");
    }

    NodeCGResult result;
    result.phase_one_mode = options.phase_one_mode;
    result.phase_two_pricing_mode = options.phase_two_pricing_mode;

    const ForwardPricingContext pricing_context =
        build_forward_pricing_context(graph, options.p, options.time_limit);

    const CGMasterOptions master_options{
        options.gurobi_log_path,
        options.gurobi_threads,
        options.solver_time_limit,
    };
    CGMasterProblem master_problem =
        build_cg_master_problem(data, graph, master_options);

    const auto global_start_time = std::chrono::steady_clock::now();

    const bool phase_one_success = run_node_phase_one(
        data,
        graph,
        pricing_context,
        master_problem,
        options,
        global_start_time,
        result,
        log_stream
    );
    if (!phase_one_success || !result.phase_one_feasible) {
        capture_column_values_if_available(master_problem, result);
        return result;
    }

    const bool phase_two_success = run_node_phase_two(
        graph,
        pricing_context,
        master_problem,
        options,
        global_start_time,
        result,
        log_stream
    );
    capture_column_values_if_available(master_problem, result);
    if (!phase_two_success) {
        return result;
    }

    if (master_problem.model->get(GRB_IntAttr_SolCount) > 0) {
        result.phase_two_objective = master_problem.model->get(GRB_DoubleAttr_ObjVal);
    }
    return result;
}

void write_node_cg_summary(
    std::ostream& out,
    const NodeCGResult& result
) {
    out << "[node-cg-summary] phase_one_mode=" << to_string(result.phase_one_mode) << '\n';
    out << "[node-cg-summary] phase_two_pricing_mode="
        << to_string(result.phase_two_pricing_mode) << '\n';
    out << "[node-cg-summary] phase_one_feasible=" << (result.phase_one_feasible ? 1 : 0) << '\n';
    out << "[node-cg-summary] phase_two_reached=" << (result.phase_two_reached ? 1 : 0) << '\n';
    out << "[node-cg-summary] solved_to_completion=" << (result.solved_to_completion ? 1 : 0) << '\n';
    out << "[node-cg-summary] phase_one_objective=" << format_double(result.phase_one_objective) << '\n';
    out << "[node-cg-summary] phase_two_objective=" << format_double(result.phase_two_objective) << '\n';
    out << "[node-cg-summary] phase_one_iterations=" << result.phase_one_iterations << '\n';
    out << "[node-cg-summary] phase_two_iterations=" << result.phase_two_iterations << '\n';
    out << "[node-cg-summary] total_columns_added=" << result.total_columns_added << '\n';
    write_cg_master_snapshot(out, result.final_snapshot);
}

void write_node_cg_solution(
    std::ostream& out,
    const NodeCGResult& result
) {
    out << "Solution:\n";

    if (!result.phase_one_feasible) {
        out << "  Node CG stopped before obtaining a Phase-I feasible RMP\n";
        out << "Solution done \n";
        return;
    }

    if (!result.solved_to_completion) {
        out << "  Node CG terminated before proving LP completion\n";
        out << "Solution done \n";
        return;
    }

    for (std::size_t idx = 0; idx < result.columns.size() && idx < result.column_values.size(); ++idx) {
        if (result.column_values[idx] <= 1e-6) {
            continue;
        }

        const CGColumn& column = result.columns[idx];
        out << "  x_" << column.id
            << " = " << format_double(result.column_values[idx])
            << " | rc(gen)=" << format_double(column.reduced_cost)
            << " | tau=" << format_double(column.tau)
            << " | time=" << format_double(column.total_time)
            << " | cost=" << format_double(column.total_cost)
            << " | edges=";

        out << "[";
        for (std::size_t edge_idx = 0; edge_idx < column.edge_ids.size(); ++edge_idx) {
            if (edge_idx > 0U) {
                out << ", ";
            }
            out << column.edge_ids[edge_idx];
        }
        out << "]\n";
    }

    out << "Solution done \n";
}

const char* to_string(NodeCGPhaseOneMode mode) {
    switch (mode) {
        case NodeCGPhaseOneMode::ExactCG:
            return "exact-cg";
        case NodeCGPhaseOneMode::HeuristicCG:
            return "heuristic-cg";
    }
    return "unknown";
}

const char* to_string(NodeCGPhaseTwoPricingMode mode) {
    switch (mode) {
        case NodeCGPhaseTwoPricingMode::ExactPricing:
            return "exact-pricing";
        case NodeCGPhaseTwoPricingMode::HeuristicPricingThenExact:
            return "heuristic-pricing-then-exact";
    }
    return "unknown";
}

const char* to_string(NodeCGPhaseTwoHeuristicEngine mode) {
    switch (mode) {
        case NodeCGPhaseTwoHeuristicEngine::Labeling:
            return "labeling";
        case NodeCGPhaseTwoHeuristicEngine::ShallowSearch:
            return "shallow-search";
    }
    return "unknown";
}

}  // namespace spdp
