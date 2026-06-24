#include "PstepBnP.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iosfwd>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <ostream>
#include <random>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
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

void append_heuristic_window_column_variants(
    const MultiDiGraph& graph,
    double time_limit,
    int p,
    const std::vector<int>& edge_ids,
    std::unordered_map<std::string, std::size_t>& seen_columns,
    std::vector<CGColumn>& columns
) {
    if (p < 1 || edge_ids.empty()) {
        return;
    }

    const std::size_t max_window_length =
        std::min<std::size_t>(static_cast<std::size_t>(p), edge_ids.size());
    for (std::size_t start = 0; start < edge_ids.size(); ++start) {
        for (std::size_t length = 1; length <= max_window_length; ++length) {
            if (start + length > edge_ids.size()) {
                break;
            }
            const std::vector<int> window_edges(
                std::next(edge_ids.begin(), static_cast<long>(start)),
                std::next(edge_ids.begin(), static_cast<long>(start + length))
            );
            append_heuristic_column_variants(
                graph,
                time_limit,
                window_edges,
                seen_columns,
                columns
            );
        }
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

struct PhaseOneSeedParams {
    std::size_t max_attempts = 0;
    std::size_t max_incumbents = 20;
    std::size_t top_l = 3;
    double weight_cost = 1.0;
    double weight_time = 0.1;
    double weight_saving = 0.5;
    unsigned int random_seed = 1;
};

struct SeedPathCandidate {
    std::vector<int> edge_ids;
    std::vector<int> request_indices;
    double total_time = 0.0;
    double total_cost = 0.0;
    double fitness = 0.0;
};

struct SeedRoute {
    std::vector<int> edge_ids;
    std::vector<int> first_block_request_indices;
    double total_time = 0.0;
    double total_cost = 0.0;
};

struct SeedIncumbent {
    std::vector<SeedRoute> routes;
    double total_cost = 0.0;
    std::string key_set;
};

double edge_path_time(const MultiDiGraph& graph, const std::vector<int>& edge_ids) {
    double total = 0.0;
    for (int edge_id : edge_ids) {
        total += graph.edges()[static_cast<std::size_t>(edge_id)].data.time;
    }
    return total;
}

double edge_path_cost(const MultiDiGraph& graph, const std::vector<int>& edge_ids) {
    double total = 0.0;
    for (int edge_id : edge_ids) {
        total += graph.edges()[static_cast<std::size_t>(edge_id)].data.cost;
    }
    return total;
}

std::string edge_sequence_key(const std::vector<int>& edge_ids) {
    std::ostringstream out;
    for (std::size_t idx = 0; idx < edge_ids.size(); ++idx) {
        if (idx > 0U) {
            out << ',';
        }
        out << edge_ids[idx];
    }
    return out.str();
}

std::string node_label(const MultiDiGraph& graph, NodeId node_id) {
    if (node_id == 0) {
        return "start";
    }
    if (node_id == graph.end_node_id()) {
        return "end";
    }
    const NodeSpec& node = graph.node(node_id);
    std::ostringstream out;
    if (node.kind == NodeSpec::Kind::Pickup) {
        out << 'P';
    } else if (node.kind == NodeSpec::Kind::Delivery) {
        out << 'D';
    } else {
        out << "node";
    }
    if (node.request_idx.has_value()) {
        out << (node.request_idx.value() + 1);
    } else {
        out << node_id;
    }
    out << "(node=" << node_id << ",loc=" << node.location << ')';
    return out.str();
}

std::string route_node_path_string(
    const MultiDiGraph& graph,
    const std::vector<int>& edge_ids
) {
    if (edge_ids.empty()) {
        return "";
    }
    std::ostringstream out;
    const EdgeRecord& first_edge = graph.edges()[static_cast<std::size_t>(edge_ids.front())];
    out << node_label(graph, first_edge.u);
    for (int edge_id : edge_ids) {
        const EdgeRecord& edge = graph.edges()[static_cast<std::size_t>(edge_id)];
        out << " -> " << node_label(graph, edge.v);
    }
    return out.str();
}

void log_best_phase_one_seed_incumbent(
    std::ostream& out,
    const MultiDiGraph& graph,
    const SeedIncumbent& incumbent,
    std::size_t generated_incumbent_count,
    std::size_t selected_incumbent_count
) {
    out << "[node-cg] phase1 heuristic-cg incumbents generated="
        << generated_incumbent_count
        << " selected=" << selected_incumbent_count << '\n';
    out << "[node-cg] phase1 heuristic-cg best total_cost="
        << format_double(incumbent.total_cost)
        << " routes=" << incumbent.routes.size() << '\n';
    for (std::size_t route_idx = 0; route_idx < incumbent.routes.size(); ++route_idx) {
        const SeedRoute& route = incumbent.routes[route_idx];
        out << "[node-cg] phase1 heuristic-cg best route=" << route_idx
            << " cost=" << format_double(route.total_cost)
            << " time=" << format_double(route.total_time)
            << " edges=" << edge_sequence_key(route.edge_ids)
            << '\n';
        out << "[node-cg] phase1 heuristic-cg best route=" << route_idx
            << " path=" << route_node_path_string(graph, route.edge_ids)
            << '\n';
    }
}

void append_unique_ints(std::vector<int>& target, const std::vector<int>& values) {
    for (int value : values) {
        if (std::find(target.begin(), target.end(), value) == target.end()) {
            target.push_back(value);
        }
    }
}

std::vector<std::vector<int>> enumerate_template_edge_paths(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    const std::vector<NodeId>& node_sequence,
    const std::vector<PiRequirement>& pi_reqs,
    const NodeStateTemplateList& state_templates,
    bool prune_pickup_symmetry,
    bool prune_delivery_symmetry
) {
    std::vector<std::vector<int>> paths;
    if (prune_pickup_symmetry && violates_pickup_symmetry_43(graph, node_sequence)) {
        return paths;
    }
    if (prune_delivery_symmetry && violates_delivery_symmetry_43(graph, node_sequence)) {
        return paths;
    }
    enumerate_pricing_paths_for_node_sequence(
        graph,
        context,
        pricing_edges_by_pair,
        node_sequence,
        pi_reqs,
        state_templates,
        [&](const std::vector<int>& edge_ids) {
            paths.push_back(edge_ids);
        }
    );
    return paths;
}

std::optional<std::vector<int>> best_path_by_cost(
    const MultiDiGraph& graph,
    const std::vector<std::vector<int>>& paths
) {
    if (paths.empty()) {
        return std::nullopt;
    }
    const auto best = std::min_element(
        paths.begin(),
        paths.end(),
        [&](const std::vector<int>& lhs, const std::vector<int>& rhs) {
            return edge_path_cost(graph, lhs) < edge_path_cost(graph, rhs);
        }
    );
    return *best;
}

std::vector<int> pickup_node_ids_for_request_count(int request_count) {
    std::vector<int> result(static_cast<std::size_t>(request_count), 0);
    for (int idx = 0; idx < request_count; ++idx) {
        result[static_cast<std::size_t>(idx)] = idx + 1;
    }
    return result;
}

NodeId pickup_node_id_for_request(int request_idx) {
    return request_idx + 1;
}

NodeId delivery_node_id_for_request(int request_count, int request_idx) {
    return request_count + request_idx + 1;
}

double single_route_reference_cost(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    int request_idx,
    const NodeCGOptions& options
) {
    const int request_count = static_cast<int>(data.requests.size());
    const Request& request = data.requests.at(static_cast<std::size_t>(request_idx));
    const std::vector<std::vector<int>> paths = enumerate_template_edge_paths(
        graph,
        context,
        pricing_edges_by_pair,
        {0,
         pickup_node_id_for_request(request_idx),
         delivery_node_id_for_request(request_count, request_idx),
         graph.end_node_id()},
        {PiRequirement::Empty, PiRequirement::NonEmpty, PiRequirement::Empty},
        {{
            std::nullopt,
            make_state(make_n_token(), make_f_token(request)),
            make_state(make_n_token(), make_n_token()),
            std::nullopt,
        }},
        options.prune_pickup_symmetry_43,
        options.prune_delivery_symmetry_43
    );
    const std::optional<std::vector<int>> best = best_path_by_cost(graph, paths);
    return best.has_value() ? edge_path_cost(graph, best.value()) : 0.0;
}

double reference_saving(
    const std::vector<double>& single_cost_by_request,
    const std::vector<int>& request_indices,
    double candidate_cost
) {
    double reference = 0.0;
    for (int request_idx : request_indices) {
        reference += single_cost_by_request.at(static_cast<std::size_t>(request_idx));
    }
    return reference - candidate_cost;
}

void score_seed_candidate(
    SeedPathCandidate& candidate,
    const PhaseOneSeedParams& params,
    const std::vector<double>& single_cost_by_request
) {
    const double saving = reference_saving(
        single_cost_by_request,
        candidate.request_indices,
        candidate.total_cost
    );
    candidate.fitness =
        params.weight_cost * candidate.total_cost +
        params.weight_time * candidate.total_time -
        params.weight_saving * saving;
}

std::optional<SeedPathCandidate> weighted_select_best_l(
    std::vector<SeedPathCandidate> candidates,
    const PhaseOneSeedParams& params,
    std::mt19937& rng
) {
    if (candidates.empty()) {
        return std::nullopt;
    }
    std::sort(
        candidates.begin(),
        candidates.end(),
        [](const SeedPathCandidate& lhs, const SeedPathCandidate& rhs) {
            return lhs.fitness < rhs.fitness;
        }
    );
    const std::size_t keep = std::min(params.top_l, candidates.size());
    double worst_kept_fitness = candidates[keep - 1U].fitness;
    std::vector<double> weights;
    weights.reserve(keep);
    for (std::size_t idx = 0; idx < keep; ++idx) {
        weights.push_back(std::max(1e-9, worst_kept_fitness - candidates[idx].fitness + 1e-9));
    }
    std::discrete_distribution<std::size_t> distribution(weights.begin(), weights.end());
    return candidates[distribution(rng)];
}

bool request_uncovered(const std::vector<bool>& covered, int request_idx) {
    return !covered.at(static_cast<std::size_t>(request_idx));
}

std::vector<int> uncovered_request_indices(const std::vector<bool>& covered) {
    std::vector<int> result;
    for (std::size_t idx = 0; idx < covered.size(); ++idx) {
        if (!covered[idx]) {
            result.push_back(static_cast<int>(idx));
        }
    }
    return result;
}

std::optional<std::vector<int>> best_close_path(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    NodeId current_node,
    State current_state
) {
    const std::vector<std::vector<int>> close_paths = enumerate_template_edge_paths(
        graph,
        context,
        pricing_edges_by_pair,
        {current_node, graph.end_node_id()},
        {PiRequirement::Empty},
        {{canonicalize_state(current_state), std::nullopt}},
        false,
        false
    );
    return best_path_by_cost(graph, close_paths);
}

std::vector<SeedPathCandidate> generate_fixed_block_candidates(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    const NodeCGOptions& options,
    const PhaseOneSeedParams& params,
    const std::vector<double>& single_cost_by_request,
    const std::vector<bool>& local_available,
    NodeId current_node,
    State current_state,
    double route_time,
    std::optional<int> target_request
) {
    const int request_count = static_cast<int>(data.requests.size());
    std::vector<SeedPathCandidate> candidates;
    for (int i = 0; i < request_count; ++i) {
        if (!request_uncovered(local_available, i)) {
            continue;
        }
        for (int j = 0; j < request_count; ++j) {
            if (j == i || !request_uncovered(local_available, j)) {
                continue;
            }
            if (target_request.has_value() && i != target_request.value() &&
                j != target_request.value()) {
                continue;
            }
            const Request& request_i = data.requests.at(static_cast<std::size_t>(i));
            const Request& request_j = data.requests.at(static_cast<std::size_t>(j));
            const NodeId pickup_i = pickup_node_id_for_request(i);
            const NodeId pickup_j = pickup_node_id_for_request(j);
            const NodeId delivery_i = delivery_node_id_for_request(request_count, i);
            const NodeId delivery_j = delivery_node_id_for_request(request_count, j);
            const std::vector<std::pair<NodeId, NodeId>> delivery_orders = {
                {delivery_i, delivery_j},
                {delivery_j, delivery_i},
            };
            for (const auto& delivery_order : delivery_orders) {
                const int delivered_first_idx =
                    request_index_for_node(graph, delivery_order.first);
                const State first_delivery_state = fixed_pair_post_first_delivery_state(
                    request_i,
                    request_j,
                    i,
                    delivered_first_idx
                );
                const std::vector<std::vector<int>> paths = enumerate_template_edge_paths(
                    graph,
                    context,
                    pricing_edges_by_pair,
                    {current_node, pickup_i, pickup_j, delivery_order.first, delivery_order.second},
                    {PiRequirement::Empty,
                     PiRequirement::Empty,
                     PiRequirement::NonEmpty,
                     PiRequirement::Empty},
                    {{
                        canonicalize_state(current_state),
                        make_state(make_n_token(), make_f_token(request_i)),
                        make_state(make_f_token(request_i), make_f_token(request_j)),
                        first_delivery_state,
                        make_state(make_n_token(), make_n_token()),
                    }},
                    options.prune_pickup_symmetry_43,
                    options.prune_delivery_symmetry_43
                );
                for (const std::vector<int>& path : paths) {
                    const double path_time = edge_path_time(graph, path);
                    const auto close_path = best_close_path(
                        graph,
                        context,
                        pricing_edges_by_pair,
                        delivery_order.second,
                        make_state(make_n_token(), make_n_token())
                    );
                    if (!close_path.has_value()) {
                        continue;
                    }
                    if (!double_less_or_equal(
                            route_time + path_time + edge_path_time(graph, close_path.value()),
                            options.time_limit
                        )) {
                        continue;
                    }
                    SeedPathCandidate candidate;
                    candidate.edge_ids = path;
                    candidate.request_indices = {i, j};
                    candidate.total_time = path_time;
                    candidate.total_cost = edge_path_cost(graph, path);
                    score_seed_candidate(candidate, params, single_cost_by_request);
                    candidates.push_back(std::move(candidate));
                }
            }
        }
    }
    return candidates;
}

std::vector<SeedPathCandidate> generate_flexible_front_candidates(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    const NodeCGOptions& options,
    const PhaseOneSeedParams& params,
    const std::vector<double>& single_cost_by_request,
    const std::vector<bool>& local_available,
    std::optional<int> target_request
) {
    const int request_count = static_cast<int>(data.requests.size());
    std::vector<SeedPathCandidate> candidates;
    for (int i = 0; i < request_count; ++i) {
        if (!request_uncovered(local_available, i)) {
            continue;
        }
        for (int j = 0; j < request_count; ++j) {
            if (j == i || !request_uncovered(local_available, j)) {
                continue;
            }
            if (target_request.has_value() && i != target_request.value() &&
                j != target_request.value()) {
                continue;
            }
            const Request& request_i = data.requests.at(static_cast<std::size_t>(i));
            const Request& request_j = data.requests.at(static_cast<std::size_t>(j));
            const std::vector<std::vector<int>> paths = enumerate_template_edge_paths(
                graph,
                context,
                pricing_edges_by_pair,
                {0, pickup_node_id_for_request(i), pickup_node_id_for_request(j)},
                {PiRequirement::Empty, PiRequirement::Empty},
                {{
                    std::nullopt,
                    make_state(make_n_token(), make_f_token(request_i)),
                    make_state(make_f_token(request_i), make_f_token(request_j)),
                }},
                options.prune_pickup_symmetry_43,
                options.prune_delivery_symmetry_43
            );
            for (const std::vector<int>& path : paths) {
                SeedPathCandidate candidate;
                candidate.edge_ids = path;
                candidate.request_indices = {i, j};
                candidate.total_time = edge_path_time(graph, path);
                candidate.total_cost = edge_path_cost(graph, path);
                score_seed_candidate(candidate, params, single_cost_by_request);
                candidates.push_back(std::move(candidate));
            }
        }
    }
    return candidates;
}

std::vector<SeedPathCandidate> generate_flexible_middle_candidates(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    const NodeCGOptions& options,
    const PhaseOneSeedParams& params,
    const std::vector<double>& single_cost_by_request,
    const std::vector<bool>& local_available,
    NodeId current_node,
    State current_state,
    double route_time,
    int exit_a,
    int exit_b
) {
    const int request_count = static_cast<int>(data.requests.size());
    std::vector<SeedPathCandidate> candidates;
    for (int i = 0; i < request_count; ++i) {
        if (!request_uncovered(local_available, i)) {
            continue;
        }
        for (int j = 0; j < request_count; ++j) {
            if (j == i || !request_uncovered(local_available, j)) {
                continue;
            }
            const Request& request_i = data.requests.at(static_cast<std::size_t>(i));
            const Request& request_j = data.requests.at(static_cast<std::size_t>(j));
            const NodeId delivery_i = delivery_node_id_for_request(request_count, i);
            const NodeId pickup_i = pickup_node_id_for_request(i);
            const NodeId delivery_j = delivery_node_id_for_request(request_count, j);
            const NodeId pickup_j = pickup_node_id_for_request(j);
            const std::vector<std::vector<int>> paths = enumerate_template_edge_paths(
                graph,
                context,
                pricing_edges_by_pair,
                {current_node, delivery_i, pickup_i, delivery_j, pickup_j},
                {PiRequirement::NonEmpty,
                 PiRequirement::Empty,
                 PiRequirement::Empty,
                 PiRequirement::Empty},
                {{
                    canonicalize_state(current_state),
                    make_state(make_n_token(), make_e_token(request_j.container_type)),
                    make_state(make_f_token(request_i), make_e_token(request_j.container_type)),
                    make_state(make_f_token(request_i), make_n_token()),
                    make_state(make_f_token(request_i), make_f_token(request_j)),
                }},
                options.prune_pickup_symmetry_43,
                options.prune_delivery_symmetry_43
            );
            for (const std::vector<int>& path : paths) {
                const double path_time = edge_path_time(graph, path);
                const State next_state = make_state(make_f_token(request_i), make_f_token(request_j));
                const std::optional<std::vector<int>> exit_path = best_close_path(
                    graph,
                    context,
                    pricing_edges_by_pair,
                    pickup_j,
                    next_state
                );
                (void)exit_path;
                const std::vector<int> exit_requests = {exit_a, exit_b};
                bool can_close = false;
                for (int first_exit : exit_requests) {
                    const int second_exit = first_exit == exit_a ? exit_b : exit_a;
                    const NodeId delivery_first =
                        delivery_node_id_for_request(request_count, first_exit);
                    const NodeId delivery_second =
                        delivery_node_id_for_request(request_count, second_exit);
                    const Request& second_request =
                        data.requests.at(static_cast<std::size_t>(second_exit));
                    const std::vector<std::vector<int>> close_paths = enumerate_template_edge_paths(
                        graph,
                        context,
                        pricing_edges_by_pair,
                        {pickup_j, delivery_first, delivery_second, graph.end_node_id()},
                        {PiRequirement::NonEmpty,
                         PiRequirement::Empty,
                         PiRequirement::Empty},
                        {{
                            next_state,
                            make_state(make_n_token(), make_e_token(second_request.container_type)),
                            make_state(make_n_token(), make_n_token()),
                            std::nullopt,
                        }},
                        options.prune_pickup_symmetry_43,
                        options.prune_delivery_symmetry_43
                    );
                    for (const std::vector<int>& close_path : close_paths) {
                        if (double_less_or_equal(
                                route_time + path_time + edge_path_time(graph, close_path),
                                options.time_limit
                            )) {
                            can_close = true;
                            break;
                        }
                    }
                    if (can_close) {
                        break;
                    }
                }
                if (!can_close) {
                    continue;
                }
                SeedPathCandidate candidate;
                candidate.edge_ids = path;
                candidate.request_indices = {i, j};
                candidate.total_time = path_time;
                candidate.total_cost = edge_path_cost(graph, path);
                score_seed_candidate(candidate, params, single_cost_by_request);
                candidates.push_back(std::move(candidate));
            }
        }
    }
    return candidates;
}

std::optional<std::vector<int>> best_flexible_exit_path(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    const NodeCGOptions& options,
    NodeId current_node,
    State current_state,
    int exit_a,
    int exit_b
) {
    const int request_count = static_cast<int>(data.requests.size());
    std::vector<std::vector<int>> candidates;
    for (int first_exit : {exit_a, exit_b}) {
        const int second_exit = first_exit == exit_a ? exit_b : exit_a;
        const Request& second_request = data.requests.at(static_cast<std::size_t>(second_exit));
        std::vector<std::vector<int>> paths = enumerate_template_edge_paths(
            graph,
            context,
            pricing_edges_by_pair,
            {current_node,
             delivery_node_id_for_request(request_count, first_exit),
             delivery_node_id_for_request(request_count, second_exit),
             graph.end_node_id()},
            {PiRequirement::NonEmpty, PiRequirement::Empty, PiRequirement::Empty},
            {{
                canonicalize_state(current_state),
                make_state(make_n_token(), make_e_token(second_request.container_type)),
                make_state(make_n_token(), make_n_token()),
                std::nullopt,
            }},
            options.prune_pickup_symmetry_43,
            options.prune_delivery_symmetry_43
        );
        candidates.insert(candidates.end(), paths.begin(), paths.end());
    }
    return best_path_by_cost(graph, candidates);
}

std::optional<SeedRoute> build_flexible_seed_route(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    const NodeCGOptions& options,
    const PhaseOneSeedParams& params,
    const std::vector<double>& single_cost_by_request,
    const std::vector<bool>& uncovered,
    std::optional<int> target_request,
    std::mt19937& rng
) {
    if (uncovered_request_indices(uncovered).size() < 4U) {
        return std::nullopt;
    }
    std::vector<SeedPathCandidate> front_candidates = generate_flexible_front_candidates(
        data,
        graph,
        context,
        pricing_edges_by_pair,
        options,
        params,
        single_cost_by_request,
        uncovered,
        target_request
    );
    std::optional<SeedPathCandidate> selected_front =
        weighted_select_best_l(std::move(front_candidates), params, rng);
    if (!selected_front.has_value()) {
        return std::nullopt;
    }

    SeedRoute route;
    route.edge_ids = selected_front->edge_ids;
    route.first_block_request_indices = selected_front->request_indices;
    route.total_time = selected_front->total_time;
    route.total_cost = selected_front->total_cost;

    std::vector<bool> local_available = uncovered;
    for (int request_idx : selected_front->request_indices) {
        local_available[static_cast<std::size_t>(request_idx)] = true;
    }
    const int exit_a = selected_front->request_indices[0];
    const int exit_b = selected_front->request_indices[1];

    NodeId current_node =
        graph.edges()[static_cast<std::size_t>(route.edge_ids.back())].v;
    State current_state =
        canonicalize_state(graph.edges()[static_cast<std::size_t>(route.edge_ids.back())]
                               .data.end_state);
    std::size_t middle_count = 0;

    while (uncovered_request_indices(local_available).size() >= 2U) {
        std::vector<SeedPathCandidate> middle_candidates = generate_flexible_middle_candidates(
            data,
            graph,
            context,
            pricing_edges_by_pair,
            options,
            params,
            single_cost_by_request,
            local_available,
            current_node,
            current_state,
            route.total_time,
            exit_a,
            exit_b
        );
        std::optional<SeedPathCandidate> selected_middle =
            weighted_select_best_l(std::move(middle_candidates), params, rng);
        if (!selected_middle.has_value()) {
            break;
        }
        route.edge_ids.insert(
            route.edge_ids.end(),
            selected_middle->edge_ids.begin(),
            selected_middle->edge_ids.end()
        );
        route.total_time += selected_middle->total_time;
        route.total_cost += selected_middle->total_cost;
        for (int request_idx : selected_middle->request_indices) {
            local_available[static_cast<std::size_t>(request_idx)] = true;
        }
        current_node =
            graph.edges()[static_cast<std::size_t>(route.edge_ids.back())].v;
        current_state =
            canonicalize_state(graph.edges()[static_cast<std::size_t>(route.edge_ids.back())]
                                   .data.end_state);
        ++middle_count;
    }

    if (middle_count == 0U) {
        return std::nullopt;
    }

    std::optional<std::vector<int>> exit_path = best_flexible_exit_path(
        data,
        graph,
        context,
        pricing_edges_by_pair,
        options,
        current_node,
        current_state,
        exit_a,
        exit_b
    );
    if (!exit_path.has_value()) {
        return std::nullopt;
    }
    const double exit_time = edge_path_time(graph, exit_path.value());
    if (!double_less_or_equal(route.total_time + exit_time, options.time_limit)) {
        return std::nullopt;
    }
    route.total_time += exit_time;
    route.total_cost += edge_path_cost(graph, exit_path.value());
    route.edge_ids.insert(route.edge_ids.end(), exit_path->begin(), exit_path->end());
    return route;
}

std::optional<SeedRoute> build_fixed_seed_route(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    const NodeCGOptions& options,
    const PhaseOneSeedParams& params,
    const std::vector<double>& single_cost_by_request,
    const std::vector<bool>& uncovered,
    std::optional<int> target_request,
    std::mt19937& rng
) {
    if (uncovered_request_indices(uncovered).size() < 2U) {
        return std::nullopt;
    }
    SeedRoute route;
    NodeId current_node = 0;
    State current_state = make_state(make_n_token(), make_n_token());
    std::vector<bool> local_available = uncovered;
    std::size_t block_count = 0;
    std::optional<int> force_target = target_request;

    while (uncovered_request_indices(local_available).size() >= 2U) {
        std::vector<SeedPathCandidate> candidates = generate_fixed_block_candidates(
            data,
            graph,
            context,
            pricing_edges_by_pair,
            options,
            params,
            single_cost_by_request,
            local_available,
            current_node,
            current_state,
            route.total_time,
            force_target
        );
        std::optional<SeedPathCandidate> selected =
            weighted_select_best_l(std::move(candidates), params, rng);
        if (!selected.has_value()) {
            break;
        }
        if (block_count == 0U) {
            route.first_block_request_indices = selected->request_indices;
        }
        route.edge_ids.insert(route.edge_ids.end(), selected->edge_ids.begin(), selected->edge_ids.end());
        route.total_time += selected->total_time;
        route.total_cost += selected->total_cost;
        for (int request_idx : selected->request_indices) {
            local_available[static_cast<std::size_t>(request_idx)] = true;
        }
        current_node =
            graph.edges()[static_cast<std::size_t>(route.edge_ids.back())].v;
        current_state =
            canonicalize_state(graph.edges()[static_cast<std::size_t>(route.edge_ids.back())]
                                   .data.end_state);
        force_target = std::nullopt;
        ++block_count;
    }

    if (block_count == 0U) {
        return std::nullopt;
    }
    std::optional<std::vector<int>> close_path = best_close_path(
        graph,
        context,
        pricing_edges_by_pair,
        current_node,
        current_state
    );
    if (!close_path.has_value()) {
        return std::nullopt;
    }
    const double close_time = edge_path_time(graph, close_path.value());
    if (!double_less_or_equal(route.total_time + close_time, options.time_limit)) {
        return std::nullopt;
    }
    route.total_time += close_time;
    route.total_cost += edge_path_cost(graph, close_path.value());
    route.edge_ids.insert(route.edge_ids.end(), close_path->begin(), close_path->end());
    return route;
}

std::optional<std::vector<int>> best_single_standalone_route(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    const NodeCGOptions& options,
    int request_idx
) {
    const int request_count = static_cast<int>(data.requests.size());
    const Request& request = data.requests.at(static_cast<std::size_t>(request_idx));
    const std::vector<std::vector<int>> paths = enumerate_template_edge_paths(
        graph,
        context,
        pricing_edges_by_pair,
        {0,
         pickup_node_id_for_request(request_idx),
         delivery_node_id_for_request(request_count, request_idx),
         graph.end_node_id()},
        {PiRequirement::Empty, PiRequirement::NonEmpty, PiRequirement::Empty},
        {{
            std::nullopt,
            make_state(make_n_token(), make_f_token(request)),
            make_state(make_n_token(), make_n_token()),
            std::nullopt,
        }},
        options.prune_pickup_symmetry_43,
        options.prune_delivery_symmetry_43
    );
    return best_path_by_cost(graph, paths);
}

struct SingleInsertionCandidate {
    std::size_t route_index = 0;
    double delta_cost = std::numeric_limits<double>::infinity();
    double total_time = 0.0;
    double total_cost = 0.0;
    std::vector<int> edge_ids;
};

std::optional<SingleInsertionCandidate> best_single_insertion_at_nn_boundary(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    const NodeCGOptions& options,
    int request_idx,
    std::size_t route_index,
    const SeedRoute& route
) {
    const int request_count = static_cast<int>(data.requests.size());
    const Request& request = data.requests.at(static_cast<std::size_t>(request_idx));
    SingleInsertionCandidate best;
    best.route_index = route_index;

    for (std::size_t edge_pos = 0; edge_pos < route.edge_ids.size(); ++edge_pos) {
        const EdgeRecord& replaced_edge =
            graph.edges()[static_cast<std::size_t>(route.edge_ids[edge_pos])];
        if (canonicalize_state(replaced_edge.data.start_state) !=
            make_state(make_n_token(), make_n_token())) {
            continue;
        }
        if (!replaced_edge.data.sequence_pi.empty()) {
            continue;
        }
        const std::vector<std::vector<int>> insertion_paths = enumerate_template_edge_paths(
            graph,
            context,
            pricing_edges_by_pair,
            {replaced_edge.u,
             pickup_node_id_for_request(request_idx),
             delivery_node_id_for_request(request_count, request_idx),
             replaced_edge.v},
            {PiRequirement::Empty, PiRequirement::NonEmpty, PiRequirement::Empty},
            {{
                canonicalize_state(replaced_edge.data.start_state),
                make_state(make_n_token(), make_f_token(request)),
                make_state(make_n_token(), make_n_token()),
                canonicalize_state(replaced_edge.data.end_state),
            }},
            options.prune_pickup_symmetry_43,
            options.prune_delivery_symmetry_43
        );
        for (const std::vector<int>& insertion : insertion_paths) {
            const double new_time =
                route.total_time - replaced_edge.data.time + edge_path_time(graph, insertion);
            if (!double_less_or_equal(new_time, options.time_limit)) {
                continue;
            }
            const double delta_cost = edge_path_cost(graph, insertion) - replaced_edge.data.cost;
            if (delta_cost < best.delta_cost) {
                best.delta_cost = delta_cost;
                best.edge_ids = route.edge_ids;
                best.edge_ids.erase(best.edge_ids.begin() + static_cast<long>(edge_pos));
                best.edge_ids.insert(
                    best.edge_ids.begin() + static_cast<long>(edge_pos),
                    insertion.begin(),
                    insertion.end()
                );
            }
        }
    }

    if (!std::isfinite(best.delta_cost)) {
        return std::nullopt;
    }
    best.total_time = edge_path_time(graph, best.edge_ids);
    best.total_cost = edge_path_cost(graph, best.edge_ids);
    return best;
}

void insert_residual_single_routes(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    const NodeCGOptions& options,
    std::vector<bool>& uncovered,
    std::vector<SeedRoute>& routes
) {
    const int request_count = static_cast<int>(data.requests.size());
    for (int request_idx = 0; request_idx < request_count; ++request_idx) {
        if (!request_uncovered(uncovered, request_idx)) {
            continue;
        }
        std::optional<SingleInsertionCandidate> best_insertion;
        for (std::size_t route_index = 0; route_index < routes.size(); ++route_index) {
            std::optional<SingleInsertionCandidate> insertion =
                best_single_insertion_at_nn_boundary(
                    data,
                    graph,
                    context,
                    pricing_edges_by_pair,
                    options,
                    request_idx,
                    route_index,
                    routes[route_index]
                );
            if (insertion.has_value() &&
                (!best_insertion.has_value() ||
                 insertion->delta_cost < best_insertion->delta_cost)) {
                best_insertion = std::move(insertion.value());
            }
        }
        if (best_insertion.has_value()) {
            SeedRoute& route = routes[best_insertion->route_index];
            route.edge_ids = std::move(best_insertion->edge_ids);
            route.total_time = best_insertion->total_time;
            route.total_cost = best_insertion->total_cost;
        } else {
            std::optional<std::vector<int>> standalone = best_single_standalone_route(
                data,
                graph,
                context,
                pricing_edges_by_pair,
                options,
                request_idx
            );
            if (!standalone.has_value()) {
                continue;
            }
            SeedRoute route;
            route.edge_ids = standalone.value();
            route.first_block_request_indices = {request_idx};
            route.total_time = edge_path_time(graph, route.edge_ids);
            route.total_cost = edge_path_cost(graph, route.edge_ids);
            routes.push_back(std::move(route));
        }
        uncovered[static_cast<std::size_t>(request_idx)] = true;
    }
}

enum class PhaseOneSeedOrder {
    FlexibleFirst,
    FixedFirst,
};

void build_routes_by_type(
    PhaseOneSeedOrder order_type,
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    const NodeCGOptions& options,
    const PhaseOneSeedParams& params,
    const std::vector<double>& single_cost_by_request,
    std::vector<bool>& uncovered,
    std::vector<SeedRoute>& routes,
    std::optional<int> target_request,
    std::vector<int>& first_block_requests,
    std::mt19937& rng
) {
    std::optional<int> forced_target = target_request;
    while (true) {
        std::optional<SeedRoute> route =
            order_type == PhaseOneSeedOrder::FlexibleFirst
                ? build_flexible_seed_route(
                      data,
                      graph,
                      context,
                      pricing_edges_by_pair,
                      options,
                      params,
                      single_cost_by_request,
                      uncovered,
                      forced_target,
                      rng
                  )
                : build_fixed_seed_route(
                      data,
                      graph,
                      context,
                      pricing_edges_by_pair,
                      options,
                      params,
                      single_cost_by_request,
                      uncovered,
                      forced_target,
                      rng
                  );
        if (!route.has_value()) {
            break;
        }
        if (first_block_requests.empty()) {
            first_block_requests = route->first_block_request_indices;
        }
        for (int request_idx : route->first_block_request_indices) {
            (void)request_idx;
        }
        for (int edge_id : route->edge_ids) {
            const EdgeRecord& edge = graph.edges()[static_cast<std::size_t>(edge_id)];
            if (graph.is_physical_service_node(edge.u)) {
                uncovered[static_cast<std::size_t>(request_index_for_node(graph, edge.u))] = true;
            }
            if (graph.is_physical_service_node(edge.v)) {
                uncovered[static_cast<std::size_t>(request_index_for_node(graph, edge.v))] = true;
            }
        }
        routes.push_back(std::move(route.value()));
        forced_target = std::nullopt;
    }
}

std::optional<SeedIncumbent> build_seed_incumbent(
    PhaseOneSeedOrder order,
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const PricingEdgePairMap& pricing_edges_by_pair,
    const NodeCGOptions& options,
    const PhaseOneSeedParams& params,
    const std::vector<double>& single_cost_by_request,
    int target_request,
    std::vector<int>& first_block_requests,
    std::mt19937& rng
) {
    std::vector<bool> uncovered(data.requests.size(), false);
    std::vector<SeedRoute> routes;

    if (order == PhaseOneSeedOrder::FlexibleFirst) {
        build_routes_by_type(
            PhaseOneSeedOrder::FlexibleFirst,
            data,
            graph,
            context,
            pricing_edges_by_pair,
            options,
            params,
            single_cost_by_request,
            uncovered,
            routes,
            target_request,
            first_block_requests,
            rng
        );
        build_routes_by_type(
            PhaseOneSeedOrder::FixedFirst,
            data,
            graph,
            context,
            pricing_edges_by_pair,
            options,
            params,
            single_cost_by_request,
            uncovered,
            routes,
            std::nullopt,
            first_block_requests,
            rng
        );
    } else {
        build_routes_by_type(
            PhaseOneSeedOrder::FixedFirst,
            data,
            graph,
            context,
            pricing_edges_by_pair,
            options,
            params,
            single_cost_by_request,
            uncovered,
            routes,
            target_request,
            first_block_requests,
            rng
        );
        build_routes_by_type(
            PhaseOneSeedOrder::FlexibleFirst,
            data,
            graph,
            context,
            pricing_edges_by_pair,
            options,
            params,
            single_cost_by_request,
            uncovered,
            routes,
            std::nullopt,
            first_block_requests,
            rng
        );
    }

    insert_residual_single_routes(
        data,
        graph,
        context,
        pricing_edges_by_pair,
        options,
        uncovered,
        routes
    );

    if (std::any_of(uncovered.begin(), uncovered.end(), [](bool value) { return !value; })) {
        return std::nullopt;
    }
    SeedIncumbent incumbent;
    incumbent.routes = std::move(routes);
    for (const SeedRoute& route : incumbent.routes) {
        incumbent.total_cost += route.total_cost;
    }
    return incumbent;
}

void append_backward_decomposed_route_columns(
    const MultiDiGraph& graph,
    const NodeCGOptions& options,
    const SeedRoute& route,
    std::unordered_map<std::string, std::size_t>& seen_columns,
    std::vector<CGColumn>& columns
) {
    if (route.edge_ids.empty()) {
        return;
    }
    std::size_t end = route.edge_ids.size();
    std::vector<std::vector<int>> segments_reversed;
    while (end > static_cast<std::size_t>(options.p)) {
        const std::size_t start = end - static_cast<std::size_t>(options.p);
        segments_reversed.emplace_back(route.edge_ids.begin() + static_cast<long>(start),
                                       route.edge_ids.begin() + static_cast<long>(end));
        end = start;
    }
    segments_reversed.emplace_back(route.edge_ids.begin(), route.edge_ids.begin() + static_cast<long>(end));
    for (auto iter = segments_reversed.rbegin(); iter != segments_reversed.rend(); ++iter) {
        append_heuristic_column_variants(
            graph,
            options.time_limit,
            *iter,
            seen_columns,
            columns
        );
    }
}

std::string incumbent_decomposition_key(
    const MultiDiGraph& graph,
    const NodeCGOptions& options,
    const SeedIncumbent& incumbent
) {
    std::vector<std::string> route_keys;
    for (const SeedRoute& route : incumbent.routes) {
        std::size_t end = route.edge_ids.size();
        while (end > static_cast<std::size_t>(options.p)) {
            const std::size_t start = end - static_cast<std::size_t>(options.p);
            std::vector<int> segment(route.edge_ids.begin() + static_cast<long>(start),
                                     route.edge_ids.begin() + static_cast<long>(end));
            route_keys.push_back(edge_sequence_key(segment));
            end = start;
        }
        std::vector<int> segment(route.edge_ids.begin(), route.edge_ids.begin() + static_cast<long>(end));
        route_keys.push_back(edge_sequence_key(segment));
    }
    (void)graph;
    std::sort(route_keys.begin(), route_keys.end());
    std::ostringstream out;
    for (const std::string& key : route_keys) {
        out << key << ';';
    }
    return out.str();
}

std::vector<CGColumn> build_phase_one_route_heuristic_columns(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const NodeCGOptions& options,
    std::ostream* log_stream
) {
    if (options.p < 1) {
        throw std::runtime_error("heuristic-cg requires p >= 1.");
    }

    const PricingEdgePairMap pricing_edges_by_pair = build_pricing_edge_pair_map(graph, context);
    PhaseOneSeedParams params;
    params.max_attempts = options.phase_one_heuristic_cg_max_attempts;
    if (params.max_attempts == 0U) {
        params.max_attempts = 2U * data.requests.size();
    }
    params.max_incumbents = options.phase_one_heuristic_cg_max_incumbents;
    params.top_l = std::max<std::size_t>(1U, options.phase_one_heuristic_cg_top_l);
    params.weight_cost = options.phase_one_heuristic_cg_weight_cost;
    params.weight_time = options.phase_one_heuristic_cg_weight_time;
    params.weight_saving = options.phase_one_heuristic_cg_weight_saving;
    params.random_seed = options.phase_one_heuristic_cg_random_seed;

    std::vector<double> single_cost_by_request(data.requests.size(), 0.0);
    for (std::size_t idx = 0; idx < data.requests.size(); ++idx) {
        single_cost_by_request[idx] = single_route_reference_cost(
            data,
            graph,
            context,
            pricing_edges_by_pair,
            static_cast<int>(idx),
            options
        );
    }

    std::mt19937 rng(params.random_seed);
    std::vector<bool> first_block_seen(data.requests.size(), false);
    std::unordered_map<std::string, SeedIncumbent> incumbent_by_key;

    for (std::size_t attempt = 0; attempt < params.max_attempts; ++attempt) {
        int target = -1;
        for (std::size_t idx = 0; idx < first_block_seen.size(); ++idx) {
            if (!first_block_seen[idx]) {
                target = static_cast<int>(idx);
                break;
            }
        }
        if (target < 0) {
            target = static_cast<int>(attempt % std::max<std::size_t>(1U, data.requests.size()));
        }

        for (PhaseOneSeedOrder order : {PhaseOneSeedOrder::FlexibleFirst,
                                        PhaseOneSeedOrder::FixedFirst}) {
            std::vector<int> first_block_requests;
            std::optional<SeedIncumbent> incumbent = build_seed_incumbent(
                order,
                data,
                graph,
                context,
                pricing_edges_by_pair,
                options,
                params,
                single_cost_by_request,
                target,
                first_block_requests,
                rng
            );
            for (int request_idx : first_block_requests) {
                if (request_idx >= 0 &&
                    request_idx < static_cast<int>(first_block_seen.size())) {
                    first_block_seen[static_cast<std::size_t>(request_idx)] = true;
                }
            }
            if (!incumbent.has_value()) {
                continue;
            }
            const std::string key = incumbent_decomposition_key(graph, options, incumbent.value());
            auto found = incumbent_by_key.find(key);
            if (found == incumbent_by_key.end() ||
                incumbent->total_cost < found->second.total_cost) {
                incumbent->key_set = key;
                incumbent_by_key[key] = std::move(incumbent.value());
            }
        }

        if (std::all_of(first_block_seen.begin(), first_block_seen.end(), [](bool value) {
                return value;
            })) {
            break;
        }
    }

    std::vector<SeedIncumbent> incumbents;
    incumbents.reserve(incumbent_by_key.size());
    for (auto& entry : incumbent_by_key) {
        incumbents.push_back(std::move(entry.second));
    }
    const std::size_t generated_incumbent_count = incumbents.size();
    std::sort(
        incumbents.begin(),
        incumbents.end(),
        [](const SeedIncumbent& lhs, const SeedIncumbent& rhs) {
            return lhs.total_cost < rhs.total_cost;
        }
    );
    if (incumbents.size() > params.max_incumbents) {
        incumbents.resize(params.max_incumbents);
    }
    if (log_stream != nullptr) {
        if (incumbents.empty()) {
            *log_stream << "[node-cg] phase1 heuristic-cg incumbents generated=0 selected=0\n";
        } else {
            log_best_phase_one_seed_incumbent(
                *log_stream,
                graph,
                incumbents.front(),
                generated_incumbent_count,
                incumbents.size()
            );
        }
    }

    std::unordered_map<std::string, std::size_t> seen_columns;
    std::vector<CGColumn> columns;
    for (const SeedIncumbent& incumbent : incumbents) {
        for (const SeedRoute& route : incumbent.routes) {
            append_backward_decomposed_route_columns(
                graph,
                options,
                route,
                seen_columns,
                columns
            );
        }
    }
    return columns;
}

std::vector<CGColumn> build_phase_one_heuristic_columns_from_patterns(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const NodeCGOptions& options,
    bool decompose_generated_paths
) {
    if (decompose_generated_paths && options.p < 1) {
        throw std::runtime_error("heuristic-cg requires p >= 1.");
    }
    if (!decompose_generated_paths && options.p != 3) {
        throw std::runtime_error("heuristic-cg-3-step currently requires p = 3.");
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
                    if (decompose_generated_paths) {
                        append_heuristic_window_column_variants(
                            graph,
                            options.time_limit,
                            options.p,
                            edge_ids,
                            seen_columns,
                            columns
                        );
                    } else {
                        append_heuristic_column_variants(
                            graph,
                            options.time_limit,
                            edge_ids,
                            seen_columns,
                            columns
                        );
                    }
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

std::vector<CGColumn> build_phase_one_heuristic_columns(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const NodeCGOptions& options,
    std::ostream* log_stream
) {
    return build_phase_one_route_heuristic_columns(
        data,
        graph,
        context,
        options,
        log_stream
    );
}

std::vector<CGColumn> build_phase_one_heuristic_columns_3_step(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const NodeCGOptions& options
) {
    return build_phase_one_heuristic_columns_from_patterns(
        data,
        graph,
        context,
        options,
        false
    );
}

struct PhaseTwoColumnPoolEntry {
    CGColumn column;
    double best_reduced_cost = std::numeric_limits<double>::infinity();
};

std::size_t compute_search_column_limit(
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
        << " heuristic_available_starts=" << log.heuristic_available_start_count
        << " heuristic_pricing_starts=" << log.heuristic_explored_start_count
        << " heuristic_pricing_cols=" << log.heuristic_found_column_count
        << " heuristic_engine=" << to_string(log.heuristic_engine)
        << " heuristic_max_starts=" << log.heuristic_max_starts
        << " heuristic_start_ratio=" << format_double(log.heuristic_start_ratio)
        << " heuristic_ladder_levels=" << log.heuristic_ladder_levels
        << " heuristic_output_cols_per_start=" << log.heuristic_output_max_columns_per_start
        << " heuristic_output_cols_total=" << log.heuristic_output_max_columns_total
        << " heuristic_search_col_ratio=" << format_double(log.heuristic_search_column_ratio)
        << " heuristic_search_cols_per_start="
        << log.heuristic_effective_search_max_columns_per_start
        << " heuristic_search_cols_total="
        << log.heuristic_effective_search_max_columns_total
        << " heuristic_total_negative_cols=" << log.heuristic_total_negative_column_count
        << " heuristic_per_start_cap_hits=" << log.heuristic_per_start_search_cap_hit_count
        << " heuristic_global_cap_hit=" << (log.heuristic_global_search_cap_hit ? 1 : 0)
        << " heuristic_labeling_top_k_next=" << log.heuristic_labeling_top_k_next
        << " heuristic_shallow_k1=" << log.heuristic_shallow_k1
        << " heuristic_shallow_k2=" << log.heuristic_shallow_k2
        << " heuristic_top_k_applied_labels=" << log.heuristic_top_k_applied_label_count
        << " heuristic_top_k_edges_before=" << log.heuristic_top_k_feasible_edges_before
        << " heuristic_top_k_edges_after=" << log.heuristic_top_k_feasible_edges_after
        << " column_pool_enabled=" << (log.column_pool_enabled ? 1 : 0)
        << " column_pool_size_before_reprice=" << log.column_pool_size_before_reprice
        << " column_pool_max_reprice=" << log.column_pool_max_reprice
        << " column_pool_found_cols=" << log.column_pool_found_column_count
        << " column_pool_deferred_added=" << log.column_pool_deferred_added_count
        << " column_pool_size_after_update=" << log.column_pool_size_after_update
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
    FullEnumerationStaticPool full_enumeration_static_pool;
    FullEnumerationPoolNodeState full_enumeration_node_state;
    const bool use_full_enumeration_pool =
        phase == CGPhase::PhaseII &&
        options.phase_two_pricing_mode == NodeCGPhaseTwoPricingMode::FullEnumeration;

    if (use_full_enumeration_pool) {
        const auto pool_build_start_time = std::chrono::steady_clock::now();
        FullEnumerationPoolBuildOptions pool_options;
        pool_options.p = options.p;
        pool_options.time_limit = options.time_limit;
        pool_options.prune_pickup_symmetry_43 = options.prune_pickup_symmetry_43;
        pool_options.prune_delivery_symmetry_43 = options.prune_delivery_symmetry_43;
        FullEnumerationPoolBuildStats static_pool_stats;
        full_enumeration_static_pool = build_full_enumeration_static_pool(
            graph,
            pool_options,
            &static_pool_stats
        );
        FullEnumerationPoolBuildStats node_state_stats = static_pool_stats;
        full_enumeration_node_state = build_full_enumeration_pool_node_state(
            full_enumeration_static_pool,
            master_problem.column_id_by_key,
            &node_state_stats
        );
        if (log_stream != nullptr) {
            const double total_pool_build_runtime =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - pool_build_start_time
                )
                    .count();
            *log_stream << "[node-cg] phase2 full-enumeration pool raw_paths="
                        << static_pool_stats.raw_path_count
                        << " compact_psteps=" << static_pool_stats.compact_pstep_count
                        << " pool_size=" << node_state_stats.inactive_column_count
                        << " inactive_paths=" << node_state_stats.inactive_path_count
                        << " skipped_master=" << node_state_stats.skipped_master_column_count
                        << " skipped_duplicate=" << static_pool_stats.skipped_duplicate_column_count
                        << " build_runtime=" << format_double(total_pool_build_runtime)
                        << '\n';
        }
    }

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
        iteration_log.heuristic_engine = options.phase_two_heuristic_engine;
        iteration_log.heuristic_available_start_count =
            pricing_context.start_node_state_indices.size();
        iteration_log.heuristic_max_starts = options.phase_two_heuristic_max_starts;
        iteration_log.heuristic_start_ratio = options.phase_two_heuristic_start_ratio;
        iteration_log.heuristic_ladder_levels = options.phase_two_heuristic_ladder_levels;
        iteration_log.heuristic_output_max_columns_per_start =
            options.phase_two_heuristic_max_columns_per_start;
        iteration_log.heuristic_output_max_columns_total =
            options.phase_two_heuristic_max_total_columns;
        iteration_log.heuristic_search_column_ratio =
            options.phase_two_heuristic_search_column_ratio;
        iteration_log.heuristic_effective_search_max_columns_per_start =
            compute_search_column_limit(
                options.phase_two_heuristic_max_columns_per_start,
                options.phase_two_heuristic_search_column_ratio
            );
        iteration_log.heuristic_effective_search_max_columns_total =
            compute_search_column_limit(
                options.phase_two_heuristic_max_total_columns,
                options.phase_two_heuristic_search_column_ratio
            );
        iteration_log.heuristic_labeling_top_k_next =
            options.phase_two_labeling_top_k_next;
        iteration_log.heuristic_shallow_k1 = options.phase_two_shallow_k1;
        iteration_log.heuristic_shallow_k2 = options.phase_two_shallow_k2;
        iteration_log.column_pool_enabled =
            options.phase_two_column_pool_enabled || use_full_enumeration_pool;
        iteration_log.column_pool_max_reprice =
            use_full_enumeration_pool
                ? full_enumeration_node_state.available_variant_count
                : options.phase_two_column_pool_max_reprice;

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
        if (use_full_enumeration_pool) {
            iteration_log.column_pool_size_before_reprice =
                full_enumeration_node_state.available_variant_count;
            pricing_result = run_full_enumeration_pool_pricing(
                full_enumeration_static_pool,
                full_enumeration_node_state,
                dual_solution,
                phase,
                options.reduced_cost_tolerance,
                options.exact_pricing_max_total_columns_per_round,
                options.full_enumeration_rc_update_mode,
                options.full_enumeration_parallel_stage1_backend,
                options.full_enumeration_parallel_stage2_backend,
                options.full_enumeration_rc_update_threads,
                options.full_enumeration_rc_detail_log,
                log_stream
            );
            iteration_log.column_pool_found_column_count = pricing_result.columns.size();
        } else if (phase == CGPhase::PhaseII &&
            options.phase_two_pricing_mode ==
                NodeCGPhaseTwoPricingMode::HeuristicPricingThenExact) {
            iteration_log.heuristic_pricing_attempted = true;
            ForwardPricingResult heuristic_stats;
            heuristic_stats.best_reduced_cost = std::numeric_limits<double>::infinity();
            ForwardPricingResult selected_heuristic_result;
            bool heuristic_found = false;

            if (options.phase_two_column_pool_enabled) {
                iteration_log.column_pool_size_before_reprice = phase_two_column_pool.size();
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
                iteration_log.column_pool_found_column_count = pool_result.columns.size();
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
                    heuristic_stats.total_negative_column_count +=
                        level_result.total_negative_column_count;
                    heuristic_stats.per_start_search_cap_hit_count +=
                        level_result.per_start_search_cap_hit_count;
                    heuristic_stats.hit_global_search_cap =
                        heuristic_stats.hit_global_search_cap ||
                        level_result.hit_global_search_cap;
                    heuristic_stats.top_k_next_applied_label_count +=
                        level_result.top_k_next_applied_label_count;
                    heuristic_stats.top_k_next_feasible_edges_before +=
                        level_result.top_k_next_feasible_edges_before;
                    heuristic_stats.top_k_next_feasible_edges_after +=
                        level_result.top_k_next_feasible_edges_after;

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
            iteration_log.heuristic_total_negative_column_count =
                heuristic_stats.total_negative_column_count;
            iteration_log.heuristic_per_start_search_cap_hit_count =
                heuristic_stats.per_start_search_cap_hit_count;
            iteration_log.heuristic_global_search_cap_hit =
                heuristic_stats.hit_global_search_cap;
            iteration_log.heuristic_top_k_applied_label_count =
                heuristic_stats.top_k_next_applied_label_count;
            iteration_log.heuristic_top_k_feasible_edges_before =
                heuristic_stats.top_k_next_feasible_edges_before;
            iteration_log.heuristic_top_k_feasible_edges_after =
                heuristic_stats.top_k_next_feasible_edges_after;

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
        if (use_full_enumeration_pool && added_columns > 0U) {
            mark_full_enumeration_pool_columns_active(
                full_enumeration_static_pool,
                full_enumeration_node_state,
                pricing_result.columns
            );
        }
        if (phase == CGPhase::PhaseII &&
            options.phase_two_column_pool_enabled &&
            !use_full_enumeration_pool &&
            !columns_to_add_to_pool.empty()) {
            iteration_log.column_pool_deferred_added_count = columns_to_add_to_pool.size();
            add_columns_to_phase_two_pool(
                columns_to_add_to_pool,
                master_problem,
                options.phase_two_column_pool_max_size,
                phase_two_column_pool,
                phase_two_pool_index_by_key
            );
        }
        if (phase == CGPhase::PhaseII &&
            (options.phase_two_column_pool_enabled || use_full_enumeration_pool)) {
            iteration_log.column_pool_size_after_update =
                use_full_enumeration_pool
                    ? full_enumeration_node_state.available_variant_count
                    : phase_two_column_pool.size();
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
    std::ostream* log_stream,
    NodeCGPhaseOneMode heuristic_mode
) {
    if (log_stream != nullptr) {
        *log_stream << "[node-cg] Phase I mode=" << to_string(heuristic_mode) << '\n';
    }

    const auto heuristic_start_time = std::chrono::steady_clock::now();
    const std::vector<CGColumn> heuristic_columns =
        heuristic_mode == NodeCGPhaseOneMode::HeuristicCG3Step
            ? build_phase_one_heuristic_columns_3_step(
                  data,
                  graph,
                  pricing_context,
                  options
              )
            : build_phase_one_heuristic_columns(
                  data,
                  graph,
                  pricing_context,
                  options,
                  log_stream
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
                log_stream,
                NodeCGPhaseOneMode::HeuristicCG
            );
        case NodeCGPhaseOneMode::HeuristicCG3Step:
            return run_node_phase_one_heuristic_cg(
                data,
                graph,
                pricing_context,
                master_problem,
                options,
                global_start_time,
                result,
                log_stream,
                NodeCGPhaseOneMode::HeuristicCG3Step
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
        case NodeCGPhaseOneMode::HeuristicCG3Step:
            return "heuristic-cg-3-step";
    }
    return "unknown";
}

const char* to_string(NodeCGPhaseTwoPricingMode mode) {
    switch (mode) {
        case NodeCGPhaseTwoPricingMode::ExactPricing:
            return "exact-pricing";
        case NodeCGPhaseTwoPricingMode::HeuristicPricingThenExact:
            return "heuristic-pricing-then-exact";
        case NodeCGPhaseTwoPricingMode::FullEnumeration:
            return "full-enumeration";
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
