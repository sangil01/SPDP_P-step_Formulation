#include "GenMultiGraph.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <ostream>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace spdp {
namespace {

constexpr StateToken N_TOKEN{'N', -1, -1};

std::size_t hash_combine(std::size_t seed, std::size_t value) {
    seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
    return seed;
}

struct StateHash {
    std::size_t operator()(const State& state) const {
        std::size_t seed = 0;
        for (const StateToken& token : state) {
            seed = hash_combine(seed, static_cast<std::size_t>(token.kind));
            seed = hash_combine(seed, std::hash<int>{}(token.container_type));
            seed = hash_combine(seed, std::hash<int>{}(token.landfill_location));
        }
        return seed;
    }
};

struct PairHash {
    std::size_t operator()(const std::pair<int, int>& value) const {
        std::size_t seed = std::hash<int>{}(value.first);
        seed = hash_combine(seed, std::hash<int>{}(value.second));
        return seed;
    }
};

struct EdgeBucketKey {
    NodeId u = 0;
    NodeId v = 0;
    State start_state{};
    State end_state{};

    bool operator==(const EdgeBucketKey& other) const {
        return u == other.u && v == other.v && start_state == other.start_state &&
               end_state == other.end_state;
    }
};

struct EdgeBucketKeyHash {
    std::size_t operator()(const EdgeBucketKey& key) const {
        StateHash state_hash;
        std::size_t seed = std::hash<int>{}(key.u);
        seed = hash_combine(seed, std::hash<int>{}(key.v));
        seed = hash_combine(seed, state_hash(key.start_state));
        seed = hash_combine(seed, state_hash(key.end_state));
        return seed;
    }
};

struct EdgeBucketEntry {
    EdgeBucketKey key;
    std::vector<EdgeData> candidates;
};

bool token_less(const StateToken& lhs, const StateToken& rhs) {
    return lhs < rhs;
}

bool state_less(const State& lhs, const State& rhs) {
    if (token_less(lhs[0], rhs[0])) {
        return true;
    }
    if (token_less(rhs[0], lhs[0])) {
        return false;
    }
    return token_less(lhs[1], rhs[1]);
}

std::string token_to_str(const StateToken& token) {
    return "('" + std::string(1, token.kind) + "', " + std::to_string(token.container_type) +
           ", " + std::to_string(token.landfill_location) + ")";
}

StateToken make_e(int container_type) {
    return StateToken{'E', container_type, -1};
}

StateToken make_f(int container_type, int landfill_location) {
    return StateToken{'F', container_type, landfill_location};
}

//container state가 순서에 상관 없으므로, 항상 일관된 순서로 상태를 가지게 함으로써, 비교, 해시, 중복 제거 등이 일관되게 함
State canonical_state(const StateToken& first, const StateToken& second) {
    State state{first, second};
    if (state[1] < state[0]) {
        std::swap(state[0], state[1]);
    }
    return state;
}

State canonical_state(const State& state) {
    return canonical_state(state[0], state[1]);
}

void log_line(std::ostream& out, const std::string& message) {
    out << message << '\n';
}

void log_graph_generation_summary(
    std::ostream& out,
    const std::vector<State>& all_states,
    int virtual_location,
    std::size_t initial_edge_count,
    std::size_t removed_infeasible_edge_count,
    std::size_t removed_dominated_edge_count,
    std::size_t removed_symmetry_40_edge_count,
    std::size_t removed_symmetry_41_edge_count,
    std::size_t kept_edge_count,
    const GraphBuildOptions& options
) {
    log_line(out, "[GenMultiGraph] State count: " + std::to_string(all_states.size()));
    log_line(out, "[GenMultiGraph] Initial edge count: " + std::to_string(initial_edge_count));

    if (options.prune_infeasible_edges) {
        log_line(
            out,
            "[GenMultiGraph] Infeasible-pruned edge count: " +
                std::to_string(removed_infeasible_edge_count)
        );
    } else {
        log_line(out, "[GenMultiGraph] Infeasible pruning disabled");
    }

    if (options.prune_dominated_edges) {
        log_line(
            out,
            "[GenMultiGraph] Dominated-pruned edge count: " +
                std::to_string(removed_dominated_edge_count)
        );
    } else {
        log_line(out, "[GenMultiGraph] Dominated pruning disabled");
    }

    // [수정] 논문의 symmetry breaking 제약 (40), (41)에 대응하는 edge pruning 통계를 로그에 남긴다.
    if (options.prune_symmetry_40) {
        log_line(
            out,
            "[GenMultiGraph] Symmetry-40-pruned edge count: " +
                std::to_string(removed_symmetry_40_edge_count)
        );
    } else {
        log_line(out, "[GenMultiGraph] Symmetry-40 pruning disabled");
    }

    // [수정] 같은 방식으로 논문의 symmetry breaking 제약 (41) pruning 통계를 기록한다.
    if (options.prune_symmetry_41) {
        log_line(
            out,
            "[GenMultiGraph] Symmetry-41-pruned edge count: " +
                std::to_string(removed_symmetry_41_edge_count)
        );
    } else {
        log_line(out, "[GenMultiGraph] Symmetry-41 pruning disabled");
    }

    log_line(out, "[GenMultiGraph] Final edge count: " + std::to_string(kept_edge_count));
}

double service_time(const NodeSpec& node, const SPDPData& data) {
    if (node.kind == NodeSpec::Kind::Pickup) {
        return data.time_pickup;
    }
    if (node.kind == NodeSpec::Kind::Delivery) {
        return data.time_delivery;
    }
    return 0.0;
}

double path_metric(
    const std::vector<std::vector<double>>& metric_matrix,
    int start_loc,
    const std::vector<int>& sequence_pi,
    int end_loc
) {
    double value = 0.0;
    int prev = start_loc;

    for (int location : sequence_pi) {
        value += metric_matrix[static_cast<std::size_t>(prev)][static_cast<std::size_t>(location)];
        prev = location;
    }

    value += metric_matrix[static_cast<std::size_t>(prev)][static_cast<std::size_t>(end_loc)];
    return value;
}

bool dominates(const EdgeData& a, const EdgeData& b) {
    const bool weak = (a.time <= b.time && a.cost <= b.cost);
    const bool strict = (a.time < b.time || a.cost < b.cost);
    return weak && strict;
}

bool has_same_objective_values(const EdgeData& a, const EdgeData& b) {
    return a.time == b.time && a.cost == b.cost;
}

std::size_t find_or_create_edge_bucket(
    const EdgeBucketKey& key,
    std::vector<EdgeBucketEntry>& edge_bucket_entries,
    std::unordered_map<EdgeBucketKey, std::size_t, EdgeBucketKeyHash>& edge_bucket_indices
) {
    const auto found = edge_bucket_indices.find(key);
    if (found != edge_bucket_indices.end()) {
        return found->second;
    }

    const std::size_t index = edge_bucket_entries.size();
    edge_bucket_indices.emplace(key, index);

    EdgeBucketEntry entry;
    entry.key = key;
    edge_bucket_entries.push_back(std::move(entry));
    return index;
}

std::size_t insert_into_pareto_bucket(
    std::vector<EdgeData>& candidates,
    EdgeData edge_data,
    bool prune_dominated_edges
) {
    for (const EdgeData& candidate : candidates) {
        if (has_same_objective_values(candidate, edge_data)) {
            return 1;
        }
    }

    if (!prune_dominated_edges) {
        candidates.push_back(std::move(edge_data));
        return 0;
    }

    for (const EdgeData& candidate : candidates) {
        if (dominates(candidate, edge_data)) {
            return 1;
        }
    }

    const std::size_t previous_size = candidates.size();
    candidates.erase(
        std::remove_if(
            candidates.begin(),
            candidates.end(),
            [&edge_data](const EdgeData& candidate) {
                return dominates(edge_data, candidate);
            }
        ),
        candidates.end()
    );
    const std::size_t removed_existing_count = previous_size - candidates.size();

    candidates.push_back(std::move(edge_data));
    return removed_existing_count;
}

std::optional<State> reconstruct_pre_service_state(const NodeSpec& node, const State& after_state) {
    if (node.kind != NodeSpec::Kind::Delivery) {
        return std::nullopt;
    }

    State before_state = after_state;
    const StateToken delivered_empty = make_e(node.container_type.value());

    for (StateToken& token : before_state) {
        if (token.kind == 'N') {
            token = delivered_empty;
            return canonical_state(before_state);
        }
    }

    throw std::runtime_error("Delivery after-state has no feasible pre-service state.");
}

bool violates_node_state_rules(
    const NodeSpec& node,
    const State& after_state,
    const std::unordered_map<State, bool, StateHash>& violates_cache
) {
    if (violates_cache.at(after_state)) {
        return true;
    }

    const std::optional<State> before_state = reconstruct_pre_service_state(node, after_state);
    if (!before_state.has_value()) {
        return false;
    }

    return violates_cache.at(before_state.value());
}

bool violates_singleton_delivery_to_own_pickup_rule(
    const NodeSpec& u,
    const NodeSpec& v,
    const std::unordered_set<int>& singleton_types
) {
    if (u.kind != NodeSpec::Kind::Delivery || v.kind != NodeSpec::Kind::Pickup) {
        return false;
    }
    if (!u.request_idx.has_value() || !v.request_idx.has_value()) {
        return false;
    }
    if (u.request_idx.value() != v.request_idx.value()) {
        return false;
    }
    if (!u.container_type.has_value()) {
        return false;
    }

    return singleton_types.find(u.container_type.value()) != singleton_types.end();
}

// [수정] 논문 식 (40a)의 compressed-edge 대응: 같은 location의 연속 pickup은 index 오름차순만 허용한다.
bool violates_symmetry_40a_rule(
    const NodeSpec& u,
    const NodeSpec& v,
    const std::vector<int>& sequence_pi
) {
    return sequence_pi.empty() && u.kind == NodeSpec::Kind::Pickup &&
           v.kind == NodeSpec::Kind::Pickup && u.location == v.location &&
           u.request_idx.has_value() && v.request_idx.has_value() &&
           u.request_idx.value() > v.request_idx.value();
}

// [수정] 논문 식 (40c)의 compressed-edge 대응: 같은 location의 연속 delivery도 index 오름차순만 허용한다.
bool violates_symmetry_40c_rule(
    const NodeSpec& u,
    const NodeSpec& v,
    const std::vector<int>& sequence_pi
) {
    return sequence_pi.empty() && u.kind == NodeSpec::Kind::Delivery &&
           v.kind == NodeSpec::Kind::Delivery && u.location == v.location &&
           u.request_idx.has_value() && v.request_idx.has_value() &&
           u.request_idx.value() > v.request_idx.value();
}

// [수정] 현재 compressed 모델에서 식 (40)에 해당하는 local 대칭성 제거 여부를 모은다.
bool violates_symmetry_40_rules(
    const NodeSpec& u,
    const NodeSpec& v,
    const std::vector<int>& sequence_pi
) {
    return violates_symmetry_40a_rule(u, v, sequence_pi) ||
           violates_symmetry_40c_rule(u, v, sequence_pi);
}

// [수정] 논문 식 (41a)의 compressed-edge 대응: 같은 location의 연속 pickup -> delivery 순서를 제거한다.
bool violates_symmetry_41a_rule(
    const NodeSpec& u,
    const NodeSpec& v,
    const std::vector<int>& sequence_pi
) {
    return sequence_pi.empty() && u.kind == NodeSpec::Kind::Pickup &&
           v.kind == NodeSpec::Kind::Delivery && u.location == v.location;
}

// [수정] 논문 식 (41b)의 compressed-edge 대응: 마지막 treatment와 도착 pickup이 같은 location이면 제거한다.
bool violates_symmetry_41b_rule(
    const NodeSpec& v,
    const std::vector<int>& sequence_pi
) {
    return !sequence_pi.empty() && v.kind == NodeSpec::Kind::Pickup &&
           sequence_pi.back() == v.location;
}

// [수정] 논문 식 (41c)의 compressed-edge 대응: 출발 delivery 직후 첫 treatment가 같은 location이면 제거한다.
bool violates_symmetry_41c_rule(
    const NodeSpec& u,
    const std::vector<int>& sequence_pi
) {
    return !sequence_pi.empty() && u.kind == NodeSpec::Kind::Delivery &&
           sequence_pi.front() == u.location;
}

// [수정] 현재 compressed 모델에서 식 (41)에 해당하는 local 대칭성 제거 여부를 모은다.
bool violates_symmetry_41_rules(
    const NodeSpec& u,
    const NodeSpec& v,
    const std::vector<int>& sequence_pi
) {
    return violates_symmetry_41a_rule(u, v, sequence_pi) ||
           violates_symmetry_41b_rule(v, sequence_pi) ||
           violates_symmetry_41c_rule(u, sequence_pi);
}

bool try_add_edge_candidate(
    const SPDPData& data,
    const NodeSpec& u,
    const NodeSpec& v,
    const State& sigma_u,
    bool sigma_u_violates,
    const std::optional<State>& sigma_v,
    const std::vector<int>& sequence_pi,
    int emptied_count,
    const std::unordered_map<State, bool, StateHash>& violates_cache,
    const std::unordered_set<int>& singleton_types,
    const GraphBuildOptions& options,
    std::vector<EdgeBucketEntry>& edge_bucket_entries,
    std::unordered_map<EdgeBucketKey, std::size_t, EdgeBucketKeyHash>& edge_bucket_indices,
    std::size_t& initial_edge_count,
    std::size_t& removed_infeasible_edge_count,
    std::size_t& removed_symmetry_40_edge_count,
    std::size_t& removed_symmetry_41_edge_count,
    std::size_t& removed_dominated_edge_count
) {
    if (!sigma_v.has_value()) {
        return false;
    }

    const double travel_time = path_metric(data.time, u.location, sequence_pi, v.location);
    double travel_cost = path_metric(data.distance, u.location, sequence_pi, v.location);

    if (u.node_id == 0 && v.kind != NodeSpec::Kind::End) {
        travel_cost += data.fixed_vehicle_cost;
    }

    const double total_time =
        travel_time + static_cast<double>(emptied_count) * data.time_empty + service_time(v, data);

    ++initial_edge_count;

    const bool sigma_v_violates = violates_node_state_rules(v, sigma_v.value(), violates_cache);
    const bool violates_state_rules = sigma_u_violates || sigma_v_violates;
    const bool violates_time_limit = total_time > data.time_limit;
    const bool violates_singleton_delivery_to_pickup =
        violates_singleton_delivery_to_own_pickup_rule(u, v, singleton_types);
    const bool violates_symmetry_40 =
        options.prune_symmetry_40 && violates_symmetry_40_rules(u, v, sequence_pi);
    const bool violates_symmetry_41 =
        options.prune_symmetry_41 && violates_symmetry_41_rules(u, v, sequence_pi);

    if (violates_symmetry_40) {
        ++removed_symmetry_40_edge_count;
        return false;
    }

    if (violates_symmetry_41) {
        ++removed_symmetry_41_edge_count;
        return false;
    }

    if (options.prune_infeasible_edges &&
        (violates_state_rules || violates_time_limit || violates_singleton_delivery_to_pickup)) {
        ++removed_infeasible_edge_count;
        return false;
    }

    EdgeData edge_data;
    edge_data.sequence_pi = sequence_pi;
    edge_data.time = total_time;
    edge_data.cost = travel_cost;
    edge_data.start_state = sigma_u;
    edge_data.end_state = sigma_v.value();

    EdgeBucketKey key;
    key.u = u.node_id;
    key.v = v.node_id;
    key.start_state = sigma_u;
    key.end_state = sigma_v.value();

    const std::size_t bucket_index =
        find_or_create_edge_bucket(key, edge_bucket_entries, edge_bucket_indices);
    std::vector<EdgeData>& candidates = edge_bucket_entries[bucket_index].candidates;

    const std::size_t previous_size = candidates.size();
    const std::size_t dominated_removed_count =
        insert_into_pareto_bucket(candidates, std::move(edge_data), options.prune_dominated_edges);

    if (candidates.size() == previous_size) {
        removed_dominated_edge_count += dominated_removed_count;
        return false;
    }

    removed_dominated_edge_count += dominated_removed_count;
    return true;
}

std::vector<NodeSpec> build_node_specs(const SPDPData& data) {
    const int n_requests = static_cast<int>(data.requests.size());
    const int end_depot_id = 2 * n_requests + 1;
    const int virtual_location = data.locations;

    std::vector<NodeSpec> nodes;
    nodes.reserve(static_cast<std::size_t>(2 * n_requests + 2));

    nodes.push_back(NodeSpec{
        0,
        NodeSpec::Kind::Start,
        virtual_location,
        std::nullopt,
        std::nullopt,
        std::nullopt,
    });

    nodes.push_back(NodeSpec{
        end_depot_id,
        NodeSpec::Kind::End,
        virtual_location,
        std::nullopt,
        std::nullopt,
        std::nullopt,
    });

    for (int idx = 0; idx < n_requests; ++idx) {
        const Request& req = data.requests[static_cast<std::size_t>(idx)];
        const int pickup_id = idx + 1;
        const int delivery_id = n_requests + idx + 1;

        nodes.push_back(NodeSpec{
            pickup_id,
            NodeSpec::Kind::Pickup,
            req.from_id,
            idx,
            req.container_type,
            req.to_id,
        });

        nodes.push_back(NodeSpec{
            delivery_id,
            NodeSpec::Kind::Delivery,
            req.from_id,
            idx,
            req.container_type,
            req.to_id,
        });
    }

    return nodes;
}

std::vector<State> generate_state_space(const SPDPData& data) {
    std::set<int> container_type_set;
    std::set<std::pair<int, int>> full_item_set;

    for (const Request& req : data.requests) {
        container_type_set.insert(req.container_type);
        full_item_set.insert({req.container_type, req.to_id});
    }

    std::vector<StateToken> atom_space;
    atom_space.reserve(1 + container_type_set.size() + full_item_set.size());
    atom_space.push_back(N_TOKEN);

    for (int container_type : container_type_set) {
        atom_space.push_back(make_e(container_type));
    }

    for (const auto& item : full_item_set) {
        atom_space.push_back(make_f(item.first, item.second));
    }

    std::vector<State> states;
    states.reserve((atom_space.size() * (atom_space.size() + 1U)) / 2U);

    for (std::size_t i = 0; i < atom_space.size(); ++i) {
        for (std::size_t j = i; j < atom_space.size(); ++j) {
            states.push_back(canonical_state(atom_space[i], atom_space[j]));
        }
    }

    return states;
}

std::optional<State> apply_service(const NodeSpec& node, const State& state) {
    State tokens = state;
    const State empty_state = canonical_state(N_TOKEN, N_TOKEN);

    if (node.kind == NodeSpec::Kind::Start || node.kind == NodeSpec::Kind::End) {
        const State candidate = canonical_state(tokens);
        if (candidate == empty_state) {
            return candidate;
        }
        return std::nullopt;
    }

    if (node.kind == NodeSpec::Kind::Pickup) {
        for (StateToken& token : tokens) {
            if (token.kind == 'N') {
                token = make_f(node.container_type.value(), node.landfill_location.value());
                return canonical_state(tokens);
            }
        }
        return std::nullopt;
    }

    if (node.kind == NodeSpec::Kind::Delivery) {
        const StateToken target = make_e(node.container_type.value());
        for (StateToken& token : tokens) {
            if (token == target) {
                token = N_TOKEN;
                return canonical_state(tokens);
            }
        }
        return std::nullopt;
    }

    throw std::runtime_error("Unknown node kind in apply_service.");
}

std::vector<std::vector<int>> candidate_sequences(const State& state) {
    std::vector<int> full_locations;
    full_locations.reserve(2);

    for (const StateToken& token : state) {
        if (token.kind == 'F') {
            full_locations.push_back(token.landfill_location);
        }
    }

    if (full_locations.empty()) {
        return {std::vector<int>{}};
    }

    if (full_locations.size() == 1U) {
        return {std::vector<int>{}, std::vector<int>{full_locations[0]}};
    }

    if (full_locations.size() == 2U) {
        const int l1 = full_locations[0];
        const int l2 = full_locations[1];

        if (l1 == l2) {
            return {std::vector<int>{}, std::vector<int>{l1}};
        }

        std::vector<std::vector<int>> seq = {
            std::vector<int>{},
            std::vector<int>{l1},
            std::vector<int>{l2},
            std::vector<int>{l1, l2},
            std::vector<int>{l2, l1},
        };
        return seq;
    }

    throw std::runtime_error("Capacity is 2, so full skip count must be <= 2");
}

std::pair<State, int> apply_emptying(const State& state, const std::vector<int>& sequence_pi) {
    State tokens = state;
    int emptied_count = 0;

    for (int treatment_location : sequence_pi) {
        for (StateToken& token : tokens) {
            if (token.kind == 'F' && token.landfill_location == treatment_location) {
                token = make_e(token.container_type);
                ++emptied_count;
            }
        }
    }

    return {canonical_state(tokens), emptied_count};
}

bool violates_single_request_state_rules(
    const State& state,
    const std::unordered_set<int>& singleton_types,
    const std::unordered_set<std::pair<int, int>, PairHash>& singleton_full_pairs
) {
    std::unordered_map<int, int> type_count;

    for (const StateToken& token : state) {
        if (token.kind == 'E' || token.kind == 'F') {
            ++type_count[token.container_type];
        }
    }

    for (const auto& entry : type_count) {
        if (entry.second >= 2 && singleton_types.find(entry.first) != singleton_types.end()) {
            return true;
        }
    }

    if (state[0].kind == 'F' && state[1].kind == 'F') {
        if (state[0].container_type == state[1].container_type &&
            state[0].landfill_location == state[1].landfill_location) {
            std::pair<int, int> pair_key{state[0].container_type, state[0].landfill_location};
            if (singleton_full_pairs.find(pair_key) != singleton_full_pairs.end()) {
                return true;
            }
        }
    }

    return false;
}

bool is_weakly_connected(
    const std::vector<NodeSpec>& node_specs,
    const std::vector<EdgeRecord>& edges
) {
    if (node_specs.empty()) {
        return true;
    }

    std::unordered_map<NodeId, std::vector<NodeId>> adjacency;
    adjacency.reserve(node_specs.size());
    for (const NodeSpec& node : node_specs) {
        adjacency[node.node_id] = {};
    }

    for (const EdgeRecord& edge : edges) {
        adjacency[edge.u].push_back(edge.v);
        adjacency[edge.v].push_back(edge.u);
    }

    std::unordered_set<NodeId> visited;
    visited.reserve(node_specs.size());
    std::vector<NodeId> stack;
    stack.reserve(node_specs.size());

    const NodeId start = node_specs.front().node_id;
    stack.push_back(start);
    visited.insert(start);

    while (!stack.empty()) {
        const NodeId current = stack.back();
        stack.pop_back();

        const auto found = adjacency.find(current);
        if (found == adjacency.end()) {
            continue;
        }
        for (NodeId next : found->second) {
            if (visited.insert(next).second) {
                stack.push_back(next);
            }
        }
    }

    return visited.size() == node_specs.size();
}

}  // namespace

bool CanonicalEdgeKey::operator==(const CanonicalEdgeKey& other) const {
    return u == other.u && v == other.v && sequence_pi == other.sequence_pi &&
           time == other.time && cost == other.cost &&
           start_state == other.start_state && end_state == other.end_state;
}

bool CanonicalEdgeKey::operator<(const CanonicalEdgeKey& other) const {
    if (u != other.u) {
        return u < other.u;
    }
    if (v != other.v) {
        return v < other.v;
    }
    if (sequence_pi != other.sequence_pi) {
        return sequence_pi < other.sequence_pi;
    }
    if (time != other.time) {
        return time < other.time;
    }
    if (cost != other.cost) {
        return cost < other.cost;
    }
    if (start_state != other.start_state) {
        return state_less(start_state, other.start_state);
    }
    if (end_state != other.end_state) {
        return state_less(end_state, other.end_state);
    }
    return false;
}

void MultiDiGraph::add_node(const NodeSpec& node) {
    nodes_[node.node_id] = node;
}

void MultiDiGraph::add_edge(NodeId u, NodeId v, const EdgeData& edge_data) {
    const std::pair<NodeId, NodeId> pair_key{u, v};
    const int key = next_key_by_pair_[pair_key]++;

    EdgeRecord record;
    record.u = u;
    record.v = v;
    record.key = key;
    record.data = edge_data;

    const std::size_t edge_index = edges_.size();
    edges_.push_back(std::move(record));
    outgoing_edge_indices_[u].push_back(edge_index);
    ingoing_edge_indices_[v].push_back(edge_index);
}

std::size_t MultiDiGraph::number_of_nodes() const {
    return nodes_.size();
}

std::size_t MultiDiGraph::number_of_edges() const {
    return edges_.size();
}

NodeId MultiDiGraph::end_node_id() const {
    NodeId end_node_id = 0;
    for (const auto& entry : nodes_) {
        end_node_id = std::max(end_node_id, entry.first);
    }
    return end_node_id;
}

bool MultiDiGraph::is_physical_service_node(NodeId node_id) const {
    return node_id != 0 && node_id != end_node_id();
}

const NodeSpec& MultiDiGraph::node(NodeId node_id) const {
    const auto found = nodes_.find(node_id);
    if (found == nodes_.end()) {
        throw std::runtime_error("Requested node does not exist in MultiDiGraph.");
    }
    return found->second;
}

const std::vector<EdgeRecord>& MultiDiGraph::edges() const {
    return edges_;
}

std::vector<NodeSpec> MultiDiGraph::nodes() const {
    std::vector<NodeSpec> result;
    result.reserve(nodes_.size());
    for (const auto& entry : nodes_) {
        result.push_back(entry.second);
    }
    std::sort(
        result.begin(),
        result.end(),
        [](const NodeSpec& left, const NodeSpec& right) {
            return left.node_id < right.node_id;
        }
    );
    return result;
}

const std::vector<std::size_t>& MultiDiGraph::outgoing_edge_indices(NodeId u) const {
    static const std::vector<std::size_t> kEmptyIndices;

    const auto found = outgoing_edge_indices_.find(u);
    if (found == outgoing_edge_indices_.end()) {
        return kEmptyIndices;
    }
    return found->second;
}

const std::vector<std::size_t>& MultiDiGraph::ingoing_edge_indices(NodeId v) const {
    static const std::vector<std::size_t> kEmptyIndices;

    const auto found = ingoing_edge_indices_.find(v);
    if (found == ingoing_edge_indices_.end()) {
        return kEmptyIndices;
    }
    return found->second;
}

std::vector<EdgeRecord> MultiDiGraph::edges_from(NodeId u) const {
    std::vector<EdgeRecord> result;

    const std::vector<std::size_t>& indices = outgoing_edge_indices(u);
    result.reserve(indices.size());
    for (std::size_t index : indices) {
        result.push_back(edges_[index]);
    }

    return result;
}

std::vector<EdgeRecord> MultiDiGraph::edges_to(NodeId v) const {
    std::vector<EdgeRecord> result;

    const std::vector<std::size_t>& indices = ingoing_edge_indices(v);
    result.reserve(indices.size());
    for (std::size_t index : indices) {
        result.push_back(edges_[index]);
    }

    return result;
}

std::vector<EdgeRecord> MultiDiGraph::edges_between(NodeId u, NodeId v) const {
    std::vector<EdgeRecord> result;

    for (std::size_t index : outgoing_edge_indices(u)) {
        const EdgeRecord& record = edges_[index];
        if (record.v == v) {
            result.push_back(record);
        }
    }

    return result;
}

std::size_t MultiDiGraph::PairHash::operator()(const std::pair<NodeId, NodeId>& value) const {
    std::size_t seed = std::hash<int>{}(value.first);
    seed = hash_combine(seed, std::hash<int>{}(value.second));
    return seed;
}

std::string state_to_str(const State& state) {
    return token_to_str(state[0]) + "|" + token_to_str(state[1]);
}

MultiDiGraph build_multigraph(
    const SPDPData& data,
    const GraphBuildOptions& options,
    std::ostream* log_stream
) {
    std::ostream& log_out = log_stream != nullptr ? *log_stream : std::cout;

    MultiDiGraph graph;
    const std::vector<NodeSpec> node_specs = build_node_specs(data);

    for (const NodeSpec& node : node_specs) {
        graph.add_node(node);
    }

    const int virtual_location = data.locations;
    
    const std::vector<State> all_states = generate_state_space(data);

    std::unordered_map<int, int> request_type_count;
    std::unordered_map<std::pair<int, int>, int, PairHash> request_full_pair_count;

    for (const Request& req : data.requests) {
        ++request_type_count[req.container_type];
        ++request_full_pair_count[{req.container_type, req.to_id}];
    }
    
    std::unordered_set<int> singleton_types;
    for (const auto& entry : request_type_count) {
        if (entry.second == 1) {
            singleton_types.insert(entry.first);
        }
    }

    std::unordered_set<std::pair<int, int>, PairHash> singleton_full_pairs;
    for (const auto& entry : request_full_pair_count) {
        if (entry.second == 1) {
            singleton_full_pairs.insert(entry.first);
        }
    }

    std::unordered_map<State, std::vector<std::vector<int>>, StateHash> sequence_cache;
    std::unordered_map<State, bool, StateHash> violates_cache;

    sequence_cache.reserve(all_states.size());
    violates_cache.reserve(all_states.size());

    for (const State& state : all_states) {
        sequence_cache.emplace(state, candidate_sequences(state));
        violates_cache.emplace(
            state,
            violates_single_request_state_rules(state, singleton_types, singleton_full_pairs)
        );
    }

    std::unordered_map<NodeId, std::vector<State>> feasible_after_states;
    feasible_after_states.reserve(node_specs.size());

    for (const NodeSpec& node : node_specs) {
        std::unordered_set<State, StateHash> states;
        states.reserve(all_states.size());

        for (const State& state : all_states) {
            const std::optional<State> after = apply_service(node, state);
            if (after.has_value()) {
                states.insert(after.value());
            }
        }

        std::vector<State> state_vector;
        state_vector.reserve(states.size());
        for (const State& state : states) {
            state_vector.push_back(state);
        }
        std::sort(state_vector.begin(), state_vector.end(), state_less);

        feasible_after_states[node.node_id] = std::move(state_vector);
    }

    std::vector<EdgeBucketEntry> edge_bucket_entries;
    std::unordered_map<EdgeBucketKey, std::size_t, EdgeBucketKeyHash> edge_bucket_indices;
    std::size_t initial_edge_count = 0;
    std::size_t removed_infeasible_edge_count = 0;
    std::size_t removed_dominated_edge_count = 0;
    std::size_t removed_symmetry_40_edge_count = 0;
    std::size_t removed_symmetry_41_edge_count = 0;

    edge_bucket_entries.reserve(4096);
    edge_bucket_indices.reserve(4096);

    //edge 기준은 start: u 서비스 이후, end: v 서비스 이후 (단, u가 virtual star인 경우는 u 서비스도 포함 (즉, fixed cost도 포함))
    for (const NodeSpec& u : node_specs) {
        if (u.kind == NodeSpec::Kind::End) {
            continue;
        }
        const auto& feasible_states_u = feasible_after_states.at(u.node_id);

        for (const NodeSpec& v : node_specs) {
            if (v.kind == NodeSpec::Kind::Start) {
                continue;
            }
            if (u.node_id == v.node_id) {
                continue;
            }

            for (const State& sigma_u : feasible_states_u) {
                const bool sigma_u_violates =
                    violates_node_state_rules(u, sigma_u, violates_cache);
                const std::vector<std::vector<int>>& sequences = sequence_cache.at(sigma_u);
                
                for (const std::vector<int>& sequence_pi : sequences) {
                    auto emptying_result = apply_emptying(sigma_u, sequence_pi);
                    const State& sigma_arrival = emptying_result.first;
                    const int emptied_count = emptying_result.second;

                    const std::optional<State> sigma_v = apply_service(v, sigma_arrival);
                    try_add_edge_candidate(
                        data,
                        u,
                        v,
                        sigma_u,
                        sigma_u_violates,
                        sigma_v,
                        sequence_pi,
                        emptied_count,
                        violates_cache,
                        singleton_types,
                        options,
                        edge_bucket_entries,
                        edge_bucket_indices,
                        initial_edge_count,
                        removed_infeasible_edge_count,
                        removed_symmetry_40_edge_count,
                        removed_symmetry_41_edge_count,
                        removed_dominated_edge_count
                    );
                }
            }
        }
    }

    std::size_t kept_edge_count = 0;

    for (EdgeBucketEntry& entry : edge_bucket_entries) {
        kept_edge_count += entry.candidates.size();

        for (EdgeData& edge_data : entry.candidates) {
            graph.add_edge(entry.key.u, entry.key.v, edge_data);
        }
    }

    log_graph_generation_summary(
        log_out,
        all_states,
        virtual_location,
        initial_edge_count,
        removed_infeasible_edge_count,
        removed_dominated_edge_count,
        removed_symmetry_40_edge_count,
        removed_symmetry_41_edge_count,
        kept_edge_count,
        options
    );

    if (!is_weakly_connected(node_specs, graph.edges())) {
        throw std::runtime_error("Generated multi-graph is disconnected.");
    }

    return graph;
}

namespace {

struct EndpointStateKey {
    NodeId u = 0;
    NodeId v = 0;
    State start_state{};
    State end_state{};

    bool operator<(const EndpointStateKey& other) const {
        if (u != other.u) {
            return u < other.u;
        }
        if (v != other.v) {
            return v < other.v;
        }
        if (start_state != other.start_state) {
            return start_state < other.start_state;
        }
        return end_state < other.end_state;
    }
};

bool is_empty_state(const State& state) {
    return std::all_of(
        state.begin(),
        state.end(),
        [](const StateToken& token) { return token.kind == 'N'; }
    );
}

bool is_empty_delivery_pickup_connector(
    const MultiDiGraph& graph,
    const EdgeRecord& edge
) {
    return graph.node(edge.u).kind == NodeSpec::Kind::Delivery &&
        graph.node(edge.v).kind == NodeSpec::Kind::Pickup &&
        is_empty_state(edge.data.start_state) &&
        edge.data.sequence_pi.empty();
}

bool better_duration_representative(
    const EdgeRecord& candidate,
    std::size_t candidate_id,
    const EdgeRecord& incumbent,
    std::size_t incumbent_id
) {
    if (candidate.data.time != incumbent.data.time) {
        return candidate.data.time < incumbent.data.time;
    }
    if (candidate.data.cost != incumbent.data.cost) {
        return candidate.data.cost < incumbent.data.cost;
    }
    if (candidate.data.sequence_pi != incumbent.data.sequence_pi) {
        return candidate.data.sequence_pi < incumbent.data.sequence_pi;
    }
    return candidate_id < incumbent_id;
}

void fnv_mix_bytes(
    std::uint64_t& hash,
    const unsigned char* bytes,
    std::size_t count
) {
    constexpr std::uint64_t kPrime = 1099511628211ULL;
    for (std::size_t i = 0; i < count; ++i) {
        hash ^= static_cast<std::uint64_t>(bytes[i]);
        hash *= kPrime;
    }
}

template <typename T>
void fnv_mix(std::uint64_t& hash, const T& value) {
    fnv_mix_bytes(
        hash,
        reinterpret_cast<const unsigned char*>(&value),
        sizeof(T)
    );
}

void fnv_mix_state(std::uint64_t& hash, const State& state) {
    for (const StateToken& token : state) {
        fnv_mix(hash, token.kind);
        fnv_mix(hash, token.container_type);
        fnv_mix(hash, token.landfill_location);
    }
}

}  // namespace

std::uint64_t multigraph_fingerprint(const MultiDiGraph& graph) {
    std::uint64_t hash = 1469598103934665603ULL;
    const std::vector<NodeSpec> nodes = graph.nodes();
    const std::size_t node_count = nodes.size();
    const std::size_t edge_count = graph.number_of_edges();
    fnv_mix(hash, node_count);
    fnv_mix(hash, edge_count);
    for (const NodeSpec& node : nodes) {
        fnv_mix(hash, node.node_id);
        const int kind = static_cast<int>(node.kind);
        fnv_mix(hash, kind);
        fnv_mix(hash, node.location);
        const int request = node.request_idx.value_or(-1);
        const int type = node.container_type.value_or(-1);
        const int landfill = node.landfill_location.value_or(-1);
        fnv_mix(hash, request);
        fnv_mix(hash, type);
        fnv_mix(hash, landfill);
    }
    for (const EdgeRecord& edge : graph.edges()) {
        fnv_mix(hash, edge.u);
        fnv_mix(hash, edge.v);
        fnv_mix(hash, edge.key);
        const std::size_t sequence_size = edge.data.sequence_pi.size();
        fnv_mix(hash, sequence_size);
        for (int location : edge.data.sequence_pi) {
            fnv_mix(hash, location);
        }
        fnv_mix(hash, edge.data.time);
        fnv_mix(hash, edge.data.cost);
        fnv_mix_state(hash, edge.data.start_state);
        fnv_mix_state(hash, edge.data.end_state);
    }
    return hash;
}

const char* graph_purpose_name(GraphPurpose purpose) {
    switch (purpose) {
        case GraphPurpose::Main:
            return "main";
        case GraphPurpose::InitialIncumbent:
            return "initial-incumbent";
        case GraphPurpose::DurationBound:
            return "duration-bound";
    }
    return "unknown";
}

DerivedMultiGraph derive_multigraph(
    const MultiDiGraph& main_graph,
    GraphPurpose purpose,
    const DerivedGraphOptions& options
) {
    if (purpose == GraphPurpose::InitialIncumbent &&
        options.prune_empty_state_connectors) {
        throw std::invalid_argument(
            "Empty-state connector pruning is unsafe for an exact-K "
            "initial-incumbent graph."
        );
    }
    const auto build_start = std::chrono::steady_clock::now();
    DerivedMultiGraph result;
    result.stats.purpose = purpose;
    result.stats.node_count = main_graph.number_of_nodes();
    result.stats.input_edge_count = main_graph.number_of_edges();
    result.main_to_local_edge.assign(main_graph.number_of_edges(), -1);

    for (const NodeSpec& node : main_graph.nodes()) {
        result.graph.add_node(node);
    }

    std::vector<bool> keep(main_graph.number_of_edges(), true);
    if (options.prune_min_duration_parallel) {
        std::map<EndpointStateKey, std::size_t> representative_by_key;
        for (std::size_t edge_id = 0; edge_id < main_graph.edges().size(); ++edge_id) {
            const EdgeRecord& edge = main_graph.edges()[edge_id];
            const EndpointStateKey key{
                edge.u,
                edge.v,
                edge.data.start_state,
                edge.data.end_state,
            };
            const auto found = representative_by_key.find(key);
            if (found == representative_by_key.end()) {
                representative_by_key.emplace(key, edge_id);
                continue;
            }
            const std::size_t incumbent_id = found->second;
            if (better_duration_representative(
                    edge,
                    edge_id,
                    main_graph.edges()[incumbent_id],
                    incumbent_id)) {
                keep[incumbent_id] = false;
                found->second = edge_id;
            } else {
                keep[edge_id] = false;
            }
            ++result.stats.removed_parallel_edge_count;
        }
    }

    if (options.prune_empty_state_connectors) {
        for (std::size_t edge_id = 0; edge_id < main_graph.edges().size(); ++edge_id) {
            if (!keep[edge_id]) {
                continue;
            }
            if (is_empty_delivery_pickup_connector(
                    main_graph,
                    main_graph.edges()[edge_id])) {
                keep[edge_id] = false;
                ++result.stats.removed_empty_connector_count;
            }
        }
    }

    for (std::size_t main_edge_id = 0;
         main_edge_id < main_graph.edges().size();
         ++main_edge_id) {
        if (!keep[main_edge_id]) {
            continue;
        }
        const EdgeRecord& edge = main_graph.edges()[main_edge_id];
        const int local_edge_id =
            static_cast<int>(result.graph.number_of_edges());
        result.graph.add_edge(edge.u, edge.v, edge.data);
        result.local_to_main_edge.push_back(main_edge_id);
        result.main_to_local_edge[main_edge_id] = local_edge_id;
    }

    result.stats.final_edge_count = result.graph.number_of_edges();
    result.stats.fingerprint = multigraph_fingerprint(result.graph);
    result.stats.build_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - build_start
    ).count();
    return result;
}

}  // namespace spdp
