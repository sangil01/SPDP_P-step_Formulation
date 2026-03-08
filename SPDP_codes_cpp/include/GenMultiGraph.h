#ifndef SPDP_GEN_MULTI_GRAPH_H
#define SPDP_GEN_MULTI_GRAPH_H

#include <array>
#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "ReadData.h"

namespace spdp {

using NodeId = int;

struct StateToken {
    char kind = 'N';
    int container_type = -1;
    int treatment_id = -1;

    bool operator==(const StateToken& other) const {
        return kind == other.kind && container_type == other.container_type &&
               treatment_id == other.treatment_id;
    }

    bool operator<(const StateToken& other) const {
        if (kind != other.kind) {
            return kind < other.kind;
        }
        if (container_type != other.container_type) {
            return container_type < other.container_type;
        }
        return treatment_id < other.treatment_id;
    }
};

using State = std::array<StateToken, 2>;

struct NodeSpec {
    enum class Kind {
        Start,
        End,
        Pickup,
        Delivery,
    };

    NodeId node_id = 0;
    Kind kind = Kind::Start;
    int location = 0;
    std::optional<int> request_idx;
    std::optional<int> container_type;
    std::optional<int> treatment_id;
};

struct EdgeData {
    std::vector<int> sequence_pi;
    double time = 0.0;
    double cost = 0.0;
    State start_state{};
    State end_state{};
    std::string start_state_str;
    std::string end_state_str;
};

struct EdgeRecord {
    NodeId u = 0;
    NodeId v = 0;
    int key = 0;
    EdgeData data;
};

struct CanonicalEdgeKey {
    NodeId u = 0;
    NodeId v = 0;
    std::vector<int> sequence_pi;
    double time = 0.0;
    double cost = 0.0;
    State start_state{};
    State end_state{};

    bool operator==(const CanonicalEdgeKey& other) const;
    bool operator<(const CanonicalEdgeKey& other) const;
};

class MultiDiGraph {
public:
    void add_node(const NodeSpec& node);
    void add_edge(NodeId u, NodeId v, const EdgeData& edge_data);

    std::size_t number_of_nodes() const;
    std::size_t number_of_edges() const;

    const std::vector<EdgeRecord>& edges() const;
    std::vector<EdgeRecord> edges_from(NodeId u) const;
    std::vector<EdgeRecord> edges_between(NodeId u, NodeId v) const;

private:
    struct PairHash {
        std::size_t operator()(const std::pair<NodeId, NodeId>& value) const;
    };

    std::unordered_map<NodeId, NodeSpec> nodes_;
    std::vector<EdgeRecord> edges_;
    std::unordered_map<std::pair<NodeId, NodeId>, int, PairHash> next_key_by_pair_;
    std::unordered_map<NodeId, std::vector<std::size_t>> outgoing_edge_indices_;
};

using LoggerFn = std::function<void(const std::string&)>;

MultiDiGraph build_multigraph(
    const SPDPData& data,
    bool prune_dominated_edges = true,
    LoggerFn logger = LoggerFn{}
);

std::vector<CanonicalEdgeKey> build_canonical_edge_keys(const MultiDiGraph& graph);

std::string state_to_str(const State& state);

}  // namespace spdp

#endif
