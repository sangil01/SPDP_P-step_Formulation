#ifndef SPDP_GEN_MULTI_GRAPH_H
#define SPDP_GEN_MULTI_GRAPH_H

#include <array>
#include <cstddef>
#include <iosfwd>
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
    int landfill_location = -1;

    bool operator==(const StateToken& other) const {
        return kind == other.kind && container_type == other.container_type &&
               landfill_location == other.landfill_location;
    }

    bool operator<(const StateToken& other) const {
        if (kind != other.kind) {
            return kind < other.kind;
        }
        if (container_type != other.container_type) {
            return container_type < other.container_type;
        }
        return landfill_location < other.landfill_location;
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
    std::optional<int> landfill_location;
};

struct EdgeData {
    std::vector<int> sequence_pi;
    double time = 0.0;
    double cost = 0.0;
    State start_state{};
    State end_state{};
};

struct EdgeRecord {
    NodeId u = 0;
    NodeId v = 0;
    int key = 0;  //같은 (u,v) 쌍에 대해 여러 edge가 있을 수 있으므로, 구분하기 위한 key
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
    std::vector<EdgeRecord> edges_to(NodeId v) const;
    std::vector<EdgeRecord> edges_between(NodeId u, NodeId v) const;

private:
    struct PairHash {
        std::size_t operator()(const std::pair<NodeId, NodeId>& value) const;
    };

    std::unordered_map<NodeId, NodeSpec> nodes_;
    std::vector<EdgeRecord> edges_;
    std::unordered_map<std::pair<NodeId, NodeId>, int, PairHash> next_key_by_pair_; //(u,v) 쌍에 대해 여러 edge가 있을 수 있으므로, 각각 구분하기 위한 추가 key 요소
    std::unordered_map<NodeId, std::vector<std::size_t>> outgoing_edge_indices_; //u에서 나가는 edge들의 edges_ 내 index 목록
    std::unordered_map<NodeId, std::vector<std::size_t>> ingoing_edge_indices_; //v로 들어오는 edge들의 edges_ 내 index 목록
};

MultiDiGraph build_multigraph(
    const SPDPData& data,
    bool prune_infeasible_edges = true,
    bool prune_dominated_edges = true,
    std::ostream* log_stream = nullptr
);

std::string state_to_str(const State& state);

}  // namespace spdp

#endif
