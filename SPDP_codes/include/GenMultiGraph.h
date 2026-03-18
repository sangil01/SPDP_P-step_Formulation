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

// 두 칸짜리 차량 상태 표현에서 한 칸을 나타내는 토큰.
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

// 두 개의 순서 있는 상태 슬롯으로 표현한 차량 상태.
using State = std::array<StateToken, 2>;

// MultiDiGraph 안의 한 노드를 설명하는 메타데이터.
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

// MultiDiGraph의 한 edge에 저장되는 속성 정보.
struct EdgeData {
    std::vector<int> sequence_pi;
    double time = 0.0;
    double cost = 0.0;
    State start_state{};
    State end_state{};
};

// 양 끝 노드와 (u,v) 내 고유 key를 포함한 실제 edge 레코드.
struct EdgeRecord {
    NodeId u = 0;
    NodeId v = 0;
    int key = 0;  //같은 (u,v) 쌍에 대해 여러 edge가 있을 수 있으므로, 구분하기 위한 key
    EdgeData data;
};

// 비교와 중복 제거에 사용하는 정규화된 edge 표현.
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

// SPDP state-space의 노드와 edge를 저장하는 방향성 MultiDiGraph.
class MultiDiGraph {
public:
    // node id를 기준으로 노드 정보를 추가한다.
    void add_node(const NodeSpec& node);
    // edge를 추가하고 같은 (u, v) 쌍 안에서 고유 key를 자동 부여한다.
    void add_edge(NodeId u, NodeId v, const EdgeData& edge_data);

    // 저장된 노드 개수를 반환한다.
    std::size_t number_of_nodes() const;
    // 저장된 edge 개수를 반환한다.
    std::size_t number_of_edges() const;
    // 그래프에서 종료 노드로 사용하는 node id를 반환한다.
    NodeId end_node_id() const;
    // 시작/종료 노드를 제외한 실제 서비스 노드인지 확인한다.
    bool is_physical_service_node(NodeId node_id) const;

    // id로 노드 정보에 접근한다.
    const NodeSpec& node(NodeId node_id) const;
    // 삽입된 순서대로 전체 edge 레코드를 반환한다.
    const std::vector<EdgeRecord>& edges() const;
    // 노드 u에서 나가는 edge들의 인덱스 목록을 반환한다.
    const std::vector<std::size_t>& outgoing_edge_indices(NodeId u) const;
    // 노드 v로 들어오는 edge들의 인덱스 목록을 반환한다.
    const std::vector<std::size_t>& ingoing_edge_indices(NodeId v) const;
    // 노드 u에서 나가는 모든 edge를 반환한다.
    std::vector<EdgeRecord> edges_from(NodeId u) const;
    // 노드 v로 들어오는 모든 edge를 반환한다.
    std::vector<EdgeRecord> edges_to(NodeId v) const;
    // 순서쌍 (u, v) 사이의 모든 edge를 반환한다.
    std::vector<EdgeRecord> edges_between(NodeId u, NodeId v) const;

private:
    // (u, v) 쌍 기준 key 맵에서 사용하는 해시 함수.
    struct PairHash {
        std::size_t operator()(const std::pair<NodeId, NodeId>& value) const;
    };

    std::unordered_map<NodeId, NodeSpec> nodes_;
    std::vector<EdgeRecord> edges_;
    std::unordered_map<std::pair<NodeId, NodeId>, int, PairHash> next_key_by_pair_; //(u,v) 쌍에 대해 여러 edge가 있을 수 있으므로, 각각 구분하기 위한 추가 key 요소
    std::unordered_map<NodeId, std::vector<std::size_t>> outgoing_edge_indices_; //u에서 나가는 edge들의 edges_ 내 index 목록
    std::unordered_map<NodeId, std::vector<std::size_t>> ingoing_edge_indices_; //v로 들어오는 edge들의 edges_ 내 index 목록
};

// 필요하면 infeasible edge와 dominated edge를 제거하며 SPDP state-space MultiDiGraph를 생성한다.
MultiDiGraph build_multigraph(
    const SPDPData& data,
    bool prune_infeasible_edges = true,
    bool prune_dominated_edges = true,
    std::ostream* log_stream = nullptr
);

// 상태를 로깅과 디버깅에 쓰기 쉬운 문자열로 변환한다.
std::string state_to_str(const State& state);

}  // namespace spdp

#endif
