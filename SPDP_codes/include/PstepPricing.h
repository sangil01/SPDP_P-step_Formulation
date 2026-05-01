#ifndef SPDP_PSTEP_PRICING_H
#define SPDP_PSTEP_PRICING_H

#include <cstddef>
#include <iosfwd>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "GenMultiGraph.h"
#include "PstepFormulation.h"

namespace spdp {

// node column generation의 현재 phase를 나타낸다.
// Phase I에서는 방문 artificial을 없애는 것이 목적이고,
// Phase II에서는 원래 compact p-step objective를 사용한다.
enum class CGPhase {
    PhaseI,
    PhaseII,
};

// node LP dual을 pricing이 읽기 쉬운 형태로 묶은 구조체.
// visit / state / time / edge-linking row dual을 모두 저장한다.
struct CGDualSolution {
    // alpha_i: physical service node i의 visit row dual.
    std::map<NodeId, double> visit_duals;

    // beta_{i,sigma}: (node,state) state-balance row dual.
    std::map<NodeStateKey, double> state_duals;

    // gamma_{i,sigma}: (node,state) time row dual.
    std::map<NodeStateKey, double> time_duals;

    // delta_e: edge-linking row dual.
    std::vector<double> edge_duals;
};

// symmetry-43용 request equivalence class.
// 같은 (pickup location, treatment location, container type) class에서
// request index가 증가하는 순서만 허용한다.
struct PricingRequestClassKey {
    int pickup_location = 0;
    int treatment_location = 0;
    int container_type = 0;

    bool operator<(const PricingRequestClassKey& other) const {
        if (pickup_location != other.pickup_location) {
            return pickup_location < other.pickup_location;
        }
        if (treatment_location != other.treatment_location) {
            return treatment_location < other.treatment_location;
        }
        return container_type < other.container_type;
    }
};

// delivery symmetry용 request equivalence class.
// 같은 (delivery location, container type) class에서
// request index가 증가하는 순서만 허용한다.
struct PricingDeliveryClassKey {
    int delivery_location = 0;
    int container_type = 0;

    bool operator<(const PricingDeliveryClassKey& other) const {
        if (delivery_location != other.delivery_location) {
            return delivery_location < other.delivery_location;
        }
        return container_type < other.container_type;
    }
};

// pricing이 생성하는 dynamic compact p-step column 표현.
// 기존 CompactPStep과 거의 같은 경로 정보를 들고 가되,
// column generation에서 바로 RMP에 추가할 sparse coefficient도 함께 저장한다.
struct CGColumn {
    // master 안에서의 동적 column id. add 시점에 master가 채운다.
    int id = -1;

    // raw path 길이 q.
    int q = 0;

    // path의 시작 / 종료 node-state.
    NodeId start_node_id = 0;
    NodeId last_node_id = 0;
    State start_state{};
    State last_state{};

    // path 자체의 실제 time/cost.
    double total_time = 0.0;
    double total_cost = 0.0;

    // 현재 dual에서 reduced cost를 최소화하는 extreme tau.
    double tau = 0.0;

    // column generation 당시의 reduced cost.
    double reduced_cost = 0.0;

    // 경로 복원 및 디버깅용 sequence 정보.
    std::vector<int> edge_ids;
    std::vector<NodeId> node_sequence;
    std::vector<State> state_sequence;

    // master row에 추가할 sparse coefficient.
    std::vector<std::pair<NodeId, int>> visit_coefficients;
    std::vector<std::pair<NodeStateKey, int>> state_coefficients;
    std::vector<std::pair<NodeStateKey, double>> time_coefficients;

    // b_{e,r}=1인 edge id 목록.
    std::vector<int> edge_incidence;
};

// start-state 고정 pricing에 앞서 graph를 전처리한 결과.
// node CG 동안 바뀌지 않는 정보만 담는다.
struct ForwardPricingContext {
    // 어떤 p-step 길이를 pricing할지.
    int p = 0;

    // instance time limit T.
    double time_limit = 0.0;

    // graph의 종료 depot node id.
    NodeId end_node_id = 0;

    // 각 node-state를 pricing에서 빠르게 참조하기 위한 정수 index화 정보.
    struct NodeStateInfo {
        NodeId node_id = 0;
        State state{};
        bool is_physical_service_node = false;

        // physical service node면 0-based dense bit index, 아니면 -1.
        int physical_bit_index = -1;

        // pickup / delivery node면 각 symmetry용 class / request index를 갖는다.
        bool is_pickup = false;
        int pickup_request_class_index = -1;

        bool is_delivery = false;
        int delivery_request_class_index = -1;

        int request_index = -1;
    };

    // pricing edge 전용 경량 표현.
    struct EdgeInfo {
        int graph_edge_id = -1;
        int from_node_state_index = -1;
        int to_node_state_index = -1;
        double time = 0.0;
        double cost = 0.0;
    };

    // graph에서 등장하는 모든 node-state pair.
    std::vector<NodeStateInfo> node_states;

    // (node,state) -> dense index lookup.
    std::map<NodeStateKey, int> node_state_index_by_key;

    // pricing에서 사용하는 모든 directed edge.
    std::vector<EdgeInfo> edges;

    // 각 node-state에서 나가는 pricing edge index 목록. forward labeling에서 빠르게 다음 edge를 찾기 위해서.
    std::vector<std::vector<int>> outgoing_edges_by_node_state;

    // 각 node-state로 들어오는 pricing edge index 목록. backward labeling에서 빠르게 이전 edge를 찾기 위해서.
    std::vector<std::vector<int>> incoming_edges_by_node_state;

    // forward labeling을 위한 feasible start class u=(s,sigma_s) 목록.
    std::vector<int> start_node_state_indices;

    // backward labeling을 위한 feasible end class 목록.
    std::vector<int> end_node_state_indices;

    // physical service node 수와 pickup / delivery symmetry class 수.
    std::size_t physical_node_count = 0;
    std::size_t pickup_request_class_count = 0;
    std::size_t delivery_request_class_count = 0;

    // physical node별 가능한 sigma rows.
    std::map<NodeId, std::vector<State>> sigma_by_node;
};

// forward pricing 제어 옵션.
struct ForwardPricingOptions {
    // 한 pricing round에서 start class마다 몇 개의 best negative column을 유지할지.
    std::size_t max_columns_per_start = 1;

    // 전체 pricing round에서 RMP에 넘길 최대 column 수.
    std::size_t max_total_columns = 256;

    // reduced cost가 이 값보다 작을 때만 column으로 인정한다.
    double reduced_cost_tolerance = -1e-6;

    // symmetry-43 pickup ordering을 적용할지 여부.
    bool prune_pickup_symmetry_43 = true;

    // delivery ordering symmetry를 적용할지 여부.
    bool prune_delivery_symmetry_43 = true;
};

enum class ForwardPricingStatus {
    ColumnsFound,
    NoNegativeColumn,
    NoCompleteLabel,
};

// forward pricing 실행 결과 요약.
struct ForwardPricingResult {
    ForwardPricingStatus status = ForwardPricingStatus::NoCompleteLabel;
    std::vector<CGColumn> columns;
    double best_reduced_cost = 0.0;

    // complete path 조건을 만족해 column 후보 평가까지 간 label 수.
    std::size_t complete_label_count = 0;

    // pricing에서 실제 생성한 root/start label 수.
    std::size_t start_label_count = 0;

    // edge 확장으로 실제 생성된 child label 수.
    std::size_t generated_label_count = 0;

    // dominance bucket에 최종적으로 남아 있는 label 수.
    std::size_t surviving_label_count = 0;

    // dominance 검사에서 버려진 label 수
    // (새 label이 기존 label에 지배되거나, 기존 label이 새 label에 의해 제거된 경우 포함).
    std::size_t dominated_label_count = 0;

    //generated_labels + start_labels = surviving_labels + dominated_labels
};

// forward pricing에 필요한 graph / state-space 전처리 결과를 생성한다.
ForwardPricingContext build_forward_pricing_context(
    const MultiDiGraph& graph,
    int p,
    double time_limit
);

// 고정된 node LP dual에 대해 start-only forward pricing을 수행한다.
// 모든 feasible start class를 훑어 negative reduced-cost compact p-step column을 반환한다.
ForwardPricingResult run_forward_pricing(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    const ForwardPricingOptions& options
);

// 현재 dual에서 나온 dynamic column을 읽기 쉬운 텍스트로 출력한다.
void write_generated_columns(std::ostream& out, const std::vector<CGColumn>& columns);

// CLI/logging용 pricing status 이름.
const char* to_string(ForwardPricingStatus status);

}  // namespace spdp

#endif
