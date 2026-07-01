#ifndef SPDP_PSTEP_PRICING_H
#define SPDP_PSTEP_PRICING_H

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
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

enum class HeuristicStartScoreMode {
    OneStepMin,
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

    // full-enumeration matrix pool에서 왔을 때의 variant row index.
    // 다른 pricing 경로에서 생성된 column이면 -1이다.
    int pool_variant_index = -1;

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
    // 한 pricing round에서 start class마다 몇 개의 best negative column을 유지할지. (heurstic/exact pricing 별로 다른 값)
    std::size_t max_columns_per_start = 1;

    // 전체 pricing round에서 RMP에 넘길 최대 column 수. (heurstic/exact pricing 별로 다른 값)
    std::size_t max_total_columns = 256;

    // reduced cost가 이 값보다 작을 때만 column으로 인정한다.
    double reduced_cost_tolerance = -1e-6;

    // symmetry-43 pickup ordering을 적용할지 여부.
    bool prune_pickup_symmetry_43 = true;

    // delivery ordering symmetry를 적용할지 여부.
    bool prune_delivery_symmetry_43 = true;

    // heuristic pricing pass를 수행할지 여부.
    bool heuristic_pricing = false;

    // heuristic pricing에서 앞에서부터 몇 개의 start class만 볼지.
    std::size_t heuristic_max_starts = 0;

    // heuristic pricing에서 탐색할 start class 수 비율.
    double heuristic_start_ratio = 0.0;

    // heuristic pricing에서 search cap을 output cap 대비 몇 배까지 허용할지.
    // 1.0이면 search cap과 output cap이 같고, 1.0보다 크면 overflow column이 deferred로 남을 수 있다.
    double heuristic_search_column_ratio = 1.0;

    // heuristic pricing에서 사용할 start ordering score 방식.
    HeuristicStartScoreMode heuristic_start_score_mode = HeuristicStartScoreMode::OneStepMin;

    // heuristic labeling에서 현재 label의 outgoing transition 중 score 기준 상위 몇 개만 확장할지.
    // 0이면 기존 exact/full forward labeling처럼 모든 outgoing transition을 확장한다.
    std::size_t labeling_top_k_next = 0;

    // 3-step shallow search용 first/second arc branching cap.
    std::size_t shallow_k1 = 8;
    std::size_t shallow_k2 = 4;

    // 비어 있지 않으면 이 start class subset만 순서대로 pricing한다.
    std::vector<int> explicit_start_node_state_indices;

    // branch-and-price node에서 금지된 original graph edge mask.
    // size = graph.number_of_edges(), 1이면 해당 edge를 사용하는 transition/column을 금지한다.
    std::vector<std::uint8_t> forbidden_edge_mask;

    // branch-and-price one-child consistency rule:
    // if required_outgoing_edge_by_node_id[i] = e, any pricing edge leaving node i must be e.
    // if required_incoming_edge_by_node_id[j] = e, any pricing edge entering node j must be e.
    std::vector<int> required_outgoing_edge_by_node_id;
    std::vector<int> required_incoming_edge_by_node_id;
};

enum class ForwardPricingStatus {
    ColumnsFound,
    NoNegativeColumn,
    NoCompleteLabel,
};

// forward pricing 실행 결과 요약.
struct ForwardPricingResult {

    // 현재 pricing round의 종료 상태: improving column을 찾았는지, complete label은 있었지만 음수 column이 없었는지, 아니면 complete label 자체가 없었는지를 나타낸다.
    ForwardPricingStatus status = ForwardPricingStatus::NoCompleteLabel;

    // 이번 pricing round에서 caller가 restricted master에 바로 추가할 최종 accepted negative reduced-cost columns이다.
    std::vector<CGColumn> columns;

    // 이번 pricing round에서 생성되었지만 global insertion cap 등으로 master에는 바로 넣지 못해 이후 재사용이나 pool 보관 대상으로 넘기는 columns이다.
    std::vector<CGColumn> deferred_columns;

    // 이번 pricing round에서 평가된 모든 complete candidate 중 가장 작은 reduced cost 값이다.
    double best_reduced_cost = 0.0;

    // 이번 pricing round에서 찾은 improving(negative reduced-cost) column 총 개수이다.
    std::size_t total_negative_column_count = 0;

    // start별 heuristic search cap에 도달해 탐색을 조기 종료한 start 수이다.
    std::size_t per_start_search_cap_hit_count = 0;

    // round 전체 heuristic search cap에 도달해 탐색을 조기 종료했는지 여부이다.
    bool hit_global_search_cap = false;

    // heuristic labeling에서 top-k next filtering을 적용한 label 수이다.
    std::size_t top_k_next_applied_label_count = 0;

    // top-k filtering 전에 feasible했던 outgoing edge 수의 누적합이다.
    std::size_t top_k_next_feasible_edges_before = 0;

    // top-k filtering 후 실제 확장 대상으로 남은 outgoing edge 수의 누적합이다.
    std::size_t top_k_next_feasible_edges_after = 0;

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

struct FullEnumerationPoolBuildOptions {
    int p = 2;
    double time_limit = 0.0;
    bool prune_pickup_symmetry_43 = true;
    bool prune_delivery_symmetry_43 = true;
};

struct FullEnumerationPoolBuildStats {
    std::size_t raw_path_count = 0;
    std::size_t compact_pstep_count = 0;
    std::size_t inactive_column_count = 0;
    std::size_t inactive_path_count = 0;
    std::size_t skipped_master_column_count = 0;
    std::size_t skipped_duplicate_column_count = 0;
    double runtime_seconds = 0.0;
};

enum class FullEnumerationRCUpdateMode {
    Sequential,
    Parallel,
};

enum class FullEnumerationRCUpdateBackend {
    Custom,
    OneMKL,
};

struct FullEnumerationStaticPool {
    struct CSRMatrix {
        std::vector<std::size_t> row_ptr;
        std::vector<int> column_index;
        std::vector<double> value;
        std::size_t column_count = 0;
    };

    double time_limit = 0.0;

    std::vector<NodeId> visit_node_ids;
    std::vector<NodeStateKey> node_state_keys;
    std::map<NodeId, int> visit_index_by_node;
    std::map<NodeStateKey, int> node_state_index_by_key;
    std::vector<int> dense_visit_row_index_by_node_id;
    std::vector<NodeId> graph_edge_source_node_id;
    std::vector<NodeId> graph_edge_target_node_id;

    std::vector<int> entry_q;
    std::vector<NodeId> entry_start_node_id;
    std::vector<NodeId> entry_last_node_id;
    std::vector<State> entry_start_state;
    std::vector<State> entry_last_state;
    std::vector<double> entry_total_time;
    std::vector<double> entry_total_cost;
    std::vector<int> entry_start_time_row_index;
    std::vector<int> entry_last_time_row_index;

    std::vector<std::size_t> entry_edge_ids_row_ptr;
    std::vector<int> entry_edge_ids;
    std::vector<std::size_t> entry_node_sequence_row_ptr;
    std::vector<NodeId> entry_node_sequence;
    std::vector<std::size_t> entry_state_sequence_row_ptr;
    std::vector<State> entry_state_sequence;

    // stage 1:
    // rows = entries, columns = [visit rows | state rows | edge rows]
    CSRMatrix stage1_matrix;
    std::size_t stage1_visit_column_count = 0;
    std::size_t stage1_state_column_count = 0;
    std::size_t stage1_edge_column_count = 0;

    // entry마다 최대 2개의 tau slot을 고정 배치한다.
    // variant row v = 2 * entry + slot, slot in {0,1}.
    std::vector<double> variant_tau;
    std::vector<std::uint8_t> variant_valid;

    // stage 2:
    // rows = variants, columns = [base_rc(entry) | time rows]
    CSRMatrix stage2_matrix;

    std::unordered_map<std::string, std::size_t> variant_index_by_key;
};

struct FullEnumerationPoolNodeState {
    std::vector<std::uint8_t> variant_active;
    std::vector<std::uint8_t> forbidden_entry;
    std::vector<std::size_t> entry_available_variant_count;
    std::size_t available_entry_count = 0;
    std::size_t available_variant_count = 0;
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

// p=3 phase-two heuristic용 shallow search pricing을 수행한다.
ForwardPricingResult run_phase_two_shallow_search(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    const ForwardPricingOptions& options
);

FullEnumerationStaticPool build_full_enumeration_static_pool(
    const MultiDiGraph& graph,
    const FullEnumerationPoolBuildOptions& options,
    FullEnumerationPoolBuildStats* stats = nullptr
);

FullEnumerationPoolNodeState build_full_enumeration_pool_node_state(
    const FullEnumerationStaticPool& static_pool,
    const std::map<std::string, int>& active_column_id_by_key,
    FullEnumerationPoolBuildStats* stats = nullptr
);

ForwardPricingResult run_full_enumeration_pool_pricing(
    const FullEnumerationStaticPool& static_pool,
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
    std::ostream* log_stream = nullptr
);

void mark_full_enumeration_pool_columns_active(
    const FullEnumerationStaticPool& static_pool,
    FullEnumerationPoolNodeState& node_state,
    const std::vector<CGColumn>& columns
);

void mark_full_enumeration_pool_entries_forbidden_by_edge_mask(
    const FullEnumerationStaticPool& static_pool,
    const std::vector<std::uint8_t>& forbidden_edge_mask,
    FullEnumerationPoolNodeState& node_state
);

void mark_full_enumeration_pool_entries_forbidden_by_required_edges(
    const FullEnumerationStaticPool& static_pool,
    const std::vector<int>& required_outgoing_edge_by_node_id,
    const std::vector<int>& required_incoming_edge_by_node_id,
    FullEnumerationPoolNodeState& node_state
);

// 현재 dual에서 heuristic start score 기준으로 start class를 정렬한다.
std::vector<int> build_heuristic_start_order(
    const MultiDiGraph& graph,
    const ForwardPricingContext& context,
    const CGDualSolution& dual_solution,
    CGPhase phase,
    HeuristicStartScoreMode mode
);

// sparse master coefficients를 이용해 column의 exact reduced cost를 계산한다.
double evaluate_column_reduced_cost(
    const CGDualSolution& dual_solution,
    CGPhase phase,
    const CGColumn& column
);

// 현재 dual에서 나온 dynamic column을 읽기 쉬운 텍스트로 출력한다.
void write_generated_columns(std::ostream& out, const std::vector<CGColumn>& columns);

// CLI/logging용 pricing status 이름.
const char* to_string(ForwardPricingStatus status);
const char* to_string(HeuristicStartScoreMode mode);
const char* to_string(FullEnumerationRCUpdateMode mode);
const char* to_string(FullEnumerationRCUpdateBackend backend);

}  // namespace spdp

#endif
