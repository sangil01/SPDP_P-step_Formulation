#ifndef SPDP_PSTEP_BNP_H
#define SPDP_PSTEP_BNP_H

#include <cstdint>
#include <cstddef>
#include <iosfwd>
#include <optional>
#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "PstepMasterCG.h"
#include "PstepPricing.h"
#include "ReadData.h"

namespace spdp {

enum class NodeCGPhaseOneMode {
    ExactCG, // Phase I RMP를 CG로 풀어서 node-feasible한 solution을 찾는다.
    HeuristicCG, // SPDP pattern path를 현재 p-step window로 쪼개 seed column을 만든다.
    HeuristicCG3Step, // 기존 3-step 전용 heuristic-cg seed generator.
};

enum class NodeCGPhaseTwoPricingMode {
    ExactPricing,
    HeuristicPricingThenExact,
    FullEnumeration,
};

enum class NodeCGPhaseTwoHeuristicEngine {
    Labeling,
    ShallowSearch,
};

enum class BnPBranchingRule {
    ClosestToHalf,
};

enum class BnPNodeSelectionRule {
    BestBound,
    DepthFirst,
};

enum class BnPTreeMode {
    RootOnly,
    FullTree,
};

// B&P node column generation solver 전체 제어 옵션.
struct NodeCGOptions {
    // pricing에 사용할 p.
    int p = 2;

    // instance time limit T.
    double time_limit = 0.0;

    // 현재 node LP solve 전체 time limit.
    double solver_time_limit = 3600.0;

    // Gurobi thread 수.
    int gurobi_threads = -1;

    // Phase I / II simplex method.
    GurobiLPMethod phase_one_lp_method = GurobiLPMethod::Primal;
    GurobiLPMethod phase_two_lp_method = GurobiLPMethod::Primal;

    // 한 phase에서 허용하는 최대 CG iteration 수.
    std::size_t max_iterations_per_phase = 1000;

    // exact pricing이 start class마다 유지할 negative column 수.
    std::size_t exact_pricing_max_columns_per_start = 1;

    // exact pricing 한 round에서 master에 추가할 최대 column 수.
    std::size_t exact_pricing_max_total_columns_per_round = 256;

    // reduced cost tolerance.
    double reduced_cost_tolerance = -1e-6;

    // pricing에서 pickup ordering symmetry를 적용할지 여부.
    bool prune_pickup_symmetry_43 = true;

    // pricing에서 delivery ordering symmetry를 적용할지 여부.
    bool prune_delivery_symmetry_43 = true;

    // Gurobi log file 경로.
    std::string gurobi_log_path;

    // Root node에만 theta-form VI 35/36/44를 추가할지 여부.
    // Full-tree B&P에서는 child node RMP를 가볍게 유지하기 위해 root input state에만 전달한다.
    bool add_root_vi_35 = false;
    bool add_root_vi_36_combined = false;
    std::size_t root_vi_36_subset_max_size = 1;  // 0 disables VI-36.
    bool add_root_vi_request_block_sec = false;
    std::size_t root_vi_request_block_sec_max_size = 2;  // 0 disables request-block SEC.
    bool add_root_vi_44 = false;
    VI44KMinOptions vi_44_k_min_options;

    // Phase I에서 node-feasible RMP를 만드는 방식.
    NodeCGPhaseOneMode phase_one_mode = NodeCGPhaseOneMode::ExactCG;

    // Phase I heuristic-cg seed incumbent construction parameters.
    // max_attempts = 0이면 2 * |request| attempts를 사용한다.
    std::size_t phase_one_heuristic_cg_max_attempts = 0;
    std::size_t phase_one_heuristic_cg_max_incumbents = 20;
    std::size_t phase_one_heuristic_cg_top_l = 3;
    double phase_one_heuristic_cg_weight_cost = 1.0;
    double phase_one_heuristic_cg_weight_time = 0.1;
    double phase_one_heuristic_cg_weight_saving = 0.5;
    unsigned int phase_one_heuristic_cg_random_seed = 1;

    // Phase II pricing 방식.
    NodeCGPhaseTwoPricingMode phase_two_pricing_mode =
        NodeCGPhaseTwoPricingMode::ExactPricing;

    // Phase II heuristic pricing에서 앞에서부터 몇 개의 start class만 볼지.
    std::size_t phase_two_heuristic_max_starts = 128;

    // Phase II heuristic pricing에서 탐색할 start class 수 비율.
    double phase_two_heuristic_start_ratio = 0.25;

    // Adaptive K_start ladder level 수.
    // 0이면 heuristic 단계에서 모든 start class를 한 번에 본다.
    std::size_t phase_two_heuristic_ladder_levels = 3;

    // Phase II heuristic pricing이 start class마다 유지할 negative column 수.
    std::size_t phase_two_heuristic_max_columns_per_start = 1;

    // Phase II heuristic pricing이 한 round에서 master에 넘길 최대 column 수.
    std::size_t phase_two_heuristic_max_total_columns = 256;

    // Phase II heuristic pricing에서 search cap을 output cap 대비 몇 배까지 허용할지.
    double phase_two_heuristic_search_column_ratio = 1.0;

    // Phase II heuristic pricing의 start score 방식.
    HeuristicStartScoreMode phase_two_heuristic_start_score_mode =
        HeuristicStartScoreMode::OneStepMin;

    // Phase II heuristic engine 선택.
    NodeCGPhaseTwoHeuristicEngine phase_two_heuristic_engine =
        NodeCGPhaseTwoHeuristicEngine::ShallowSearch;

    // Phase II heuristic labeling에서 현재 label마다 확장할 top-k next transition 수.
    // 0이면 exact/full forward labeling과 동일하게 모든 outgoing transition을 확장한다.
    std::size_t phase_two_labeling_top_k_next = 0;

    // 3-step shallow search branching cap.
    std::size_t phase_two_shallow_k1 = 8;
    std::size_t phase_two_shallow_k2 = 4;

    // Column pool 사용 여부와 제어 파라미터.
    bool phase_two_column_pool_enabled = false;
    std::size_t phase_two_column_pool_max_size = 5000;
    std::size_t phase_two_column_pool_max_reprice = 256;

    // Full-enumeration pool RC update engine and worker count.
    FullEnumerationRCUpdateMode full_enumeration_rc_update_mode =
        FullEnumerationRCUpdateMode::Sequential;
    FullEnumerationRCUpdateBackend full_enumeration_parallel_stage1_backend =
        FullEnumerationRCUpdateBackend::Custom;
    FullEnumerationRCUpdateBackend full_enumeration_parallel_stage2_backend =
        FullEnumerationRCUpdateBackend::Custom;
    std::size_t full_enumeration_rc_update_threads = 0;
    bool full_enumeration_rc_detail_log = true;
};

struct NodeCGInputState {
    std::vector<CGColumn> initial_columns;
    std::vector<std::uint8_t> forbidden_edge_mask;
    std::vector<std::uint8_t> forced_edge_mask;
    std::vector<int> fixed_theta_value_by_edge;
    std::vector<int> required_outgoing_edge_by_node_id;
    std::vector<int> required_incoming_edge_by_node_id;
    bool branch_fixing_infeasible = false;
    bool enable_phase_one_seed_generation = true;
    bool add_vi_35 = false;
    bool add_vi_36_combined = false;
    std::size_t vi_36_subset_max_size = 1;  // 0 disables VI-36.
    bool add_vi_request_block_sec = false;
    std::size_t vi_request_block_sec_max_size = 2;  // 0 disables request-block SEC.
    bool add_vi_44 = false;
    const FullEnumerationStaticPool* shared_full_enumeration_static_pool = nullptr;
};

// 한 CG iteration의 phase / LP / pricing 결과 요약.
struct NodeCGIterationLog {
    CGPhase phase = CGPhase::PhaseI;
    std::size_t iteration_index = 0;
    double lp_objective_value = 0.0;
    double artificial_sum = 0.0;
    double best_reduced_cost = 0.0;
    ForwardPricingStatus pricing_status = ForwardPricingStatus::NoCompleteLabel;
    double lp_runtime_seconds = 0.0;
    double pricing_runtime_seconds = 0.0;
    std::size_t added_column_count = 0;
    std::size_t complete_label_count = 0;
    std::size_t start_label_count = 0;
    std::size_t generated_label_count = 0;
    std::size_t surviving_label_count = 0;
    std::size_t dominated_label_count = 0;
    bool heuristic_pricing_attempted = false;
    bool exact_pricing_fallback_used = false;
    std::size_t heuristic_explored_start_count = 0;
    std::size_t heuristic_found_column_count = 0;
    std::size_t heuristic_available_start_count = 0;
    NodeCGPhaseTwoHeuristicEngine heuristic_engine =
        NodeCGPhaseTwoHeuristicEngine::ShallowSearch;
    std::size_t heuristic_max_starts = 0;
    double heuristic_start_ratio = 0.0;
    std::size_t heuristic_ladder_levels = 0;
    std::size_t heuristic_output_max_columns_per_start = 0;
    std::size_t heuristic_output_max_columns_total = 0;
    double heuristic_search_column_ratio = 1.0;
    std::size_t heuristic_effective_search_max_columns_per_start = 0;
    std::size_t heuristic_effective_search_max_columns_total = 0;
    std::size_t heuristic_total_negative_column_count = 0;
    std::size_t heuristic_per_start_search_cap_hit_count = 0;
    bool heuristic_global_search_cap_hit = false;
    std::size_t heuristic_labeling_top_k_next = 0;
    std::size_t heuristic_shallow_k1 = 0;
    std::size_t heuristic_shallow_k2 = 0;
    std::size_t heuristic_top_k_applied_label_count = 0;
    std::size_t heuristic_top_k_feasible_edges_before = 0;
    std::size_t heuristic_top_k_feasible_edges_after = 0;
    bool column_pool_enabled = false;
    std::size_t column_pool_size_before_reprice = 0;
    std::size_t column_pool_max_reprice = 0;
    std::size_t column_pool_found_column_count = 0;
    std::size_t column_pool_deferred_added_count = 0;
    std::size_t column_pool_size_after_update = 0;
};

// B&P node column generation 종료 결과.
struct NodeCGResult {
    NodeCGPhaseOneMode phase_one_mode = NodeCGPhaseOneMode::ExactCG;
    NodeCGPhaseTwoPricingMode phase_two_pricing_mode = NodeCGPhaseTwoPricingMode::ExactPricing;
    bool phase_one_feasible = false;
    bool solved_to_completion = false;
    bool phase_two_reached = false;

    double phase_one_objective = 0.0;
    double phase_two_objective = 0.0;

    std::size_t phase_one_iterations = 0;
    std::size_t phase_two_iterations = 0;
    std::size_t total_columns_added = 0;

    CGLPSnapshot final_snapshot;
    std::vector<NodeCGIterationLog> iteration_logs;

    // 최종 RMP에 존재하는 모든 column과 그 LP 값.
    std::vector<CGColumn> columns;
    std::vector<double> column_values;
    std::vector<double> theta_values;
};

struct BranchAndPriceOptions {
    NodeCGOptions node_cg_options;
    BnPTreeMode tree_mode = BnPTreeMode::FullTree;
    BnPBranchingRule branching_rule = BnPBranchingRule::ClosestToHalf;
    BnPNodeSelectionRule node_selection_rule = BnPNodeSelectionRule::BestBound;
    double theta_integrality_tolerance = 1e-6;
    double gap_tolerance = 1e-6;
    std::string instance_name;
    double initial_upper_bound = -1.0;
    std::vector<CGColumn> initial_incumbent_columns;
    std::optional<RecoveredSolution> initial_incumbent_solution;
    double initial_incumbent_value = -1.0;
};

struct BranchAndPriceResult {
    BnPTreeMode tree_mode = BnPTreeMode::FullTree;
    bool solved_to_optimality = false;
    bool hit_time_limit = false;
    bool has_incumbent = false;
    bool initial_incumbent_used = false;
    bool incumbent_updated = false;

    double incumbent_value = 0.0;
    double best_global_lower_bound = 0.0;
    double global_gap = 0.0;

    std::size_t processed_node_count = 0;
    std::size_t open_node_count = 0;
    std::size_t infeasible_prune_count = 0;
    std::size_t bound_prune_count = 0;
    std::size_t integer_prune_count = 0;

    std::size_t phase_one_iterations = 0;
    std::size_t phase_two_iterations = 0;

    double phase_one_lp_runtime_seconds = 0.0;
    double phase_two_lp_runtime_seconds = 0.0;
    double phase_one_pricing_runtime_seconds = 0.0;
    double phase_two_pricing_runtime_seconds = 0.0;
    double total_lp_runtime_seconds = 0.0;
    double total_pricing_runtime_seconds = 0.0;

    std::vector<double> incumbent_theta_values;
    std::vector<CGColumn> incumbent_columns;
    std::vector<double> incumbent_column_values;
    std::optional<RecoveredSolution> incumbent_recovered_solution;
};

// Uses the same backward route decomposition and extreme-tau construction as
// heuristic-cg to turn complete incumbent routes into root-RMP seed columns.
std::vector<CGColumn> build_initial_incumbent_cg_columns(
    const MultiDiGraph& graph,
    int p,
    double time_limit,
    const RecoveredSolution& incumbent
);

// two-phase column generation으로 한 B&P node의 LP relaxation을 푼다.
NodeCGResult solve_node_column_generation(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const NodeCGOptions& options,
    std::ostream* log_stream = nullptr
);

NodeCGResult solve_node_column_generation(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const NodeCGOptions& options,
    const NodeCGInputState& input_state,
    std::ostream* log_stream = nullptr
);

BranchAndPriceResult solve_branch_and_price(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const BranchAndPriceOptions& options,
    std::ostream* log_stream = nullptr
);

// node CG 결과를 요약 출력한다.
void write_node_cg_summary(std::ostream& out, const NodeCGResult& result);

// node LP에서 양의 값을 갖는 x-column과 theta edge들을 출력한다.
void write_node_cg_solution(
    std::ostream& out,
    const NodeCGResult& result,
    const MultiDiGraph& graph
);
void write_branch_and_price_summary(std::ostream& out, const BranchAndPriceResult& result);
void write_branch_and_price_solution(
    std::ostream& out,
    const BranchAndPriceResult& result,
    const MultiDiGraph& graph
);

// CLI/logging용 phase-one mode 이름.
const char* to_string(NodeCGPhaseOneMode mode);
const char* to_string(NodeCGPhaseTwoPricingMode mode);
const char* to_string(NodeCGPhaseTwoHeuristicEngine mode);
const char* to_string(BnPTreeMode mode);
const char* to_string(BnPBranchingRule rule);
const char* to_string(BnPNodeSelectionRule rule);

}  // namespace spdp

#endif
