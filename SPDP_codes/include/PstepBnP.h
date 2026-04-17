#ifndef SPDP_PSTEP_BNP_H
#define SPDP_PSTEP_BNP_H

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "PstepMasterCG.h"
#include "PstepPricing.h"
#include "ReadData.h"

namespace spdp {

enum class NodeCGPhaseOneMode {
    ExactCG, // Phase I RMP를 CG로 풀어서 node-feasible한 solution을 찾는다.
    HeuristicSeed, // Phase I RMP를 푸는게 아니라 heuristic한 seed column들을 추가해서 빠르게 node-feasible한 solution을 찾는다. 즉 artificial variable을 아예 안 만든다.
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

    // 한 phase에서 허용하는 최대 CG iteration 수.
    std::size_t max_iterations_per_phase = 1000;

    // pricing이 start class마다 유지할 negative column 수.
    std::size_t max_columns_per_start = 1;

    // 한 pricing round에서 master에 추가할 최대 column 수.
    std::size_t max_total_columns_per_round = 256;

    // reduced cost tolerance.
    double reduced_cost_tolerance = -1e-6;

    // Gurobi log file 경로.
    std::string gurobi_log_path;

    // Phase I에서 node-feasible RMP를 만드는 방식.
    NodeCGPhaseOneMode phase_one_mode = NodeCGPhaseOneMode::ExactCG;
};

// 한 CG iteration의 phase / LP / pricing 결과 요약.
struct NodeCGIterationLog {
    CGPhase phase = CGPhase::PhaseI;
    std::size_t iteration_index = 0;
    double lp_objective_value = 0.0;
    double artificial_sum = 0.0;
    double best_reduced_cost = 0.0;
    ForwardPricingStatus pricing_status = ForwardPricingStatus::NoCompleteLabel;
    double pricing_runtime_seconds = 0.0;
    std::size_t added_column_count = 0;
    std::size_t complete_label_count = 0;
    std::size_t start_label_count = 0;
    std::size_t generated_label_count = 0;
    std::size_t surviving_label_count = 0;
    std::size_t dominated_label_count = 0;
};

// B&P node column generation 종료 결과.
struct NodeCGResult {
    NodeCGPhaseOneMode phase_one_mode = NodeCGPhaseOneMode::ExactCG;
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
};

// two-phase column generation으로 한 B&P node의 LP relaxation을 푼다.
NodeCGResult solve_node_column_generation(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const NodeCGOptions& options,
    std::ostream* log_stream = nullptr
);

// node CG 결과를 요약 출력한다.
void write_node_cg_summary(std::ostream& out, const NodeCGResult& result);

// node LP에서 양의 값을 갖는 x-column들을 출력한다.
void write_node_cg_solution(std::ostream& out, const NodeCGResult& result);

// CLI/logging용 phase-one mode 이름.
const char* to_string(NodeCGPhaseOneMode mode);

}  // namespace spdp

#endif
