#ifndef SPDP_PSTEP_FORMULATION_H
#define SPDP_PSTEP_FORMULATION_H

#include <cstddef>
#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "gurobi_c++.h"
#include "GenMultiGraph.h"
#include "ReadData.h"

namespace spdp {

// compact 집계 전 단계의 feasible raw p-step 경로.
struct RawPStepPath {
    int id = 0;
    std::vector<int> edge_ids;
    int q = 0; // raw p-step 길이 (edge 개수), <= p
    NodeId start_node_id = 0;
    NodeId last_node_id = 0;
    State start_state{};
    State last_state{};
    double total_time = 0.0;
    double total_cost = 0.0;
    std::vector<NodeId> node_sequence;
    std::vector<State> state_sequence;
};

// master problem에서 사용하는 compact 형태의 p-step 표현.
struct CompactPStep {
    int id = 0;
    int raw_path_id = 0;
    std::vector<int> edge_ids;
    int q = 0;
    NodeId start_node_id = 0;
    NodeId last_node_id = 0;
    State start_state{};
    State last_state{};
    double total_time = 0.0;
    double total_cost = 0.0;
    double tau = 0.0;
    std::vector<NodeId> node_sequence;
    std::vector<State> state_sequence;
};

// 특정 노드와 하나의 차량 상태를 함께 식별하는 key.
struct NodeStateKey {
    NodeId node_id = 0;
    State state{};

    bool operator<(const NodeStateKey& other) const;
    bool operator==(const NodeStateKey& other) const;
};

// 모델 생성을 위해 compact p-step들로부터 만든 희소 계수 테이블.
struct CompactPStepCoefficients {
    std::vector<NodeId> physical_nodes; // 실제 서비스 노드 id 목록 (1, 2, ..., end_node_id-1)
    std::map<NodeId, std::vector<State>> sigma_by_node; // 각 physical node에서 등장하는 상태 집합 sigma_i

    // p-step 관점의 계수 표현.
    // 즉 "주어진 p-step r가 어떤 node / (node,state) / edge에 어떤 계수를 가지는가"를 저장한다.
    // col 기반 표현으로, p-step id 순서대로 벡터에 저장한다.
    std::vector<std::map<NodeId, int>> visit_coefficients_by_pstep;  // a_{i,r} 계수, vector index = p-step id
    std::vector<std::map<NodeStateKey, int>> state_coefficients_by_pstep; // s_{i,sigma,r} 계수, vector index = p-step id
    std::vector<std::map<NodeStateKey, double>> time_coefficients_by_pstep; // q_{i,sigma,r} 계수, vector index = p-step id
    std::vector<std::vector<int>> edge_incidence_by_pstep; // b_{e,r} 계수, vector index = p-step id, inner vector는 해당 p-step에 포함된 edge id 목록


    // row(node 또는 node,state) 관점의 희소 표현.
    // "주어진 제약식 row에 어떤 p-step들이 어떤 계수로 등장하는가"를 저장한다.
    std::map<NodeId, std::vector<std::pair<int, int>>> visit_rows;  //각 node에 대해, a_{i,r} 계수가 0이 아닌 (p-step id, 계수) 쌍 목록
    std::map<NodeStateKey, std::vector<std::pair<int, int>>> state_rows; //각 (node,state)에 대해, s_{i,sigma,r} 계수가 0이 아닌 (p-step id, 계수) 쌍 목록
    std::map<NodeStateKey, std::vector<std::pair<int, double>>> time_rows; //각 (node,state)에 대해, q_{i,sigma,r} 계수가 0이 아닌 (p-step id, 계수) 쌍 목록
    std::vector<std::vector<int>> edge_rows; //각 edge_id(=_edges의 index)에 대해, b_{e,r} 계수가 0이 아닌(1인) p-step id 목록
};

// 이후 최적화 코드에서 사용하는 p-step 생성 결과 묶음.
struct CompactPStepArtifacts {
    int p = 0;
    double time_limit = 0.0;
    std::vector<RawPStepPath> raw_paths;
    std::vector<CompactPStep> compact_psteps;
    CompactPStepCoefficients coefficients;
};

// p-step 생성과 검증, 로깅 동작을 제어하는 옵션.
struct CompactPStepOptions {
    int p = 2;
    double time_limit = -1.0;
    std::size_t dump_limit = 10;
    bool validate = true;
};

// Gurobi master problem 객체와 변수 핸들 모음.
struct CompactMasterProblem {
    std::unique_ptr<GRBEnv> env;
    std::unique_ptr<GRBModel> model;
    std::vector<GRBVar> x_vars;
    std::vector<GRBVar> theta_vars;
};

// 사람이 읽기 쉬운 해 복원 결과에서의 한 route action.
struct RouteSolutionAction {
    int request_id = 0;
    int location = 0;
    int pattern = 0;
    int load = 0;
};

// incumbent master problem 해에서 복원한 하나의 route 정보.
struct RecoveredRouteSolution {
    int route_id = 0;
    double total_time = 0.0;
    double total_cost = 0.0;
    std::vector<int> edge_ids;
    std::vector<RouteSolutionAction> actions;
};

// 목적함수값, 실행 시간, route 상세를 포함한 복원 해 전체 정보.
struct RecoveredSolution {
    bool has_incumbent = false;
    double objective_value = 0.0;
    double runtime_seconds = 0.0;
    std::vector<int> active_edge_ids;
    std::vector<RecoveredRouteSolution> routes;
};

// 시간 제한 아래에서 길이 p 이하의 feasible raw p-step 경로를 모두 나열한다.
std::vector<RawPStepPath> enumerate_feasible_raw_psteps(
    const MultiDiGraph& graph,
    int p,
    double time_limit
);

// raw path를 tau 값과 종료 노드 조건을 반영한 compact p-step으로 변환한다.
std::vector<CompactPStep> build_compact_psteps(
    const std::vector<RawPStepPath>& raw_paths,
    double time_limit,
    NodeId end_node_id
);

// compact master formulation에 필요한 희소 계수 구조를 생성한다.
CompactPStepCoefficients build_compact_pstep_coefficients(
    const MultiDiGraph& graph,
    const std::vector<CompactPStep>& compact_psteps
);

// raw path, compact p-step, 계수 테이블을 한 번에 생성한다.
CompactPStepArtifacts build_compact_pstep_artifacts(
    const MultiDiGraph& graph,
    const CompactPStepOptions& options,
    std::ostream* log_stream = nullptr
);

// graph와 p-step 결과를 바탕으로 Gurobi compact master problem을 만든다.
CompactMasterProblem build_compact_master_problem(
    const MultiDiGraph& graph,
    const CompactPStepArtifacts& artifacts,
    const std::string& log_path = ""
);

// raw p-step 나열 결과가 graph와 입력 제한을 만족하는지 검증한다.
void validate_raw_psteps(
    const MultiDiGraph& graph,
    const std::vector<RawPStepPath>& raw_paths,
    int p,
    double time_limit
);

// compact p-step 생성 결과를 검증한다.
void validate_compact_psteps(
    const std::vector<CompactPStep>& compact_psteps,
    double time_limit,
    NodeId end_node_id
);

// compact p-step으로부터 만든 계수 테이블을 검증한다.
void validate_compact_pstep_coefficients(
    const MultiDiGraph& graph,
    const std::vector<CompactPStep>& compact_psteps,
    const CompactPStepCoefficients& coefficients
);

// 생성된 compact p-step 목록의 앞부분을 출력한다.
void dump_compact_psteps(
    std::ostream& out,
    const std::vector<CompactPStep>& compact_psteps,
    std::size_t limit
);

// compact master problem을 제자리에서 최적화한다.
void solve_compact_master_problem(CompactMasterProblem& problem);

// 최적화된 master problem으로부터 route 수준 해 정보를 복원한다.
RecoveredSolution recover_incumbent_solution(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const CompactPStepArtifacts& artifacts,
    const CompactMasterProblem& problem
);

// 복원된 incumbent 해를 읽기 쉬운 텍스트 형식으로 출력한다.
void write_recovered_solution(std::ostream& out, const RecoveredSolution& solution);

}  // namespace spdp

#endif
