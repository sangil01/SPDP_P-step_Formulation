#ifndef SPDP_PSTEP_MASTER_CG_H
#define SPDP_PSTEP_MASTER_CG_H

#include <cstdint>
#include <cstddef>
#include <iosfwd>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "gurobi_c++.h"
#include "GenMultiGraph.h"
#include "PstepPricing.h"
#include "ReadData.h"

namespace spdp {

enum class GurobiLPMethod {
    Automatic = -1,
    Primal = 0,
    Dual = 1,
    Barrier = 2,
    Concurrent = 3,
};

// B&P node LP master construction / solve 옵션.
struct CGMasterOptions {
    // Gurobi 로그 파일 경로. 비어 있으면 기본 Gurobi 로그를 사용한다.
    std::string gurobi_log_path;

    // Gurobi thread 수. 음수면 Gurobi 기본값을 사용한다.
    int gurobi_threads = -1;

    // 현재 node LP solve에 허용할 시간 제한.
    double solver_time_limit = 3600.0;

    // 초기 LP method. 이후 solve 시점에 override 가능하다.
    GurobiLPMethod initial_lp_method = GurobiLPMethod::Primal;

    // optional theta fixing:
    // fixed_theta_value_by_edge[e] = -1이면 free, 0이면 theta_e = 0, 1이면 theta_e = 1.
    std::vector<int> fixed_theta_value_by_edge;
};

// node RMP에서 고정으로 존재하는 모든 row / variable handle 묶음.
// column generation 동안 x-column만 계속 추가되고, 나머지는 그대로 유지된다.
struct CGMasterProblem {
    // Gurobi environment / model.
    std::unique_ptr<GRBEnv> env;
    std::unique_ptr<GRBModel> model;

    // edge-linking용 theta_e 변수.
    // node LP relaxation이므로 binary가 아니라 [0,1] continuous로 둔다.
    std::vector<GRBVar> theta_vars;

    // visit row artificial z_i.
    // Phase I에서는 아직 real column이 못 메운 방문량을 대신 채운다.
    // Phase II에서는 upper bound 0으로 고정한다.
    std::vector<GRBVar> artificial_visit_vars;

    // state equality artificial:
    //   sum_r s_{i,sigma}^r x_r + z^+_{i,sigma} - z^-_{i,sigma} = 0
    // Phase I에서만 사용하고 Phase II에서 0으로 고정한다.
    std::map<NodeStateKey, GRBVar> artificial_state_pos_vars;
    std::map<NodeStateKey, GRBVar> artificial_state_neg_vars;

    // time inequality artificial:
    //   sum_r q_{i,sigma}^r x_r + h_{i,sigma} >= 0
    // real column 조합이 아직 부족할 때 음의 lhs를 임시로 보정한다.
    std::map<NodeStateKey, GRBVar> artificial_time_vars;

    // dynamic x-column 변수 목록.
    std::vector<GRBVar> x_vars;

    // master에 실제로 추가된 column 데이터.
    std::vector<CGColumn> columns;

    // duplicate column 삽입을 막기 위한 canonical key.
    std::map<std::string, int> column_id_by_key;

    // visit/state/time/linking row handle.
    std::map<NodeId, GRBConstr> visit_rows;
    std::map<NodeStateKey, GRBConstr> state_rows;
    std::map<NodeStateKey, GRBConstr> time_rows;
    std::vector<GRBConstr> edge_rows;

    // graph 전체에서 가능한 sigma 집합. state/time row를 upfront로 고정 생성할 때 사용한다.
    std::map<NodeId, std::vector<State>> sigma_by_node;

    // theta fixing과 phase-I artificial.
    std::vector<int> fixed_theta_value_by_edge;
    std::map<std::size_t, GRBVar> artificial_fixed_theta_vars;

    // 현재 phase.
    CGPhase phase = CGPhase::PhaseI;
};

// 한 번의 LP optimize 이후 필요한 요약 정보.
struct CGLPSnapshot {
    int gurobi_status = 0;
    bool has_primal_solution = false;
    double objective_value = 0.0;
    double objective_bound = 0.0;
    double runtime_seconds = 0.0;
    double artificial_sum = 0.0;
    int positive_x_count = 0;
    int active_theta_count = 0;
};

// graph와 data만으로 node RMP를 빈 column 집합 상태에서 생성한다.
// Phase I를 위해 visit/state/time row에는 artificial을 넣고,
// edge-linking row는 원래 formulation 그대로 둔다.
CGMasterProblem build_cg_master_problem(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const CGMasterOptions& options
);

// 현재 phase objective에 맞춰 dynamic compact p-step column들을 RMP에 추가한다.
// 이미 존재하는 duplicate column은 건너뛴다.
std::size_t add_columns_to_cg_master(
    CGMasterProblem& problem,
    const std::vector<CGColumn>& columns
);

// 현재 RMP를 node LP relaxation으로 최적화한다.
void solve_cg_master_lp(CGMasterProblem& problem);
void solve_cg_master_lp(CGMasterProblem& problem, GurobiLPMethod method);

// 현재 LP dual을 pricing이 바로 사용할 수 있게 추출한다.
CGDualSolution extract_cg_master_duals(const CGMasterProblem& problem);

// 현재 LP 상태를 요약한다.
CGLPSnapshot capture_cg_master_snapshot(const CGMasterProblem& problem);

// Phase I가 끝난 뒤 artificial visit variable을 0으로 고정하고 원래 objective로 전환한다.
void switch_cg_master_to_phase_two(CGMasterProblem& problem);

// node LP 상태를 텍스트로 출력한다.
void write_cg_master_snapshot(std::ostream& out, const CGLPSnapshot& snapshot);

}  // namespace spdp

#endif
