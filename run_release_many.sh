#!/usr/bin/env bash
set -u

export GUROBI_HOME="/opt/gurobi1303/linux64"
export PATH="$GUROBI_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$GUROBI_HOME/lib"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/SPDP_codes/build-release/SPDP_P_step"
SUMMARY_EXE="$SCRIPT_DIR/generate_root_cg_summary_xlsx.py"
OUTPUT_DIR=SPDP_output

if [[ ! -x "$EXE" ]]; then
    echo "Release executable not found or not executable: $EXE"
    exit 1
fi

# ============================================
# INPUT PARAMETERS
# ============================================
P=0 # 0 uses the direct two-index model when SOLVER_MODE=enumeration
SOLVER_MODE=enumeration # Options: enumeration, branch-and-price
ENUMERATION_SOS1_MODE=default # Options: default, sos1-auto, sos1-native (used only with SOLVER_MODE=enumeration)
ENUMERATION_MODEL_TYPE=ip # Options: ip, lp (lp is supported for P=0 direct two-index enumeration)
ENUMERATION_OBJECTIVE=original-cost  # Options: original-cost, travel-cost-only, duration, duration-plus-fixed (used only with SOLVER_MODE=enumeration)
ADD_FIXED_VEHICLE_NUMBER=0 # Uses the selected k_min; supported for P=0 direct two-index enumeration
BNP_TREE_MODE=root-only # Options: root-only, full-tree
NODE_CG_PHASE1_MODE=heuristic-cg # Options: exact-cg, heuristic-cg, heuristic-cg-3-step
NODE_CG_PHASE2_PRICING_MODE=full-enumeration # Options: exact-pricing, heuristic-pricing-then-exact, full-enumeration
SOLVER_TIME_LIMIT=3600 # seconds
INITIAL_INCUMBENT_ENABLE=0
INITIAL_INCUMBENT_TIME_LIMIT=1200 # seconds; 0 means no time limit
INITIAL_INCUMBENT_MAX_K_INCREMENTS=0 # additional K values after the initial lower bound
INITIAL_INCUMBENT_TIMEOUT_ACTION=advance # Options: stop, advance
INITIAL_INCUMBENT_BACKEND=two-index-milp # Options: two-index-milp, cp-sat
INITIAL_INCUMBENT_MILP_MODE=duration # Options: duration, feasibility, makespan
INITIAL_INCUMBENT_MILP_MAKESPAN_HORIZON_FACTOR=1.5 # Used only by makespan
INITIAL_INCUMBENT_MILP_DIRECT_DFF_IDENTITY_ENABLE=0
INITIAL_INCUMBENT_MILP_DIRECT_DFF_FS_ENABLE=0
DFF_FS_LAMBDA_LIST="0.1,0.2,0.3,0.4" # Each value must satisfy 0 < lambda < 0.5; lambda=0 is identity.
INITIAL_INCUMBENT_CP_GRAPH_MODE=original-graph # Options: original-graph, multigraph
INITIAL_INCUMBENT_CP_WORKERS=0
INITIAL_INCUMBENT_CP_MODE=duration # Options: satisfaction, duration, threshold-optimization
INITIAL_INCUMBENT_CP_THRESHOLD_HORIZON_FACTOR=1.5 # Used only by threshold-optimization
INITIAL_INCUMBENT_CP_REDUNDANT_TERMINAL_BALANCE=1
INITIAL_INCUMBENT_CP_REDUNDANT_FULL_RESERVOIR=1
INITIAL_INCUMBENT_CP_REDUNDANT_CONTAINER_WORKLOAD=1
INITIAL_INCUMBENT_CP_REDUNDANT_AGGREGATE_DURATION=1
INITIAL_INCUMBENT_CP_SYMMETRY_FIRST_PICKUP=1
INITIAL_INCUMBENT_CP_SYMMETRY_43=1
GUROBI_THREADS=0
NODE_CG_PHASE1_LP_METHOD=automatic # Options: automatic, primal, dual, barrier, concurrent
NODE_CG_PHASE2_LP_METHOD=primal # Options: automatic, primal, dual, barrier, concurrent
BNP_BRANCHING_RULE=closest-to-half # Options: closest-to-half
BNP_NODE_SELECTION_RULE=dfs # Options: best-bound, dfs
BNP_THETA_INTEGRALITY_TOLERANCE=1e-6
BNP_GAP_TOLERANCE=1e-6
BNP_INITIAL_UPPER_BOUND=-1 # 음수면 instance별 baseline UB 사용
VI_FORMULATION=theta # Options: theta, x
DUMP_PSTEPS=10
VALIDATE_PSTEPS=0
SOLVE_MODEL=1
PRUNE_INFEASIBLE_EDGES=1
PRUNE_DOMINATED_EDGES=1
PRUNE_SYMMETRY_40=1
PRUNE_SYMMETRY_41=1
PRUNE_PICKUP_SYMMETRY_43=1
PRUNE_DELIVERY_SYMMETRY_43=1
DURATION_GRAPH_PRUNE_MIN_TIME_PARALLEL=1
DURATION_GRAPH_PRUNE_EMPTY_STATE_CONNECTORS=1
INITIAL_INCUMBENT_GRAPH_PRUNE_MIN_TIME_PARALLEL=1
ADD_VI_35=0
ADD_VI_36_COMBINED=0
VI_36_SUBSET_MAX_SIZE=0 # 0이면 비활성화, 1이면 기존 singleton Eq.(36)과 동일, 7이면 현재 instance들에 대해 전체 set까지 포함
ADD_VI_Request_BLOCK_SEC=0
VI_Request_BLOCK_SEC_MAX_SIZE=0 # 0이면 비활성화, 활성화할 때는 2 이상
ADD_VI_44=1
# ============================================
# CUT-STUDY OPTIONS (added on top of main)
# Every cut-management knob implemented in the solver is exposed here.
# ============================================

# --- overall budget ---
INSTANCE_TIME_LIMIT=0 # seconds; wall-clock budget of one instance shared by preprocessing, the VI44 duration IP, the travel-time cover duration IPs, the incumbent MILP, the pre-MIP root LP and the main MIP (0 disables the global budget)

# --- main MILP duration modelling (independent; at least one must be active) ---
ADD_TIME_FLOW_FORMULATION=0 # 0: off, 1: node conservation (TF1), 2: node-state conservation (TF2)
ADD_BIG_M_TIME_CONSTRAINTS=1 # 0: drop the big-M route-duration rows, 1: keep them

# --- root cut processing ---
ROOT_CUT_MODE=callback # Options: pre-mip (cut loop on a standalone LP, rows copied into the MIP), callback (separate at the MIP root node)
PRE_MIP_ROOT_LP_TIME_LIMIT=30 # seconds; budget of the pre-mip LP phase, always deducted from the main MIP budget
PRE_MIP_METHOD=-1 # Gurobi Method of the pre-mip LP: -1 automatic, 0 primal simplex, 1 dual simplex, 2 barrier, 4 deterministic concurrent
PRE_MIP_CROSSOVER=-1 # -1 Gurobi default, 0 no crossover (allowed only with ROOT_CUT_MODE=pre-mip and PRE_MIP_METHOD=2)
MIP_METHOD=-1 # Gurobi Method of the main MIP root relaxation: -1, 0, 1, 2, 4
MIP_NODE_METHOD=-1 # Gurobi NodeMethod of the main MIP: -1, 0, 1, 2
MIP_CROSSOVER=-1 # -1 Gurobi default, 0 no crossover (allowed only with MIP_METHOD=2 and MIP_NODE_METHOD=2)

# --- capacity blossom cuts: candidate generation ---
ADD_CAPACITY_BLOSSOM_CUTS=0 # 1 separates odd-set capacity cuts for pickup/delivery sets (exact Padberg-Rao) as user cuts
CAPACITY_BLOSSOM_ROW_FORM=internal # Options: internal, inbound
CAPACITY_BLOSSOM_SUPPORT_TOLERANCE=1e-6 # edges below this weight are dropped from the separation support graph
# --- capacity blossom cuts: cut management ---
CAPACITY_BLOSSOM_SCOPE=adaptive-tree # Options: root-only, adaptive-tree, full-tree
CAPACITY_BLOSSOM_ROOT_MAX_ROUNDS=20
CAPACITY_BLOSSOM_ROOT_MAX_CUTS=256
CAPACITY_BLOSSOM_ROOT_MAX_PER_ROUND=32
CAPACITY_BLOSSOM_TREE_MAX_PER_ROUND=8
CAPACITY_BLOSSOM_TREE_DENSE_NODE_LIMIT=50
CAPACITY_BLOSSOM_TREE_FREQUENCY=100
CAPACITY_BLOSSOM_MAX_TOTAL=2000
CAPACITY_BLOSSOM_MIN_VIOLATION=1e-4
CAPACITY_BLOSSOM_TIME_FRACTION=0.15 # stop tree separation above this fraction of the solver runtime (0 disables)
CAPACITY_BLOSSOM_MAX_OVERLAP_JACCARD=0.9 # skip a candidate overlapping an accepted one of the same round above this index
CAPACITY_BLOSSOM_ROOT_NO_CUT_ROUND_LIMIT=2 # close the root after this many consecutive rounds without a cut (0 disables)
CAPACITY_BLOSSOM_ROOT_LOW_VIOLATION=0.02
CAPACITY_BLOSSOM_ROOT_LOW_VIOLATION_ROUND_LIMIT=2 # close the root after this many consecutive rounds below the threshold (0 disables)

# --- type-closed travel-time cover cuts: candidate generation ---
ADD_TYPE_TRAVEL_TIME_COVER_CUTS=0 # 1 adds sum_{delta^-(S_H)} y >= rho(H) for type-closed sets S_H (certified duration IPs on metric-closed restricted instances)
TYPE_TRAVEL_TIME_COVER_MODE=screened # Options: static (all certified rows up front), screened (add violated rows as user cuts)
TYPE_TRAVEL_TIME_COVER_MAX_TYPE_SET_SIZE=0 # 0: every union of types
TYPE_TRAVEL_TIME_COVER_SUBPROBLEM_TIME_LIMIT=20 # seconds per restricted duration IP
TYPE_TRAVEL_TIME_COVER_TOTAL_TIME_LIMIT=120 # seconds for the whole preprocessing
TYPE_TRAVEL_TIME_COVER_MIN_RHO=1 # rows with rho(H) below this value are dropped
TYPE_TRAVEL_TIME_COVER_THREADS=1
# --- type-closed travel-time cover cuts: cut management (screened mode only) ---
TYPE_TRAVEL_TIME_COVER_SCOPE=root-only # Options: root-only, adaptive-tree, full-tree
TYPE_TRAVEL_TIME_COVER_ROOT_MAX_ROUNDS=20
TYPE_TRAVEL_TIME_COVER_ROOT_MAX_CUTS=256
TYPE_TRAVEL_TIME_COVER_ROOT_MAX_PER_ROUND=8
TYPE_TRAVEL_TIME_COVER_TREE_MAX_PER_ROUND=4
TYPE_TRAVEL_TIME_COVER_TREE_DENSE_NODE_LIMIT=50
TYPE_TRAVEL_TIME_COVER_TREE_FREQUENCY=100
TYPE_TRAVEL_TIME_COVER_MAX_TOTAL=2000
TYPE_TRAVEL_TIME_COVER_MIN_VIOLATION=1e-4
TYPE_TRAVEL_TIME_COVER_TIME_FRACTION=0.15
TYPE_TRAVEL_TIME_COVER_MAX_OVERLAP_JACCARD=0.9
TYPE_TRAVEL_TIME_COVER_ROOT_NO_CUT_ROUND_LIMIT=2
TYPE_TRAVEL_TIME_COVER_ROOT_LOW_VIOLATION=0.02
TYPE_TRAVEL_TIME_COVER_ROOT_LOW_VIOLATION_ROUND_LIMIT=2

# ============ end of cut-study options ============
VI_44_K_MIN_USE_COR=1
VI_44_K_MIN_USE_SUBPROBLEM=1
VI_44_K_MIN_USE_VEHICLE_ASSIGNMENT=0
VI_44_SUBPROBLEM_TYPE=ip # Options: lp, ip
VI_44_SUBPROBLEM_TIME_LIMIT=600 # seconds; 0 means no time limit
VI_44_DURATION_IP_ROUNDED_BOUND_STOP=1
VI_44_SUBPROBLEM_DFF_FS_ENABLE=0
VI_44_VEHICLE_ASSIGNMENT_ADD_TSP_BOUND=1
VI_44_VEHICLE_ASSIGNMENT_ADD_CONTAINER_BOUND=1
VI_44_VEHICLE_ASSIGNMENT_TIME_LIMIT=120 # seconds; 0 means no time limit
CG_MAX_ITERATIONS_PER_PHASE=1000
EXACT_PRICING_MAX_COLUMNS_PER_START=64
EXACT_PRICING_MAX_TOTAL_COLUMNS_PER_ROUND=4096
CG_REDUCED_COST_TOLERANCE=-1e-6
PHASE1_HEURISTIC_CG_MAX_ATTEMPTS=0 # 0이면 2 * request_count
PHASE1_HEURISTIC_CG_MAX_INCUMBENTS=20
PHASE1_HEURISTIC_CG_TOP_L=3
PHASE1_HEURISTIC_CG_WEIGHT_COST=1.0
PHASE1_HEURISTIC_CG_WEIGHT_TIME=0.1
PHASE1_HEURISTIC_CG_WEIGHT_SAVING=0.5
PHASE1_HEURISTIC_CG_RANDOM_SEED=1
PHASE2_HEURISTIC_MAX_STARTS=4096  #1024면 모든 인스턴스에서 |x,\sigma(s)|보다 큼 (A11에서 701개임)
PHASE2_HEURISTIC_START_RATIO=1
PHASE2_HEURISTIC_LADDER_LEVELS=0 #0이면 모든 starte (node,state) 조합을 탐색. (위에 MAX_STARTS와 상관없이)
PHASE2_HEURISTIC_MAX_COLUMNS_PER_START=8
PHASE2_HEURISTIC_MAX_COLUMNS_TOTAL=512
PHASE2_HEURISTIC_SEARCH_COLUMN_RATIO=1.0
PHASE2_HEURISTIC_START_SCORE_MODE=one-step-min # Options: one-step-min
PHASE2_HEURISTIC_ENGINE=labeling # Options: labeling, shallow-search
PHASE2_LABELING_TOP_K_NEXT=4 #0이면 기존 exact/full forward labeling처럼 모든 outgoing transition을 확장한다. 4면 상위 4개만 확장한다.
PHASE2_SHALLOW_K1=4
PHASE2_SHALLOW_K2=2
PHASE2_COLUMN_POOL_ENABLE=0
PHASE2_COLUMN_POOL_MAX_SIZE=100000
PHASE2_COLUMN_POOL_MAX_REPRICE=10000
FULL_ENUMERATION_RC_UPDATE_MODE=sequential # Options: sequential, parallel
FULL_ENUMERATION_PARALLEL_STAGE1_BACKEND=custom # Options: custom, onemkl
FULL_ENUMERATION_PARALLEL_STAGE2_BACKEND=custom # Options: custom, onemkl
FULL_ENUMERATION_RC_UPDATE_THREADS=0 # 0이면 hardware_concurrency 사용
FULL_ENUMERATION_RC_DETAIL_LOG=0 # 0이면 entry/variant RC 상세 로그 비활성화
DATA_LIST=(
    #=========Request 20 이하=========#
    "RecDep_day_A1.dat"
    "RecDep_day_A2.dat"
    "RecDep_day_A3.dat"
    "RecDep_day_A4.dat"
    "RecDep_day_A5.dat"
    "RecDep_day_A6.dat"
    "RecDep_day_A7.dat"
    "RecDep_day_A8.dat"
    "RecDep_day_A9.dat"
    "RecDep_day_A10.dat"
    "RecDep_day_A11.dat"
    "RecDep_day_B1.dat"
    "RecDep_day_B2.dat"
    "RecDep_day_C1.dat"
    "RecDep_day_C2.dat"
    "RecDep_day_C3.dat"
    "RecDep_day_C4.dat"
    #=================================#
    #=========Request 50 이하=========#
    "RecDep_day_A12.dat"
    "RecDep_day_A13.dat"
    "RecDep_day_A14.dat"
    "RecDep_day_A15.dat"
    "RecDep_day_A16.dat"
    "RecDep_day_A17.dat"
    "RecDep_day_A18.dat"
    "RecDep_day_A19.dat"
    "RecDep_day_A20.dat"
    "RecDep_day_B3.dat"
    "RecDep_day_B4.dat"
    "RecDep_day_B5.dat"
    "RecDep_day_B6.dat"
    "RecDep_day_B7.dat"
    "RecDep_day_B8.dat"
    "RecDep_day_B9.dat"
    "RecDep_day_B10.dat"
    "RecDep_day_B11.dat"
    "RecDep_day_B12.dat"
    "RecDep_day_B13.dat"
    "RecDep_day_B14.dat"
    "RecDep_day_C5.dat"
    "RecDep_day_C6.dat"
    "RecDep_day_C7.dat"
    "RecDep_day_C8.dat"
    "RecDep_day_C9.dat"
    "RecDep_day_C10.dat"
    "RecDep_day_C11.dat"
    "RecDep_day_C12.dat"
    "RecDep_day_D1.dat"
    "RecDep_day_D2.dat"
    #=================================#
    #=====Request 50 초과 100 이하=====#
    "RecDep_day_B15.dat"
    "RecDep_day_B16.dat"
    "RecDep_day_B17.dat"
    "RecDep_day_B18.dat"
    "RecDep_day_B19.dat"
    "RecDep_day_B20.dat"
    "RecDep_day_C13.dat"
    "RecDep_day_C14.dat"
    "RecDep_day_C15.dat"
    "RecDep_day_C16.dat"
    "RecDep_day_C17.dat"
    "RecDep_day_C18.dat"
    "RecDep_day_C19.dat"
    "RecDep_day_C20.dat"
    "RecDep_day_D3.dat"
    "RecDep_day_D4.dat"
    "RecDep_day_D5.dat"
    "RecDep_day_D6.dat"
    "RecDep_day_D7.dat"
    #=================================#
    #====Request 100 초과 200 이하====#
    '''"RecDep_day_D8.dat"
    "RecDep_day_D9.dat"
    "RecDep_day_D10.dat"
    "RecDep_day_D11.dat"
    "RecDep_day_D12.dat"
    "RecDep_day_D13.dat"
    "RecDep_day_D14.dat"
    "RecDep_day_D15.dat"
    "RecDep_day_D16.dat"
    "RecDep_day_D17.dat"
    "RecDep_day_D18.dat"
    "RecDep_day_D19.dat"
    "RecDep_day_D20.dat"'''
    #=================================#
)
# Put one data file name per line in DATA_LIST.
# ============================================

RUN_TAG_SUFFIX="${RUN_TAG_SUFFIX:-}"
ENUMERATION_SOS1_TAG=""
if [[ "$SOLVER_MODE" == "enumeration" ]]; then
    ENUMERATION_SOS1_TAG="_${ENUMERATION_SOS1_MODE}_${ENUMERATION_MODEL_TYPE}_${ENUMERATION_OBJECTIVE}-fk${ADD_FIXED_VEHICLE_NUMBER}"
fi
RUN_TAG="P${P}_${SOLVER_MODE}_${BNP_TREE_MODE}_${SOLVER_TIME_LIMIT}s_${NODE_CG_PHASE1_MODE}_${NODE_CG_PHASE2_PRICING_MODE}${ENUMERATION_SOS1_TAG}${RUN_TAG_SUFFIX}"
RUN_TAG+="_dgmp${DURATION_GRAPH_PRUNE_MIN_TIME_PARALLEL}-dges${DURATION_GRAPH_PRUNE_EMPTY_STATE_CONNECTORS}-igmp${INITIAL_INCUMBENT_GRAPH_PRUNE_MIN_TIME_PARALLEL}-rbs${VI_44_DURATION_IP_ROUNDED_BOUND_STOP}"
RUN_TAG+="-vifs${VI_44_SUBPROBLEM_DFF_FS_ENABLE}"
RUN_TAG+="-tf${ADD_TIME_FLOW_FORMULATION}-bigm${ADD_BIG_M_TIME_CONSTRAINTS}-cb${ADD_CAPACITY_BLOSSOM_CUTS}-ttc${ADD_TYPE_TRAVEL_TIME_COVER_CUTS}"
RUN_TAG+="-rc${ROOT_CUT_MODE}-pm${PRE_MIP_METHOD}-px${PRE_MIP_CROSSOVER}-mm${MIP_METHOD}-nm${MIP_NODE_METHOD}-mx${MIP_CROSSOVER}"
if [[ "$INITIAL_INCUMBENT_ENABLE" == "1" ]]; then
    RUN_TAG+="_initial-incumbent-${INITIAL_INCUMBENT_BACKEND}-${INITIAL_INCUMBENT_TIME_LIMIT}s-kinc${INITIAL_INCUMBENT_MAX_K_INCREMENTS}-timeout-${INITIAL_INCUMBENT_TIMEOUT_ACTION}"
    if [[ "$INITIAL_INCUMBENT_BACKEND" == "two-index-milp" ]]; then
        RUN_TAG+="-milpm${INITIAL_INCUMBENT_MILP_MODE}"
        RUN_TAG+="-did${INITIAL_INCUMBENT_MILP_DIRECT_DFF_IDENTITY_ENABLE}-dfs${INITIAL_INCUMBENT_MILP_DIRECT_DFF_FS_ENABLE}"
        if [[ "$INITIAL_INCUMBENT_MILP_MODE" == "makespan" ]]; then
            RUN_TAG+="-milpf${INITIAL_INCUMBENT_MILP_MAKESPAN_HORIZON_FACTOR}"
        fi
    fi
    if [[ "$INITIAL_INCUMBENT_BACKEND" == "cp-sat" ]]; then
        RUN_TAG+="-cpg${INITIAL_INCUMBENT_CP_GRAPH_MODE}"
        RUN_TAG+="-cpm${INITIAL_INCUMBENT_CP_MODE}"
        if [[ "$INITIAL_INCUMBENT_CP_MODE" == "threshold-optimization" ]]; then
            RUN_TAG+="-cpf${INITIAL_INCUMBENT_CP_THRESHOLD_HORIZON_FACTOR}"
        fi
        RUN_TAG+="-cpw${INITIAL_INCUMBENT_CP_WORKERS}-rtb${INITIAL_INCUMBENT_CP_REDUNDANT_TERMINAL_BALANCE}-rfr${INITIAL_INCUMBENT_CP_REDUNDANT_FULL_RESERVOIR}-rcw${INITIAL_INCUMBENT_CP_REDUNDANT_CONTAINER_WORKLOAD}-rad${INITIAL_INCUMBENT_CP_REDUNDANT_AGGREGATE_DURATION}-sfp${INITIAL_INCUMBENT_CP_SYMMETRY_FIRST_PICKUP}-s43${INITIAL_INCUMBENT_CP_SYMMETRY_43}"
    fi
fi
if [[ "$INITIAL_INCUMBENT_MILP_DIRECT_DFF_FS_ENABLE" == "1" || "$VI_44_SUBPROBLEM_DFF_FS_ENABLE" == "1" ]]; then
    DFF_FS_TAG="${DFF_FS_LAMBDA_LIST//,/_}"
    RUN_TAG+="-fsl${DFF_FS_TAG}"
fi
if [[ "$OUTPUT_DIR" = /* ]]; then
    OUTPUT_PATH="$OUTPUT_DIR"
else
    OUTPUT_PATH="$SCRIPT_DIR/$OUTPUT_DIR"
fi
mkdir -p "$OUTPUT_PATH"
RUNNER_LOG="$OUTPUT_PATH/${RUN_TAG}.log"
SUMMARY_XLSX="$OUTPUT_PATH/${RUN_TAG}.xlsx"

FAILED=0
TOTAL=${#DATA_LIST[@]}
CURRENT=0

describe_exit_code() {
    local exit_code=$1
    if (( exit_code >= 128 )); then
        local signal=$((exit_code - 128))
        if command -v kill >/dev/null 2>&1; then
            local signal_name
            signal_name="$(kill -l "$signal" 2>/dev/null || true)"
            if [[ -n "$signal_name" ]]; then
                printf 'signal %s (%s)' "$signal" "$signal_name"
                return
            fi
        fi
        printf 'signal %s' "$signal"
        return
    fi

    printf 'exit code %s' "$exit_code"
}

for data_name in "${DATA_LIST[@]}"; do
    CURRENT=$((CURRENT + 1))

    echo "=================================================="
    echo "Running data [$CURRENT/$TOTAL]: $data_name"
    echo "Runner log: $RUNNER_LOG"
    {
        echo "=================================================="
        echo "Started at: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "Running data [$CURRENT/$TOTAL]: $data_name"
        echo "p: $P"
        echo "output-dir: $OUTPUT_DIR"
        echo "solver-mode: $SOLVER_MODE"
        echo "enumeration-sos1-mode: $ENUMERATION_SOS1_MODE"
        echo "enumeration-model-type: $ENUMERATION_MODEL_TYPE"
        echo "enumeration-objective: $ENUMERATION_OBJECTIVE"
        echo "add-fixed-vehicle-number: $ADD_FIXED_VEHICLE_NUMBER"
        echo "bnp-tree-mode: $BNP_TREE_MODE"
        echo "node-cg-phase1-mode: $NODE_CG_PHASE1_MODE"
        echo "node-cg-phase2-pricing-mode: $NODE_CG_PHASE2_PRICING_MODE"
        echo "node-cg-phase1-lp-method: $NODE_CG_PHASE1_LP_METHOD"
        echo "node-cg-phase2-lp-method: $NODE_CG_PHASE2_LP_METHOD"
        echo "bnp-branching-rule: $BNP_BRANCHING_RULE"
        echo "bnp-node-selection-rule: $BNP_NODE_SELECTION_RULE"
        echo "bnp-theta-integrality-tolerance: $BNP_THETA_INTEGRALITY_TOLERANCE"
        echo "bnp-gap-tolerance: $BNP_GAP_TOLERANCE"
        echo "bnp-initial-upper-bound: $BNP_INITIAL_UPPER_BOUND"
        echo "initial-incumbent-enable: $INITIAL_INCUMBENT_ENABLE"
        echo "initial-incumbent-time-limit: $INITIAL_INCUMBENT_TIME_LIMIT"
        echo "initial-incumbent-max-k-increments: $INITIAL_INCUMBENT_MAX_K_INCREMENTS"
        echo "initial-incumbent-timeout-action: $INITIAL_INCUMBENT_TIMEOUT_ACTION"
        echo "initial-incumbent-backend: $INITIAL_INCUMBENT_BACKEND"
        echo "initial-incumbent-milp-mode: $INITIAL_INCUMBENT_MILP_MODE"
        echo "initial-incumbent-milp-makespan-horizon-factor: $INITIAL_INCUMBENT_MILP_MAKESPAN_HORIZON_FACTOR"
        echo "initial-incumbent-milp-direct-dff-identity-enable: $INITIAL_INCUMBENT_MILP_DIRECT_DFF_IDENTITY_ENABLE"
        echo "initial-incumbent-milp-direct-dff-fs-enable: $INITIAL_INCUMBENT_MILP_DIRECT_DFF_FS_ENABLE"
        echo "dff-fs-lambda-list: $DFF_FS_LAMBDA_LIST"
        echo "initial-incumbent-cp-graph-mode: $INITIAL_INCUMBENT_CP_GRAPH_MODE"
        echo "initial-incumbent-cp-workers: $INITIAL_INCUMBENT_CP_WORKERS"
        echo "initial-incumbent-cp-mode: $INITIAL_INCUMBENT_CP_MODE"
        echo "initial-incumbent-cp-threshold-horizon-factor: $INITIAL_INCUMBENT_CP_THRESHOLD_HORIZON_FACTOR"
        echo "initial-incumbent-cp-redundant-terminal-balance: $INITIAL_INCUMBENT_CP_REDUNDANT_TERMINAL_BALANCE"
        echo "initial-incumbent-cp-redundant-full-reservoir: $INITIAL_INCUMBENT_CP_REDUNDANT_FULL_RESERVOIR"
        echo "initial-incumbent-cp-redundant-container-workload: $INITIAL_INCUMBENT_CP_REDUNDANT_CONTAINER_WORKLOAD"
        echo "initial-incumbent-cp-redundant-aggregate-duration: $INITIAL_INCUMBENT_CP_REDUNDANT_AGGREGATE_DURATION"
        echo "initial-incumbent-cp-symmetry-first-pickup: $INITIAL_INCUMBENT_CP_SYMMETRY_FIRST_PICKUP"
        echo "initial-incumbent-cp-symmetry-43: $INITIAL_INCUMBENT_CP_SYMMETRY_43"
        echo "vi-formulation: $VI_FORMULATION"
        echo "cg-max-iterations-per-phase: $CG_MAX_ITERATIONS_PER_PHASE"
        echo "exact-pricing-max-columns-per-start: $EXACT_PRICING_MAX_COLUMNS_PER_START"
        echo "exact-pricing-max-total-columns-per-round: $EXACT_PRICING_MAX_TOTAL_COLUMNS_PER_ROUND"
        echo "cg-reduced-cost-tolerance: $CG_REDUCED_COST_TOLERANCE"
        echo "phase1-heuristic-cg-max-attempts: $PHASE1_HEURISTIC_CG_MAX_ATTEMPTS"
        echo "phase1-heuristic-cg-max-incumbents: $PHASE1_HEURISTIC_CG_MAX_INCUMBENTS"
        echo "phase1-heuristic-cg-top-l: $PHASE1_HEURISTIC_CG_TOP_L"
        echo "phase1-heuristic-cg-weight-cost: $PHASE1_HEURISTIC_CG_WEIGHT_COST"
        echo "phase1-heuristic-cg-weight-time: $PHASE1_HEURISTIC_CG_WEIGHT_TIME"
        echo "phase1-heuristic-cg-weight-saving: $PHASE1_HEURISTIC_CG_WEIGHT_SAVING"
        echo "phase1-heuristic-cg-random-seed: $PHASE1_HEURISTIC_CG_RANDOM_SEED"
        echo "phase2-heuristic-max-starts: $PHASE2_HEURISTIC_MAX_STARTS"
        echo "phase2-heuristic-start-ratio: $PHASE2_HEURISTIC_START_RATIO"
        echo "phase2-heuristic-ladder-levels: $PHASE2_HEURISTIC_LADDER_LEVELS"
        echo "phase2-heuristic-max-columns-per-start: $PHASE2_HEURISTIC_MAX_COLUMNS_PER_START"
        echo "phase2-heuristic-max-columns-total: $PHASE2_HEURISTIC_MAX_COLUMNS_TOTAL"
        echo "phase2-heuristic-search-column-ratio: $PHASE2_HEURISTIC_SEARCH_COLUMN_RATIO"
        echo "phase2-heuristic-start-score-mode: $PHASE2_HEURISTIC_START_SCORE_MODE"
        echo "phase2-heuristic-engine: $PHASE2_HEURISTIC_ENGINE"
        echo "phase2-labeling-top-k-next: $PHASE2_LABELING_TOP_K_NEXT"
        echo "phase2-shallow-k1: $PHASE2_SHALLOW_K1"
        echo "phase2-shallow-k2: $PHASE2_SHALLOW_K2"
        echo "phase2-column-pool-enable: $PHASE2_COLUMN_POOL_ENABLE"
        echo "phase2-column-pool-max-size: $PHASE2_COLUMN_POOL_MAX_SIZE"
        echo "phase2-column-pool-max-reprice: $PHASE2_COLUMN_POOL_MAX_REPRICE"
        echo "full-enumeration-rc-update-mode: $FULL_ENUMERATION_RC_UPDATE_MODE"
        echo "full-enumeration-parallel-stage1-backend: $FULL_ENUMERATION_PARALLEL_STAGE1_BACKEND"
        echo "full-enumeration-parallel-stage2-backend: $FULL_ENUMERATION_PARALLEL_STAGE2_BACKEND"
        echo "full-enumeration-rc-update-threads: $FULL_ENUMERATION_RC_UPDATE_THREADS"
        echo "full-enumeration-rc-detail-log: $FULL_ENUMERATION_RC_DETAIL_LOG"
        echo "prune-pickup-symmetry-43: $PRUNE_PICKUP_SYMMETRY_43"
        echo "prune-delivery-symmetry-43: $PRUNE_DELIVERY_SYMMETRY_43"
        echo "duration-graph-prune-min-time-parallel: $DURATION_GRAPH_PRUNE_MIN_TIME_PARALLEL"
        echo "duration-graph-prune-empty-state-connectors: $DURATION_GRAPH_PRUNE_EMPTY_STATE_CONNECTORS"
        echo "initial-incumbent-graph-prune-min-time-parallel: $INITIAL_INCUMBENT_GRAPH_PRUNE_MIN_TIME_PARALLEL"
        echo "add-vi-35: $ADD_VI_35"
        echo "add-vi-36-combined: $ADD_VI_36_COMBINED"
        echo "vi-36-subset-max-size: $VI_36_SUBSET_MAX_SIZE"
        echo "add-vi-request-block-sec: $ADD_VI_Request_BLOCK_SEC"
        echo "vi-request-block-sec-max-size: $VI_Request_BLOCK_SEC_MAX_SIZE"
        echo "add-vi-44: $ADD_VI_44"
        echo "--- cut-study options ---"
        echo "instance-time-limit: $INSTANCE_TIME_LIMIT"
        echo "add-time-flow-formulation: $ADD_TIME_FLOW_FORMULATION"
        echo "add-big-m-time-constraints: $ADD_BIG_M_TIME_CONSTRAINTS"
        echo "root-cut-mode: $ROOT_CUT_MODE"
        echo "pre-mip-root-lp-time-limit: $PRE_MIP_ROOT_LP_TIME_LIMIT"
        echo "pre-mip-method: $PRE_MIP_METHOD"
        echo "pre-mip-crossover: $PRE_MIP_CROSSOVER"
        echo "mip-method: $MIP_METHOD"
        echo "mip-node-method: $MIP_NODE_METHOD"
        echo "mip-crossover: $MIP_CROSSOVER"
        echo "add-capacity-blossom-cuts: $ADD_CAPACITY_BLOSSOM_CUTS"
        echo "capacity-blossom-row-form: $CAPACITY_BLOSSOM_ROW_FORM"
        echo "capacity-blossom-support-tolerance: $CAPACITY_BLOSSOM_SUPPORT_TOLERANCE"
        echo "capacity-blossom-scope: $CAPACITY_BLOSSOM_SCOPE"
        echo "capacity-blossom-root-max-rounds: $CAPACITY_BLOSSOM_ROOT_MAX_ROUNDS"
        echo "capacity-blossom-root-max-cuts: $CAPACITY_BLOSSOM_ROOT_MAX_CUTS"
        echo "capacity-blossom-root-max-per-round: $CAPACITY_BLOSSOM_ROOT_MAX_PER_ROUND"
        echo "capacity-blossom-tree-max-per-round: $CAPACITY_BLOSSOM_TREE_MAX_PER_ROUND"
        echo "capacity-blossom-tree-dense-node-limit: $CAPACITY_BLOSSOM_TREE_DENSE_NODE_LIMIT"
        echo "capacity-blossom-tree-frequency: $CAPACITY_BLOSSOM_TREE_FREQUENCY"
        echo "capacity-blossom-max-total: $CAPACITY_BLOSSOM_MAX_TOTAL"
        echo "capacity-blossom-min-violation: $CAPACITY_BLOSSOM_MIN_VIOLATION"
        echo "capacity-blossom-time-fraction: $CAPACITY_BLOSSOM_TIME_FRACTION"
        echo "capacity-blossom-max-overlap-jaccard: $CAPACITY_BLOSSOM_MAX_OVERLAP_JACCARD"
        echo "capacity-blossom-root-no-cut-round-limit: $CAPACITY_BLOSSOM_ROOT_NO_CUT_ROUND_LIMIT"
        echo "capacity-blossom-root-low-violation: $CAPACITY_BLOSSOM_ROOT_LOW_VIOLATION"
        echo "capacity-blossom-root-low-violation-round-limit: $CAPACITY_BLOSSOM_ROOT_LOW_VIOLATION_ROUND_LIMIT"
        echo "add-type-travel-time-cover-cuts: $ADD_TYPE_TRAVEL_TIME_COVER_CUTS"
        echo "type-travel-time-cover-mode: $TYPE_TRAVEL_TIME_COVER_MODE"
        echo "type-travel-time-cover-max-type-set-size: $TYPE_TRAVEL_TIME_COVER_MAX_TYPE_SET_SIZE"
        echo "type-travel-time-cover-subproblem-time-limit: $TYPE_TRAVEL_TIME_COVER_SUBPROBLEM_TIME_LIMIT"
        echo "type-travel-time-cover-total-time-limit: $TYPE_TRAVEL_TIME_COVER_TOTAL_TIME_LIMIT"
        echo "type-travel-time-cover-min-rho: $TYPE_TRAVEL_TIME_COVER_MIN_RHO"
        echo "type-travel-time-cover-threads: $TYPE_TRAVEL_TIME_COVER_THREADS"
        echo "type-travel-time-cover-scope: $TYPE_TRAVEL_TIME_COVER_SCOPE"
        echo "type-travel-time-cover-root-max-rounds: $TYPE_TRAVEL_TIME_COVER_ROOT_MAX_ROUNDS"
        echo "type-travel-time-cover-root-max-cuts: $TYPE_TRAVEL_TIME_COVER_ROOT_MAX_CUTS"
        echo "type-travel-time-cover-root-max-per-round: $TYPE_TRAVEL_TIME_COVER_ROOT_MAX_PER_ROUND"
        echo "type-travel-time-cover-tree-max-per-round: $TYPE_TRAVEL_TIME_COVER_TREE_MAX_PER_ROUND"
        echo "type-travel-time-cover-tree-dense-node-limit: $TYPE_TRAVEL_TIME_COVER_TREE_DENSE_NODE_LIMIT"
        echo "type-travel-time-cover-tree-frequency: $TYPE_TRAVEL_TIME_COVER_TREE_FREQUENCY"
        echo "type-travel-time-cover-max-total: $TYPE_TRAVEL_TIME_COVER_MAX_TOTAL"
        echo "type-travel-time-cover-min-violation: $TYPE_TRAVEL_TIME_COVER_MIN_VIOLATION"
        echo "type-travel-time-cover-time-fraction: $TYPE_TRAVEL_TIME_COVER_TIME_FRACTION"
        echo "type-travel-time-cover-max-overlap-jaccard: $TYPE_TRAVEL_TIME_COVER_MAX_OVERLAP_JACCARD"
        echo "type-travel-time-cover-root-no-cut-round-limit: $TYPE_TRAVEL_TIME_COVER_ROOT_NO_CUT_ROUND_LIMIT"
        echo "type-travel-time-cover-root-low-violation: $TYPE_TRAVEL_TIME_COVER_ROOT_LOW_VIOLATION"
        echo "type-travel-time-cover-root-low-violation-round-limit: $TYPE_TRAVEL_TIME_COVER_ROOT_LOW_VIOLATION_ROUND_LIMIT"
        echo "--- end of cut-study options ---"
        echo "vi-44-k-min-use-cor: $VI_44_K_MIN_USE_COR"
        echo "vi-44-k-min-use-subproblem: $VI_44_K_MIN_USE_SUBPROBLEM"
        echo "vi-44-k-min-use-vehicle-assignment: $VI_44_K_MIN_USE_VEHICLE_ASSIGNMENT"
        echo "vi-44-subproblem-type: $VI_44_SUBPROBLEM_TYPE"
        echo "vi-44-subproblem-time-limit: $VI_44_SUBPROBLEM_TIME_LIMIT"
        echo "vi-44-duration-ip-rounded-bound-stop: $VI_44_DURATION_IP_ROUNDED_BOUND_STOP"
        echo "vi-44-subproblem-dff-fs-enable: $VI_44_SUBPROBLEM_DFF_FS_ENABLE"
        echo "vi-44-vehicle-assignment-add-tsp-bound: $VI_44_VEHICLE_ASSIGNMENT_ADD_TSP_BOUND"
        echo "vi-44-vehicle-assignment-add-container-bound: $VI_44_VEHICLE_ASSIGNMENT_ADD_CONTAINER_BOUND"
        echo "vi-44-vehicle-assignment-time-limit: $VI_44_VEHICLE_ASSIGNMENT_TIME_LIMIT"
        echo "----------------------------------------"
    } >> "$RUNNER_LOG"

    if "$EXE" "$data_name" \
        --p "$P" \
        --output-dir "$OUTPUT_DIR" \
        --solver-mode "$SOLVER_MODE" \
        --enumeration-sos1-mode "$ENUMERATION_SOS1_MODE" \
        --enumeration-model-type "$ENUMERATION_MODEL_TYPE" \
        --enumeration-objective "$ENUMERATION_OBJECTIVE" \
        --add-fixed-vehicle-number "$ADD_FIXED_VEHICLE_NUMBER" \
        --bnp-tree-mode "$BNP_TREE_MODE" \
        --node-cg-phase1-mode "$NODE_CG_PHASE1_MODE" \
        --node-cg-phase2-pricing-mode "$NODE_CG_PHASE2_PRICING_MODE" \
        --node-cg-phase1-lp-method "$NODE_CG_PHASE1_LP_METHOD" \
        --node-cg-phase2-lp-method "$NODE_CG_PHASE2_LP_METHOD" \
        --bnp-branching-rule "$BNP_BRANCHING_RULE" \
        --bnp-node-selection-rule "$BNP_NODE_SELECTION_RULE" \
        --bnp-theta-integrality-tolerance "$BNP_THETA_INTEGRALITY_TOLERANCE" \
        --bnp-gap-tolerance "$BNP_GAP_TOLERANCE" \
        --bnp-initial-upper-bound "$BNP_INITIAL_UPPER_BOUND" \
        --solver-time-limit "$SOLVER_TIME_LIMIT" \
        --initial-incumbent-enable "$INITIAL_INCUMBENT_ENABLE" \
        --initial-incumbent-time-limit "$INITIAL_INCUMBENT_TIME_LIMIT" \
        --initial-incumbent-max-k-increments "$INITIAL_INCUMBENT_MAX_K_INCREMENTS" \
        --initial-incumbent-timeout-action "$INITIAL_INCUMBENT_TIMEOUT_ACTION" \
        --initial-incumbent-backend "$INITIAL_INCUMBENT_BACKEND" \
        --initial-incumbent-milp-mode "$INITIAL_INCUMBENT_MILP_MODE" \
        --initial-incumbent-milp-makespan-horizon-factor "$INITIAL_INCUMBENT_MILP_MAKESPAN_HORIZON_FACTOR" \
        --initial-incumbent-milp-direct-dff-identity-enable "$INITIAL_INCUMBENT_MILP_DIRECT_DFF_IDENTITY_ENABLE" \
        --initial-incumbent-milp-direct-dff-fs-enable "$INITIAL_INCUMBENT_MILP_DIRECT_DFF_FS_ENABLE" \
        --dff-fs-lambda-list "$DFF_FS_LAMBDA_LIST" \
        --initial-incumbent-cp-graph-mode "$INITIAL_INCUMBENT_CP_GRAPH_MODE" \
        --initial-incumbent-cp-workers "$INITIAL_INCUMBENT_CP_WORKERS" \
        --initial-incumbent-cp-mode "$INITIAL_INCUMBENT_CP_MODE" \
        --initial-incumbent-cp-threshold-horizon-factor "$INITIAL_INCUMBENT_CP_THRESHOLD_HORIZON_FACTOR" \
        --initial-incumbent-cp-redundant-terminal-balance "$INITIAL_INCUMBENT_CP_REDUNDANT_TERMINAL_BALANCE" \
        --initial-incumbent-cp-redundant-full-reservoir "$INITIAL_INCUMBENT_CP_REDUNDANT_FULL_RESERVOIR" \
        --initial-incumbent-cp-redundant-container-workload "$INITIAL_INCUMBENT_CP_REDUNDANT_CONTAINER_WORKLOAD" \
        --initial-incumbent-cp-redundant-aggregate-duration "$INITIAL_INCUMBENT_CP_REDUNDANT_AGGREGATE_DURATION" \
        --initial-incumbent-cp-symmetry-first-pickup "$INITIAL_INCUMBENT_CP_SYMMETRY_FIRST_PICKUP" \
        --initial-incumbent-cp-symmetry-43 "$INITIAL_INCUMBENT_CP_SYMMETRY_43" \
        --gurobi-threads "$GUROBI_THREADS" \
        --vi-formulation "$VI_FORMULATION" \
        --dump-psteps "$DUMP_PSTEPS" \
        --validate-psteps "$VALIDATE_PSTEPS" \
        --solve "$SOLVE_MODEL" \
        --prune-infeasible-edges "$PRUNE_INFEASIBLE_EDGES" \
        --prune-dominated-edges "$PRUNE_DOMINATED_EDGES" \
        --prune-symmetry-40 "$PRUNE_SYMMETRY_40" \
        --prune-symmetry-41 "$PRUNE_SYMMETRY_41" \
        --duration-graph-prune-min-time-parallel "$DURATION_GRAPH_PRUNE_MIN_TIME_PARALLEL" \
        --duration-graph-prune-empty-state-connectors "$DURATION_GRAPH_PRUNE_EMPTY_STATE_CONNECTORS" \
        --initial-incumbent-graph-prune-min-time-parallel "$INITIAL_INCUMBENT_GRAPH_PRUNE_MIN_TIME_PARALLEL" \
        --prune-pickup-symmetry-43 "$PRUNE_PICKUP_SYMMETRY_43" \
        --prune-delivery-symmetry-43 "$PRUNE_DELIVERY_SYMMETRY_43" \
        --add-vi-35 "$ADD_VI_35" \
        --add-vi-36-combined "$ADD_VI_36_COMBINED" \
        --vi-36-subset-max-size "$VI_36_SUBSET_MAX_SIZE" \
        --add-vi-request-block-sec "$ADD_VI_Request_BLOCK_SEC" \
        --vi-request-block-sec-max-size "$VI_Request_BLOCK_SEC_MAX_SIZE" \
        --add-vi-44 "$ADD_VI_44" \
        --instance-time-limit "$INSTANCE_TIME_LIMIT" \
        --add-time-flow-formulation "$ADD_TIME_FLOW_FORMULATION" \
        --add-big-m-time-constraints "$ADD_BIG_M_TIME_CONSTRAINTS" \
        --root-cut-mode "$ROOT_CUT_MODE" \
        --pre-mip-root-lp-time-limit "$PRE_MIP_ROOT_LP_TIME_LIMIT" \
        --pre-mip-method "$PRE_MIP_METHOD" \
        --pre-mip-crossover "$PRE_MIP_CROSSOVER" \
        --mip-method "$MIP_METHOD" \
        --mip-node-method "$MIP_NODE_METHOD" \
        --mip-crossover "$MIP_CROSSOVER" \
        --add-capacity-blossom-cuts "$ADD_CAPACITY_BLOSSOM_CUTS" \
        --capacity-blossom-row-form "$CAPACITY_BLOSSOM_ROW_FORM" \
        --capacity-blossom-support-tolerance "$CAPACITY_BLOSSOM_SUPPORT_TOLERANCE" \
        --capacity-blossom-scope "$CAPACITY_BLOSSOM_SCOPE" \
        --capacity-blossom-root-max-rounds "$CAPACITY_BLOSSOM_ROOT_MAX_ROUNDS" \
        --capacity-blossom-root-max-cuts "$CAPACITY_BLOSSOM_ROOT_MAX_CUTS" \
        --capacity-blossom-root-max-per-round "$CAPACITY_BLOSSOM_ROOT_MAX_PER_ROUND" \
        --capacity-blossom-tree-max-per-round "$CAPACITY_BLOSSOM_TREE_MAX_PER_ROUND" \
        --capacity-blossom-tree-dense-node-limit "$CAPACITY_BLOSSOM_TREE_DENSE_NODE_LIMIT" \
        --capacity-blossom-tree-frequency "$CAPACITY_BLOSSOM_TREE_FREQUENCY" \
        --capacity-blossom-max-total "$CAPACITY_BLOSSOM_MAX_TOTAL" \
        --capacity-blossom-min-violation "$CAPACITY_BLOSSOM_MIN_VIOLATION" \
        --capacity-blossom-time-fraction "$CAPACITY_BLOSSOM_TIME_FRACTION" \
        --capacity-blossom-max-overlap-jaccard "$CAPACITY_BLOSSOM_MAX_OVERLAP_JACCARD" \
        --capacity-blossom-root-no-cut-round-limit "$CAPACITY_BLOSSOM_ROOT_NO_CUT_ROUND_LIMIT" \
        --capacity-blossom-root-low-violation "$CAPACITY_BLOSSOM_ROOT_LOW_VIOLATION" \
        --capacity-blossom-root-low-violation-round-limit "$CAPACITY_BLOSSOM_ROOT_LOW_VIOLATION_ROUND_LIMIT" \
        --add-type-travel-time-cover-cuts "$ADD_TYPE_TRAVEL_TIME_COVER_CUTS" \
        --type-travel-time-cover-mode "$TYPE_TRAVEL_TIME_COVER_MODE" \
        --type-travel-time-cover-max-type-set-size "$TYPE_TRAVEL_TIME_COVER_MAX_TYPE_SET_SIZE" \
        --type-travel-time-cover-subproblem-time-limit "$TYPE_TRAVEL_TIME_COVER_SUBPROBLEM_TIME_LIMIT" \
        --type-travel-time-cover-total-time-limit "$TYPE_TRAVEL_TIME_COVER_TOTAL_TIME_LIMIT" \
        --type-travel-time-cover-min-rho "$TYPE_TRAVEL_TIME_COVER_MIN_RHO" \
        --type-travel-time-cover-threads "$TYPE_TRAVEL_TIME_COVER_THREADS" \
        --type-travel-time-cover-scope "$TYPE_TRAVEL_TIME_COVER_SCOPE" \
        --type-travel-time-cover-root-max-rounds "$TYPE_TRAVEL_TIME_COVER_ROOT_MAX_ROUNDS" \
        --type-travel-time-cover-root-max-cuts "$TYPE_TRAVEL_TIME_COVER_ROOT_MAX_CUTS" \
        --type-travel-time-cover-root-max-per-round "$TYPE_TRAVEL_TIME_COVER_ROOT_MAX_PER_ROUND" \
        --type-travel-time-cover-tree-max-per-round "$TYPE_TRAVEL_TIME_COVER_TREE_MAX_PER_ROUND" \
        --type-travel-time-cover-tree-dense-node-limit "$TYPE_TRAVEL_TIME_COVER_TREE_DENSE_NODE_LIMIT" \
        --type-travel-time-cover-tree-frequency "$TYPE_TRAVEL_TIME_COVER_TREE_FREQUENCY" \
        --type-travel-time-cover-max-total "$TYPE_TRAVEL_TIME_COVER_MAX_TOTAL" \
        --type-travel-time-cover-min-violation "$TYPE_TRAVEL_TIME_COVER_MIN_VIOLATION" \
        --type-travel-time-cover-time-fraction "$TYPE_TRAVEL_TIME_COVER_TIME_FRACTION" \
        --type-travel-time-cover-max-overlap-jaccard "$TYPE_TRAVEL_TIME_COVER_MAX_OVERLAP_JACCARD" \
        --type-travel-time-cover-root-no-cut-round-limit "$TYPE_TRAVEL_TIME_COVER_ROOT_NO_CUT_ROUND_LIMIT" \
        --type-travel-time-cover-root-low-violation "$TYPE_TRAVEL_TIME_COVER_ROOT_LOW_VIOLATION" \
        --type-travel-time-cover-root-low-violation-round-limit "$TYPE_TRAVEL_TIME_COVER_ROOT_LOW_VIOLATION_ROUND_LIMIT" \
        --vi-44-k-min-use-cor "$VI_44_K_MIN_USE_COR" \
        --vi-44-k-min-use-subproblem "$VI_44_K_MIN_USE_SUBPROBLEM" \
        --vi-44-k-min-use-vehicle-assignment "$VI_44_K_MIN_USE_VEHICLE_ASSIGNMENT" \
        --vi-44-subproblem-type "$VI_44_SUBPROBLEM_TYPE" \
        --vi-44-subproblem-time-limit "$VI_44_SUBPROBLEM_TIME_LIMIT" \
        --vi-44-duration-ip-rounded-bound-stop "$VI_44_DURATION_IP_ROUNDED_BOUND_STOP" \
        --vi-44-subproblem-dff-fs-enable "$VI_44_SUBPROBLEM_DFF_FS_ENABLE" \
        --vi-44-vehicle-assignment-add-tsp-bound "$VI_44_VEHICLE_ASSIGNMENT_ADD_TSP_BOUND" \
        --vi-44-vehicle-assignment-add-container-bound "$VI_44_VEHICLE_ASSIGNMENT_ADD_CONTAINER_BOUND" \
        --vi-44-vehicle-assignment-time-limit "$VI_44_VEHICLE_ASSIGNMENT_TIME_LIMIT" \
        --cg-max-iterations-per-phase "$CG_MAX_ITERATIONS_PER_PHASE" \
        --exact-pricing-max-columns-per-start "$EXACT_PRICING_MAX_COLUMNS_PER_START" \
        --exact-pricing-max-total-columns-per-round "$EXACT_PRICING_MAX_TOTAL_COLUMNS_PER_ROUND" \
        --cg-reduced-cost-tolerance "$CG_REDUCED_COST_TOLERANCE" \
        --phase1-heuristic-cg-max-attempts "$PHASE1_HEURISTIC_CG_MAX_ATTEMPTS" \
        --phase1-heuristic-cg-max-incumbents "$PHASE1_HEURISTIC_CG_MAX_INCUMBENTS" \
        --phase1-heuristic-cg-top-l "$PHASE1_HEURISTIC_CG_TOP_L" \
        --phase1-heuristic-cg-weight-cost "$PHASE1_HEURISTIC_CG_WEIGHT_COST" \
        --phase1-heuristic-cg-weight-time "$PHASE1_HEURISTIC_CG_WEIGHT_TIME" \
        --phase1-heuristic-cg-weight-saving "$PHASE1_HEURISTIC_CG_WEIGHT_SAVING" \
        --phase1-heuristic-cg-random-seed "$PHASE1_HEURISTIC_CG_RANDOM_SEED" \
        --phase2-heuristic-max-starts "$PHASE2_HEURISTIC_MAX_STARTS" \
        --phase2-heuristic-start-ratio "$PHASE2_HEURISTIC_START_RATIO" \
        --phase2-heuristic-ladder-levels "$PHASE2_HEURISTIC_LADDER_LEVELS" \
        --phase2-heuristic-max-columns-per-start "$PHASE2_HEURISTIC_MAX_COLUMNS_PER_START" \
        --phase2-heuristic-max-columns-total "$PHASE2_HEURISTIC_MAX_COLUMNS_TOTAL" \
        --phase2-heuristic-search-column-ratio "$PHASE2_HEURISTIC_SEARCH_COLUMN_RATIO" \
        --phase2-heuristic-start-score-mode "$PHASE2_HEURISTIC_START_SCORE_MODE" \
        --phase2-heuristic-engine "$PHASE2_HEURISTIC_ENGINE" \
        --phase2-labeling-top-k-next "$PHASE2_LABELING_TOP_K_NEXT" \
        --phase2-shallow-k1 "$PHASE2_SHALLOW_K1" \
        --phase2-shallow-k2 "$PHASE2_SHALLOW_K2" \
        --phase2-column-pool-enable "$PHASE2_COLUMN_POOL_ENABLE" \
        --phase2-column-pool-max-size "$PHASE2_COLUMN_POOL_MAX_SIZE" \
        --phase2-column-pool-max-reprice "$PHASE2_COLUMN_POOL_MAX_REPRICE" \
        --full-enumeration-rc-update-mode "$FULL_ENUMERATION_RC_UPDATE_MODE" \
        --full-enumeration-parallel-stage1-backend "$FULL_ENUMERATION_PARALLEL_STAGE1_BACKEND" \
        --full-enumeration-parallel-stage2-backend "$FULL_ENUMERATION_PARALLEL_STAGE2_BACKEND" \
        --full-enumeration-rc-update-threads "$FULL_ENUMERATION_RC_UPDATE_THREADS" \
        --full-enumeration-rc-detail-log "$FULL_ENUMERATION_RC_DETAIL_LOG" >> "$RUNNER_LOG" 2>&1; then
        echo "Completed: $data_name"
        {
            echo "Finished at: $(date '+%Y-%m-%d %H:%M:%S')"
            echo "Result: success"
        } >> "$RUNNER_LOG"
    else
        exit_code=$?
        exit_summary="$(describe_exit_code "$exit_code")"
        echo "Failed: $data_name ($exit_summary)"
        {
            echo "Finished at: $(date '+%Y-%m-%d %H:%M:%S')"
            echo "Result: failure ($exit_summary)"
        } >> "$RUNNER_LOG"
        FAILED=$((FAILED + 1))
    fi
done

if [[ "$SOLVER_MODE" == "branch-and-price" ]]; then
    echo "=================================================="
    echo "Generating summary workbook: $SUMMARY_XLSX"
    {
        echo "=================================================="
        echo "Generating summary workbook: $SUMMARY_XLSX"
    } >> "$RUNNER_LOG"

    if python3 "$SUMMARY_EXE" \
        --output "$SUMMARY_XLSX" \
        --log-dir "$OUTPUT_PATH" \
        --p "$P" \
        "${DATA_LIST[@]}" >> "$RUNNER_LOG" 2>&1; then
        echo "Summary workbook created: $SUMMARY_XLSX"
        {
            echo "Summary workbook created: $SUMMARY_XLSX"
        } >> "$RUNNER_LOG"
    else
        echo "Summary workbook generation failed."
        {
            echo "Summary workbook generation failed."
        } >> "$RUNNER_LOG"
        FAILED=$((FAILED + 1))
    fi
fi

echo "=================================================="
echo "Finished. Failed runs: $FAILED"
exit "$FAILED"
