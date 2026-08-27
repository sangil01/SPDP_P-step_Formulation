#!/usr/bin/env bash
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/SPDP_codes/build-release/SPDP_P_step"
SUMMARY_EXE="$SCRIPT_DIR/generate_root_cg_summary_xlsx.py"

if [[ ! -x "$EXE" ]]; then
    echo "Release executable not found or not executable: $EXE"
    exit 1
fi

# ============================================
# INPUT PARAMETERS
# ============================================
P=0 # 0 uses the direct binary two-index IP when SOLVER_MODE=enumeration
SOLVER_MODE=enumeration # Options: enumeration, branch-and-price
ENUMERATION_SOS1_MODE=default # Options: default, sos1-auto, sos1-native (used only with SOLVER_MODE=enumeration)
ENUMERATION_OBJECTIVE=original-cost  # Options: original-cost, travel-cost-only, duration, duration-plus-fixed (used only with SOLVER_MODE=enumeration)
BNP_TREE_MODE=root-only # Options: root-only, full-tree
NODE_CG_PHASE1_MODE=heuristic-cg # Options: exact-cg, heuristic-cg, heuristic-cg-3-step
NODE_CG_PHASE2_PRICING_MODE=full-enumeration # Options: exact-pricing, heuristic-pricing-then-exact, full-enumeration
SOLVER_TIME_LIMIT=60 # seconds
INITIAL_INCUMBENT_ENABLE=1
INITIAL_INCUMBENT_TIME_LIMIT=1200 # seconds; 0 means no time limit
INITIAL_INCUMBENT_MAX_K_INCREMENTS=2 # additional K values after the initial lower bound
INITIAL_INCUMBENT_TIMEOUT_ACTION=advance # Options: stop, advance
INITIAL_INCUMBENT_BACKEND=two-index-milp # Options: two-index-milp, cp-sat
INITIAL_INCUMBENT_CP_WORKERS=0
INITIAL_INCUMBENT_CP_RANDOM_SEED=1
INITIAL_INCUMBENT_CP_LOG_PROGRESS=0
INITIAL_INCUMBENT_CP_TERMINAL_BALANCE=1
INITIAL_INCUMBENT_CP_FULL_RESERVOIR=1
INITIAL_INCUMBENT_CP_CONTAINER_WORKLOAD=0
INITIAL_INCUMBENT_CP_AGGREGATE_DURATION=1
INITIAL_INCUMBENT_CP_FIRST_PICKUP_SYMMETRY=0
INITIAL_INCUMBENT_CP_SYMMETRY_43=0
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
ADD_VI_35=0
ADD_VI_36_COMBINED=0
VI_36_SUBSET_MAX_SIZE=0 # 0이면 비활성화, 1이면 기존 singleton Eq.(36)과 동일, 7이면 현재 instance들에 대해 전체 set까지 포함
ADD_VI_Request_BLOCK_SEC=0
VI_Request_BLOCK_SEC_MAX_SIZE=0 # 0이면 비활성화, 활성화할 때는 2 이상
ADD_VI_44=1
VI_44_K_MIN_USE_COR=1
VI_44_K_MIN_USE_SUBPROBLEM=1
VI_44_K_MIN_USE_VEHICLE_ASSIGNMENT=0
VI_44_SUBPROBLEM_TYPE=ip # Options: lp, ip
VI_44_SUBPROBLEM_TIME_LIMIT=600 # seconds; 0 means no time limit
VI_44_SUBPROBLEM_ADD_TIME_CONSTRAINTS=1 # 1 adds state-time B variables and route-duration constraints
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
    '''"RecDep_day_A1.dat"
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
    "RecDep_day_C4.dat"'''
    #=================================#
    #=========Request 50 이하=========#
    '''"RecDep_day_A12.dat"
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
    "RecDep_day_D2.dat"'''
    #=================================#
    #=====Request 50 초과 100 이하=====#
    '''"RecDep_day_B15.dat"
    "RecDep_day_B16.dat"
    "RecDep_day_B17.dat"
    "RecDep_day_B18.dat"
    "RecDep_day_B19.dat"
    "RecDep_day_B20.dat"
    "RecDep_day_C13.dat"
    "RecDep_day_C14.dat"
    "RecDep_day_C15.dat"'''
    "RecDep_day_C16.dat"
    "RecDep_day_C17.dat"
    "RecDep_day_C18.dat"
    "RecDep_day_C19.dat"
    "RecDep_day_C20.dat"
    '''"RecDep_day_D3.dat"
    "RecDep_day_D4.dat"'''
    "RecDep_day_D5.dat"
    "RecDep_day_D6.dat"
    #"RecDep_day_D7.dat"
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
    ENUMERATION_SOS1_TAG="_${ENUMERATION_SOS1_MODE}_${ENUMERATION_OBJECTIVE}"
fi
RUN_TAG="P${P}_${SOLVER_MODE}_${BNP_TREE_MODE}_${SOLVER_TIME_LIMIT}s_${NODE_CG_PHASE1_MODE}_${NODE_CG_PHASE2_PRICING_MODE}${ENUMERATION_SOS1_TAG}${RUN_TAG_SUFFIX}"
if [[ "$INITIAL_INCUMBENT_ENABLE" == "1" ]]; then
    RUN_TAG+="_initial-incumbent-${INITIAL_INCUMBENT_BACKEND}-${INITIAL_INCUMBENT_TIME_LIMIT}s-kinc${INITIAL_INCUMBENT_MAX_K_INCREMENTS}-timeout-${INITIAL_INCUMBENT_TIMEOUT_ACTION}"
    if [[ "$INITIAL_INCUMBENT_BACKEND" == "cp-sat" ]]; then
        RUN_TAG+="-cpw${INITIAL_INCUMBENT_CP_WORKERS}-seed${INITIAL_INCUMBENT_CP_RANDOM_SEED}-tb${INITIAL_INCUMBENT_CP_TERMINAL_BALANCE}-fr${INITIAL_INCUMBENT_CP_FULL_RESERVOIR}-cw${INITIAL_INCUMBENT_CP_CONTAINER_WORKLOAD}-ad${INITIAL_INCUMBENT_CP_AGGREGATE_DURATION}-fp${INITIAL_INCUMBENT_CP_FIRST_PICKUP_SYMMETRY}-s43${INITIAL_INCUMBENT_CP_SYMMETRY_43}"
    fi
fi
RUNNER_LOG="$SCRIPT_DIR/${RUN_TAG}.log"
SUMMARY_XLSX="$SCRIPT_DIR/${RUN_TAG}.xlsx"

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
        echo "solver-mode: $SOLVER_MODE"
        echo "enumeration-sos1-mode: $ENUMERATION_SOS1_MODE"
        echo "enumeration-objective: $ENUMERATION_OBJECTIVE"
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
        echo "initial-incumbent-cp-workers: $INITIAL_INCUMBENT_CP_WORKERS"
        echo "initial-incumbent-cp-random-seed: $INITIAL_INCUMBENT_CP_RANDOM_SEED"
        echo "initial-incumbent-cp-log-progress: $INITIAL_INCUMBENT_CP_LOG_PROGRESS"
        echo "initial-incumbent-cp-terminal-balance: $INITIAL_INCUMBENT_CP_TERMINAL_BALANCE"
        echo "initial-incumbent-cp-full-reservoir: $INITIAL_INCUMBENT_CP_FULL_RESERVOIR"
        echo "initial-incumbent-cp-container-workload: $INITIAL_INCUMBENT_CP_CONTAINER_WORKLOAD"
        echo "initial-incumbent-cp-aggregate-duration: $INITIAL_INCUMBENT_CP_AGGREGATE_DURATION"
        echo "initial-incumbent-cp-first-pickup-symmetry: $INITIAL_INCUMBENT_CP_FIRST_PICKUP_SYMMETRY"
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
        echo "add-vi-35: $ADD_VI_35"
        echo "add-vi-36-combined: $ADD_VI_36_COMBINED"
        echo "vi-36-subset-max-size: $VI_36_SUBSET_MAX_SIZE"
        echo "add-vi-request-block-sec: $ADD_VI_Request_BLOCK_SEC"
        echo "vi-request-block-sec-max-size: $VI_Request_BLOCK_SEC_MAX_SIZE"
        echo "add-vi-44: $ADD_VI_44"
        echo "vi-44-k-min-use-cor: $VI_44_K_MIN_USE_COR"
        echo "vi-44-k-min-use-subproblem: $VI_44_K_MIN_USE_SUBPROBLEM"
        echo "vi-44-k-min-use-vehicle-assignment: $VI_44_K_MIN_USE_VEHICLE_ASSIGNMENT"
        echo "vi-44-subproblem-type: $VI_44_SUBPROBLEM_TYPE"
        echo "vi-44-subproblem-time-limit: $VI_44_SUBPROBLEM_TIME_LIMIT"
        echo "vi-44-subproblem-add-time-constraints: $VI_44_SUBPROBLEM_ADD_TIME_CONSTRAINTS"
        echo "vi-44-vehicle-assignment-add-tsp-bound: $VI_44_VEHICLE_ASSIGNMENT_ADD_TSP_BOUND"
        echo "vi-44-vehicle-assignment-add-container-bound: $VI_44_VEHICLE_ASSIGNMENT_ADD_CONTAINER_BOUND"
        echo "vi-44-vehicle-assignment-time-limit: $VI_44_VEHICLE_ASSIGNMENT_TIME_LIMIT"
        echo "----------------------------------------"
    } >> "$RUNNER_LOG"

    if "$EXE" "$data_name" \
        --p "$P" \
        --solver-mode "$SOLVER_MODE" \
        --enumeration-sos1-mode "$ENUMERATION_SOS1_MODE" \
        --enumeration-objective "$ENUMERATION_OBJECTIVE" \
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
        --initial-incumbent-cp-workers "$INITIAL_INCUMBENT_CP_WORKERS" \
        --initial-incumbent-cp-random-seed "$INITIAL_INCUMBENT_CP_RANDOM_SEED" \
        --initial-incumbent-cp-log-progress "$INITIAL_INCUMBENT_CP_LOG_PROGRESS" \
        --initial-incumbent-cp-terminal-balance "$INITIAL_INCUMBENT_CP_TERMINAL_BALANCE" \
        --initial-incumbent-cp-full-reservoir "$INITIAL_INCUMBENT_CP_FULL_RESERVOIR" \
        --initial-incumbent-cp-container-workload "$INITIAL_INCUMBENT_CP_CONTAINER_WORKLOAD" \
        --initial-incumbent-cp-aggregate-duration "$INITIAL_INCUMBENT_CP_AGGREGATE_DURATION" \
        --initial-incumbent-cp-first-pickup-symmetry "$INITIAL_INCUMBENT_CP_FIRST_PICKUP_SYMMETRY" \
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
        --prune-pickup-symmetry-43 "$PRUNE_PICKUP_SYMMETRY_43" \
        --prune-delivery-symmetry-43 "$PRUNE_DELIVERY_SYMMETRY_43" \
        --add-vi-35 "$ADD_VI_35" \
        --add-vi-36-combined "$ADD_VI_36_COMBINED" \
        --vi-36-subset-max-size "$VI_36_SUBSET_MAX_SIZE" \
        --add-vi-request-block-sec "$ADD_VI_Request_BLOCK_SEC" \
        --vi-request-block-sec-max-size "$VI_Request_BLOCK_SEC_MAX_SIZE" \
        --add-vi-44 "$ADD_VI_44" \
        --vi-44-k-min-use-cor "$VI_44_K_MIN_USE_COR" \
        --vi-44-k-min-use-subproblem "$VI_44_K_MIN_USE_SUBPROBLEM" \
        --vi-44-k-min-use-vehicle-assignment "$VI_44_K_MIN_USE_VEHICLE_ASSIGNMENT" \
        --vi-44-subproblem-type "$VI_44_SUBPROBLEM_TYPE" \
        --vi-44-subproblem-time-limit "$VI_44_SUBPROBLEM_TIME_LIMIT" \
        --vi-44-subproblem-add-time-constraints "$VI_44_SUBPROBLEM_ADD_TIME_CONSTRAINTS" \
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
        --log-dir "$SCRIPT_DIR/SPDP_output" \
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
