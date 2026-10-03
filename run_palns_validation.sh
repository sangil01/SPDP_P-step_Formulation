#!/usr/bin/env bash
# PALNS incumbent validation runner.
# Every solver parameter and the instance list are read from run_release_many.sh,
# so only the incumbent backend, its time limit and the stop point differ.
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
export GUROBI_HOME="/opt/gurobi1303/linux64"
export PATH="$GUROBI_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$GUROBI_HOME/lib"
EXE="$SCRIPT_DIR/SPDP_codes/build-release/SPDP_P_step"

eval "$(sed -n '/^# INPUT PARAMETERS/,/^)$/p' "$SCRIPT_DIR/run_release_many.sh")"

INITIAL_INCUMBENT_BACKEND=palns
INITIAL_INCUMBENT_TIME_LIMIT="${PALNS_TIME_LIMIT:-1200}"
OUTPUT_DIR="${PALNS_OUTPUT_DIR:-SPDP_output_palns}"
STOP_AFTER="${STOP_AFTER:-mip-start}"          # none | incumbent | mip-start
JOBS="${JOBS:-4}"
BASELINE_DIR="${BASELINE_DIR:-SPDP_output_TF2_capacity_cut_travel_time_cut_ver2}"
TARGET_MODE="${TARGET_MODE:-1}"                # 1: stop at the baseline incumbent cost
PALNS_SEED="${PALNS_SEED:-$INITIAL_INCUMBENT_PALNS_SEED}"
PALNS_THREADS="${PALNS_THREADS:-$INITIAL_INCUMBENT_PALNS_THREADS}"
PALNS_EXTRA="${PALNS_EXTRA:-}"
if [[ -n "${INSTANCES:-}" ]]; then
    read -r -a DATA_LIST <<< "$INSTANCES"
fi
mkdir -p "$OUTPUT_DIR"

run_one() {
    local data_name="$1"
    local stem="${data_name%.dat}"
    local inst_log="$OUTPUT_DIR/${stem}_runner.log"
    local PALNS_TARGET_ARGS=""
    if [[ "$TARGET_MODE" == "1" ]]; then
        local base_log="$BASELINE_DIR/${stem}_p0_log.txt"
        if [[ -f "$base_log" ]]; then
            local target
            target="$(grep -m1 '^\[initial-incumbent\] original_cost=' "$base_log" | sed 's/.*=//')"
            if [[ -n "$target" ]]; then
                PALNS_TARGET_ARGS="--initial-incumbent-palns-target-cost $target"
            fi
        fi
    fi
    "$EXE" "$data_name" \
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
        --initial-incumbent-palns-composite-kappa "$INITIAL_INCUMBENT_PALNS_COMPOSITE_KAPPA" \
        --initial-incumbent-palns-composite-skeletons "$INITIAL_INCUMBENT_PALNS_COMPOSITE_SKELETONS" \
        --initial-incumbent-palns-sim-weights "$INITIAL_INCUMBENT_PALNS_SIM_WEIGHTS" \
        --initial-incumbent-palns-construct-starts "$INITIAL_INCUMBENT_PALNS_CONSTRUCT_STARTS" \
        --initial-incumbent-palns-construct-composite-filter "$INITIAL_INCUMBENT_PALNS_CONSTRUCT_COMPOSITE_FILTER" \
        --initial-incumbent-palns-construct-order "$INITIAL_INCUMBENT_PALNS_CONSTRUCT_ORDER" \
        --initial-incumbent-palns-eta "$INITIAL_INCUMBENT_PALNS_ETA" \
        --initial-incumbent-palns-kmax-ratio "$INITIAL_INCUMBENT_PALNS_KMAX_RATIO" \
        --initial-incumbent-palns-kmax-min "$INITIAL_INCUMBENT_PALNS_KMAX_MIN" \
        --initial-incumbent-palns-ls-frequency "$INITIAL_INCUMBENT_PALNS_LS_FREQUENCY" \
        --initial-incumbent-palns-ls-every-m "$INITIAL_INCUMBENT_PALNS_LS_EVERY_M" \
        --initial-incumbent-palns-ls-first-improvement "$INITIAL_INCUMBENT_PALNS_LS_FIRST_IMPROVEMENT" \
        --initial-incumbent-palns-ls-deep-ratio "$INITIAL_INCUMBENT_PALNS_LS_DEEP_RATIO" \
        --initial-incumbent-palns-phase2-enable "$INITIAL_INCUMBENT_PALNS_PHASE2_ENABLE" \
        --initial-incumbent-palns-phase2-max-noimp "$INITIAL_INCUMBENT_PALNS_PHASE2_MAX_NOIMP" \
        --initial-incumbent-palns-sa-start-ratio "$INITIAL_INCUMBENT_PALNS_SA_START_RATIO" \
        --initial-incumbent-palns-sa-start-ratio-phase1 "$INITIAL_INCUMBENT_PALNS_SA_START_RATIO_PHASE1" \
        --initial-incumbent-palns-sa-end-ratio "$INITIAL_INCUMBENT_PALNS_SA_END_RATIO" \
        --initial-incumbent-palns-sa-cycle-iters "$INITIAL_INCUMBENT_PALNS_SA_CYCLE_ITERS" \
        --initial-incumbent-palns-sa-cycle-iters-phase2 "$INITIAL_INCUMBENT_PALNS_SA_CYCLE_ITERS_PHASE2" \
        --initial-incumbent-palns-sa-restart-best "$INITIAL_INCUMBENT_PALNS_SA_RESTART_BEST" \
        --initial-incumbent-palns-phase1-eta-max "$INITIAL_INCUMBENT_PALNS_PHASE1_ETA_MAX" \
        --initial-incumbent-palns-portfolio-mixed-cycles "$INITIAL_INCUMBENT_PALNS_PORTFOLIO_MIXED_CYCLES" \
        --initial-incumbent-palns-portfolio-diverse "$INITIAL_INCUMBENT_PALNS_PORTFOLIO_DIVERSE" \
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
        --pre-mip-separate-main-mip-root "$PRE_MIP_SEPARATE_MAIN_MIP_ROOT" \
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
        --full-enumeration-rc-detail-log "$FULL_ENUMERATION_RC_DETAIL_LOG" \
        --stop-after-initial-incumbent "$STOP_AFTER" \
        --initial-incumbent-palns-seed "$PALNS_SEED" \
        --initial-incumbent-palns-threads "$PALNS_THREADS" \
        $PALNS_TARGET_ARGS $PALNS_EXTRA > "$inst_log" 2>&1
    echo "done $data_name exit=$?"
}

running=0
for data_name in "${DATA_LIST[@]}"; do
    run_one "$data_name" &
    running=$((running + 1))
    if (( running >= JOBS )); then
        wait -n
        running=$((running - 1))
    fi
done
wait
