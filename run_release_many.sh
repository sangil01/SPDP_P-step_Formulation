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
P=3
SOLVER_MODE=root-cg # Options: enumeration, root-cg
NODE_CG_PHASE1_MODE=heuristic-cg # Options: exact-cg, heuristic-cg
NODE_CG_PHASE2_PRICING_MODE=heuristic-pricing-then-exact # Options: exact-pricing, heuristic-pricing-then-exact
SOLVER_TIME_LIMIT=120
GUROBI_THREADS=0
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
ADD_VI_36=0
ADD_VI_44=0
CG_MAX_ITERATIONS_PER_PHASE=1000
EXACT_PRICING_MAX_COLUMNS_PER_START=64
EXACT_PRICING_MAX_TOTAL_COLUMNS_PER_ROUND=4096
CG_REDUCED_COST_TOLERANCE=-1e-6
PHASE2_HEURISTIC_MAX_STARTS=4096  #1024면 모든 인스턴스에서 |x,\sigma(s)|보다 큼 (A11에서 701개임)
PHASE2_HEURISTIC_START_RATIO=1.0
PHASE2_HEURISTIC_LADDER_LEVELS=0 #0이면 모든 starte (node,state) 조합을 탐색. (위에 MAX_STARTS와 상관없이)
PHASE2_HEURISTIC_MAX_COLUMNS_PER_START=8
PHASE2_HEURISTIC_MAX_COLUMNS_TOTAL=512
PHASE2_HEURISTIC_SEARCH_COLUMN_RATIO=2.0
PHASE2_HEURISTIC_START_SCORE_MODE=one-step-min # Options: one-step-min
PHASE2_HEURISTIC_ENGINE=labeling # Options: labeling, shallow-search
PHASE2_LABELING_TOP_K_NEXT=4
PHASE2_SHALLOW_K1=4
PHASE2_SHALLOW_K2=2
PHASE2_COLUMN_POOL_ENABLE=1
PHASE2_COLUMN_POOL_MAX_SIZE=10000
PHASE2_COLUMN_POOL_MAX_REPRICE=1024
DATA_LIST=(
    "RecDep_day_B1.dat"
    "RecDep_day_B2.dat"
    "RecDep_day_C1.dat"
    "RecDep_day_C2.dat"
    "RecDep_day_C3.dat"
    "RecDep_day_C4.dat"
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
)
# Put one data file name per line in DATA_LIST.
# ============================================

RUNNER_LOG="$SCRIPT_DIR/P${P}_root-cg_${SOLVER_TIME_LIMIT}s_${NODE_CG_PHASE1_MODE}_${NODE_CG_PHASE2_PRICING_MODE}.log"
SUMMARY_XLSX="$SCRIPT_DIR/P${P}_root-cg_${SOLVER_TIME_LIMIT}s_${NODE_CG_PHASE1_MODE}_${NODE_CG_PHASE2_PRICING_MODE}.xlsx"

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
        echo "node-cg-phase1-mode: $NODE_CG_PHASE1_MODE"
        echo "node-cg-phase2-pricing-mode: $NODE_CG_PHASE2_PRICING_MODE"
        echo "vi-formulation: $VI_FORMULATION"
        echo "cg-max-iterations-per-phase: $CG_MAX_ITERATIONS_PER_PHASE"
        echo "exact-pricing-max-columns-per-start: $EXACT_PRICING_MAX_COLUMNS_PER_START"
        echo "exact-pricing-max-total-columns-per-round: $EXACT_PRICING_MAX_TOTAL_COLUMNS_PER_ROUND"
        echo "cg-reduced-cost-tolerance: $CG_REDUCED_COST_TOLERANCE"
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
        echo "prune-pickup-symmetry-43: $PRUNE_PICKUP_SYMMETRY_43"
        echo "prune-delivery-symmetry-43: $PRUNE_DELIVERY_SYMMETRY_43"
        echo "----------------------------------------"
    } >> "$RUNNER_LOG"

    if "$EXE" "$data_name" \
        --p "$P" \
        --solver-mode "$SOLVER_MODE" \
        --node-cg-phase1-mode "$NODE_CG_PHASE1_MODE" \
        --node-cg-phase2-pricing-mode "$NODE_CG_PHASE2_PRICING_MODE" \
        --solver-time-limit "$SOLVER_TIME_LIMIT" \
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
        --add-vi-36 "$ADD_VI_36" \
        --add-vi-44 "$ADD_VI_44" \
        --cg-max-iterations-per-phase "$CG_MAX_ITERATIONS_PER_PHASE" \
        --exact-pricing-max-columns-per-start "$EXACT_PRICING_MAX_COLUMNS_PER_START" \
        --exact-pricing-max-total-columns-per-round "$EXACT_PRICING_MAX_TOTAL_COLUMNS_PER_ROUND" \
        --cg-reduced-cost-tolerance "$CG_REDUCED_COST_TOLERANCE" \
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
        --phase2-column-pool-max-reprice "$PHASE2_COLUMN_POOL_MAX_REPRICE" >> "$RUNNER_LOG" 2>&1; then
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

if [[ "$SOLVER_MODE" == "root-cg" ]]; then
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
