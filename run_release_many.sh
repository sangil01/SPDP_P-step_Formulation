#!/usr/bin/env bash
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/SPDP_codes/build-release/SPDP_P_step"

if [[ ! -x "$EXE" ]]; then
    echo "Release executable not found or not executable: $EXE"
    exit 1
fi

# ============================================
# INPUT PARAMETERS
# ============================================
P=2
SOLVER_TIME_LIMIT=1800
DUMP_PSTEPS=10
VALIDATE_PSTEPS=0
SOLVE_MODEL=1
PRUNE_INFEASIBLE_EDGES=1
PRUNE_DOMINATED_EDGES=1
PRUNE_SYMMETRY_40=1
PRUNE_SYMMETRY_41=1
PRUNE_SYMMETRY_43=1
ADD_VI_35=0
ADD_VI_36=0
ADD_VI_44=0
DATA_LIST=(
    "A0.dat"
    "A1.dat"
    "A2.dat"
    "A3.dat"
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

FAILED=0

for data_name in "${DATA_LIST[@]}"; do
    echo "=================================================="
    echo "Running data: $data_name"

    "$EXE" "$data_name" \
        --p "$P" \
        --solver-time-limit "$SOLVER_TIME_LIMIT" \
        --dump-psteps "$DUMP_PSTEPS" \
        --validate-psteps "$VALIDATE_PSTEPS" \
        --solve "$SOLVE_MODEL" \
        --prune-infeasible-edges "$PRUNE_INFEASIBLE_EDGES" \
        --prune-dominated-edges "$PRUNE_DOMINATED_EDGES" \
        --prune-symmetry-40 "$PRUNE_SYMMETRY_40" \
        --prune-symmetry-41 "$PRUNE_SYMMETRY_41" \
        --prune-symmetry-43 "$PRUNE_SYMMETRY_43" \
        --add-vi-35 "$ADD_VI_35" \
        --add-vi-36 "$ADD_VI_36" \
        --add-vi-44 "$ADD_VI_44"

    exit_code=$?
    if [[ $exit_code -ne 0 ]]; then
        echo "Failed: $data_name"
        FAILED=$((FAILED + 1))
    else
        echo "Completed: $data_name"
    fi
done

echo "=================================================="
echo "Finished. Failed runs: $FAILED"
exit "$FAILED"
