#!/usr/bin/env bash
set -euo pipefail

# Inputs (edit if needed)
SCRIPT="./real_world.R"
DATA_P="in_scores_000200.csv"
DATA_Q="out_scores_000200.csv"
CORES=250
TRIALS=250
BURN=50

# R: seq(0.5, 1.3, length.out=5)
MU_LIST=(1.35 1.4)

for mu in "${MU_LIST[@]}"; do
  echo "=== Running mu=${mu} ==="
  "${SCRIPT}" "${mu}" \
    --dataP="${DATA_P}" \
    --dataQ="${DATA_Q}" \
    --cores="${CORES}" \
    --trials="${TRIALS}" \
    --burn="${BURN}"
done
