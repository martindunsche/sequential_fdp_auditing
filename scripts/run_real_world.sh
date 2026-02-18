#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../experiments"

DATA_P="../data/in_scores_000200.csv"
DATA_Q="../data/out_scores_000200.csv"
CORES=250
TRIALS=250
BURN=50

MU_LIST=(1.35 1.4)

for mu in "${MU_LIST[@]}"; do
  echo "=== Running mu=${mu} ==="
  ./real_world.R "${mu}" \
    --dataP="${DATA_P}" \
    --dataQ="${DATA_Q}" \
    --cores="${CORES}" \
    --trials="${TRIALS}" \
    --burn="${BURN}"
done
