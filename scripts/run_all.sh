#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/../experiments"

CORES=250
TRIALS=1000
BURN=50

for mu in 0.5 0.8 1.0; do
    ./laplace.R "$mu" \
    --cores="$CORES" \
    --trials="$TRIALS" \
    --burn="$BURN"
done

for tau in 5 7 10; do
    ./dpsgd.R "$tau" \
    --cores="$CORES" \
    --trials="$TRIALS" \
    --burn="$BURN"
done

for mu in 0.5 0.8 1.0; do
    ./gaussian.R "$mu" \
    --cores="$CORES" \
    --trials="$TRIALS" \
    --burn="$BURN"
done

for mu in 0.5 0.8 1.0; do
    ./gaussian_parametric.R "$mu" \
    --cores="$CORES" \
    --trials="$TRIALS" \
    --burn="$BURN"
done
