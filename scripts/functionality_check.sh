#!/usr/bin/env bash
# functionality_check.sh
# ---------------------------------------------------------------
# Smoke-tests every experiment script with tiny settings to verify
# that the restructured paths and imports all work correctly,
# without triggering heavy computation.
#
# Usage (from any directory):
#   bash scripts/functionality_check.sh
# ---------------------------------------------------------------
set -uo pipefail

# Always run R scripts from experiments/ so ../R/ and ../results/ resolve correctly
cd "$(dirname "$0")/../experiments"

# --- Minimal resource settings ---
CORES=1
TRIALS=1
BURN=10
QSIMS=100    # simulate_gaussian_sup_quantile sims (default 10000 — slow)
QKMAX=100    # simulate_gaussian_sup_quantile k_max (default 10000 — slow)

PASS=0
FAIL=0
TMPLOG=$(mktemp /tmp/fdp_check_XXXXXX.log)

# ---------------------------------------------------------------
run_check() {
  local label="$1"
  local script="$2"
  shift 2
  printf "  %-45s " "$label"
  if Rscript "$script" "$@" > "$TMPLOG" 2>&1; then
    echo "PASS"
    PASS=$((PASS + 1))
  else
    echo "FAIL"
    echo "    --- last 8 lines of output ---"
    tail -8 "$TMPLOG" | sed 's/^/    /'
    FAIL=$((FAIL + 1))
  fi
}

# ---------------------------------------------------------------
echo "=== Functionality Check ==="
echo "Settings: trials=$TRIALS, burn=$BURN, cores=$CORES"
echo ""

echo "--- Simulated mechanisms ---"

run_check "gaussian (mu=0.5)" gaussian.R 0.5 \
  --cores=$CORES --trials=$TRIALS --burn=$BURN \
  --qalpha_sims=$QSIMS --qalpha_kmax=$QKMAX

run_check "gaussian_parametric (mu=0.5)" gaussian_parametric.R 0.5 \
  --cores=$CORES --trials=$TRIALS --burn=$BURN \
  --qalpha_sims=$QSIMS --qalpha_kmax=$QKMAX

run_check "laplace (mu=0.5)" laplace.R 0.5 \
  --cores=$CORES --trials=$TRIALS --burn=$BURN \
  --qalpha_sims=$QSIMS --qalpha_kmax=$QKMAX

run_check "dpsgd (tau=5)" dpsgd.R 5 \
  --cores=$CORES --trials=$TRIALS --burn=$BURN \
  --qalpha_sims=$QSIMS --qalpha_kmax=$QKMAX

echo ""
echo "--- Comparison experiments ---"

run_check "gaussian_comp DP (mu=0.1)" gaussian_comp.R 0.1 DP \
  --cores=$CORES --trials=$TRIALS --burn=$BURN \
  --qalpha_sims=$QSIMS --qalpha_kmax=$QKMAX

run_check "gaussian_comp NonDP (mu=0.1)" gaussian_comp.R 0.1 NonDP \
  --cores=$CORES --trials=$TRIALS --burn=$BURN \
  --qalpha_sims=$QSIMS --qalpha_kmax=$QKMAX

run_check "laplace_comp DP (mu=0.1)" laplace_comp.R 0.1 DP \
  --cores=$CORES --trials=$TRIALS --burn=$BURN \
  --qalpha_sims=$QSIMS --qalpha_kmax=$QKMAX

run_check "laplace_comp NonDP (mu=0.1)" laplace_comp.R 0.1 NonDP \
  --cores=$CORES --trials=$TRIALS --burn=$BURN \
  --qalpha_sims=$QSIMS --qalpha_kmax=$QKMAX

echo ""
echo "--- Real-world (precomputed scores) ---"

run_check "real_world (mu=1.0)" real_world.R 1.0 \
  --dataP=../data/in_scores_000200.csv \
  --dataQ=../data/out_scores_000200.csv \
  --cores=$CORES --trials=$TRIALS --burn=$BURN \
  --max_iter=50 --qalpha_sims=200 --qalpha_kmax=200

echo ""
echo "--- Python plotting ---"

printf "  %-45s " "plot_power_curves.py"
if python3 ../plotting/plot_power_curves.py > "$TMPLOG" 2>&1; then
  echo "PASS"
  PASS=$((PASS + 1))
else
  echo "FAIL"
  tail -8 "$TMPLOG" | sed 's/^/    /'
  FAIL=$((FAIL + 1))
fi

rm -f "$TMPLOG"

echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
[ $FAIL -eq 0 ] && exit 0 || exit 1
