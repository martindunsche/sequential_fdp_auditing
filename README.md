# Sequential Auditing for f-Differential Privacy

This repository provides the artifact for the submission **"Sequential Auditing for f-Differential Privacy"**.

## Overview

We present **APT (Advanced Privacy Testing)**, a sequential auditing framework for f-Differential Privacy that:
- Detects violations across the full tradeoff curve.
- Adaptively determines near-optimal sample sizes, avoiding the excessively large samples common in auditing studies
- Supports both **whitebox** and **blackbox** settings
- Can be executed in single-run frameworks for DP-SGD

## Repository Structure

```
sequential_fdp_auditing/
├── R/                             # Core auditing library
│   ├── audit_engine.R             # APT algorithm and sequential test
│   ├── classifiers.R              # Classifier factories
│   ├── mechanisms.R               # Mechanism implementations
│   └── kde_estimator.R            # KDE-based score estimation
├── experiments/                   # Experiment runner scripts
│   ├── gaussian.R                 # Blackbox Gaussian mechanism
│   ├── gaussian_parametric.R      # Whitebox Gaussian mechanism
│   ├── laplace.R                  # Blackbox Laplace mechanism
│   ├── dpsgd.R                    # Blackbox DP-SGD
│   ├── gaussian_comp.R            # Comparison with prior work (Gaussian)
│   ├── laplace_comp.R             # Comparison with prior work (Laplace)
│   └── real_world.R               # Real-world DP-SGD auditing (Figure 5)
├── scripts/
│   ├── run_all.sh                 # Run all main experiments
│   ├── run_real_world.sh          # Run real-world experiment
│   └── functionality_check.sh     # Smoke-test all scripts (fast)
├── plotting/
│   └── plot_power_curves.py       # Generate power curve plots from results
├── data/                          # Pre-computed audit scores (CSV)
├── results/                       # Pre-computed experimental results
│   ├── gaussian/
│   ├── gaussian_parametric/
│   ├── laplace/
│   ├── dpsgd/
│   ├── real_world/
│   ├── gaussian_comp/             # Comparison with prior work
│   └── laplace_comp/              # Comparison with prior work
└── python/                        # Real-world DP-SGD score generation (Python)
    ├── scripts/
    │   └── gen_scores_DP_whitebox.py
    ├── src/
    └── requirements.txt
```

## Requirements

### R Dependencies (for main experiments)

The R scripts require:
- R >= 4.0
- `parallel` (built-in)
- `KernSmooth` (for kernel density estimation)
- `rmutil` (for Laplace mechanism experiments)

### Python Dependencies (for DP-SGD score generation)

See `python/requirements.txt`.

## Quick Start

Run a smoke-test to verify all scripts work correctly (completes in about a minute):

```bash
bash scripts/functionality_check.sh
```

Run all main experiments (heavy computation, requires many cores):

```bash
bash scripts/run_all.sh
```

Regenerate all plots from pre-computed results:

```bash
python3 plotting/plot_power_curves.py
```

## Running Individual Experiments

All experiment scripts are run from any directory and accept common flags:

```bash
# Positional argument: privacy parameter (mu or tau)
# --cores=N      number of parallel workers (default 250)
# --trials=N     number of independent trials (default 250)
# --burn=N       burn-in period (default 100)
# --qalpha_sims=N  quantile simulation count (default 10000)
# --qalpha_kmax=N  quantile simulation k_max (default 10000)

Rscript experiments/gaussian.R 0.5 --cores=8 --trials=250 --burn=50
Rscript experiments/dpsgd.R 5 --cores=8 --trials=250 --burn=50
Rscript experiments/real_world.R 1.0 \
    --dataP=data/in_scores_000200.csv \
    --dataQ=data/out_scores_000200.csv
```

## Pre-computed Results

Pre-computed experimental results are available in `results/`:
- `gaussian/`: Blackbox Gaussian mechanism results
- `gaussian_parametric/`: Whitebox Gaussian mechanism results
- `laplace/`: Laplace mechanism results
- `dpsgd/`: DP-SGD results
- `real_world/`: Real-world DP-SGD auditing results
- `gaussian_comp/`: Comparison with prior work (Gaussian)
- `laplace_comp/`: Comparison with prior work (Laplace)
