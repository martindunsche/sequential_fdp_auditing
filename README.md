# Sequential Auditing for f-Differential Privacy

This repository provides the artifact for the submission **"Sequential Auditing for f-Differential Privacy"**.

## Overview

We present **APT (Advanced Privacy Testing)**, a sequential auditing framework for f-Differential Privacy that:
- Detects violations across the full tradeoff curve.
- Adaptively determines near-optimal sample sizes, avoiding the excessively large samples common in auditing studies.
- Supports both **whitebox** and **blackbox** settings.
- Can be executed in single-run frameworks for DP-SGD.

## Repository Structure

```
sequential_fdp_auditing/
├── src_log/                         # Core R source files (APT algorithm)
│   ├── functions_functional.R
│   ├── classifiers.R
│   ├── mechanisms.R
│   ├── KDE_estimator.R
│   └── sequential23_mmd.py          # MMD-surrogate utilities (Section 4.3)
├── experiments_log/
│   ├── gaussian.R                   # Blackbox Gaussian mechanism
│   ├── gaussian_parametric.R        # Whitebox Gaussian mechanism
│   ├── laplace.R                    # Blackbox Laplace mechanism
│   ├── dpsgd.R                      # Blackbox DP-SGD
│   ├── real_world.R                 # Real-world DP-SGD auditing
│   ├── gaussian_comp_other_paper.R  # Comparison with prior work (Gaussian)
│   ├── laplace_comp_other_paper.R   # Comparison with prior work (Laplace)
│   ├── nondp_gaussian_auditor_sweep.R   # Non-DP Gaussian sweep over epsilon
│   ├── nondp_laplace_auditor_sweep.R    # Non-DP Laplace sweep over epsilon
│   ├── compute_mmd_eps_star_sweep.py    # MMD diagnostic sweep
│   ├── run_all.sh                   # Runs all main experiments
│   ├── real_world.sh                # Runs real-world experiments
│   ├── plot.py                      # Power-curve plots
│   ├── plot_power_curve_wilson_ribbon.py  # Power-curve plots with Wilson 95% CI ribbons
│   ├── plot_nondp_gaussian_auditor_sweep.py
│   ├── real_world_plot.py
│   └── results/                     # Pre-computed experimental results
├── scripts/local_scripts/
│   ├── setup_r_environment.sh       # R environment setup (cluster-compatible)
│   └── setup_r_environment_pace.sh  # PACE cluster variant
└── sgd_exp_code/                    # Python code for DP-SGD score generation
    ├── scripts/gen_scores_DP_whitebox.py
    ├── src/
    └── requirements.txt
```

## Requirements

### R Dependencies

- R ≥ 4.0
- `KernSmooth`, `rmutil`, `parallel` (built-in)

Run `source scripts/local_scripts/setup_r_environment.sh` to set up the environment and install packages automatically.

### Python Dependencies

```bash
pip install pandas numpy matplotlib scipy
pip install -r sgd_exp_code/requirements.txt  # for DP-SGD score generation
```

## Running the Experiments

```bash
cd experiments_log
./run_all.sh          # main experiments (Gaussian, Laplace, DP-SGD)
./real_world.sh       # real-world DP-SGD auditing
```

The non-DP sweeps and MMD diagnostic can be run individually:

```bash
Rscript nondp_gaussian_auditor_sweep.R
Rscript nondp_laplace_auditor_sweep.R
python compute_mmd_eps_star_sweep.py
```

See `sgd_exp_code/README.md` for instructions on regenerating the DP-SGD audit scores.

## Pre-computed Results

Pre-computed results are in `experiments_log/results/`:

| Directory | Description |
|-----------|-------------|
| `gaussian/` | Blackbox Gaussian mechanism |
| `gaussian_parametric/` | Whitebox Gaussian mechanism |
| `laplace/` | Blackbox Laplace mechanism |
| `dpsgd/` | Blackbox DP-SGD |
| `real_world/` | Real-world DP-SGD auditing |
| `gaussian_comp_other_paper/` | Comparison with prior work (Gaussian) |
| `laplace_comp_other_paper/` | Comparison with prior work (Laplace) |
| `nondp_gaussian_auditor_sweep/` | Non-DP Gaussian sweep |
| `nondp_laplace_auditor_sweep/` | Non-DP Laplace sweep |
| `mmd_eps_star_sweep/` | MMD diagnostic sweep |
