# Sequential f-DP Auditing

This repository contains experiments, simulation scripts, and utilities for sequential auditing of f-DP mechanisms. It supports parametric and non-parametric classifiers, sequential evaluations, and comparisons to baselines from other papers.

## Repository Structure

### experiments/
- dpsgd.R
- gaussian.R
- gaussian_parametric.R 
- gaussian_comp_other_paper.R: Gaussian comparison experiments
- laplace.R 
- laplace_comp_other_paper.R: Laplace comparison experiments
- real_world_sgd.R: audits using trained models (still to do)
- results/: generated outputs (stopping time CSV files and power curve PNGs in subfolders)
- run_all.sh: batch launcher for multiple configurations

### src/
- functions.R: core sequential audit implementation
- KDE_estimator.R: KDE-based density estimation and classifier fitting
- mechanisms.R: mechanism wrappers (Gaussian, Laplace, DPSGD, etc.)

### scripts/


## Running an Experiment

Example (Gaussian mechanism, mean shift mu = 0.01, NONDP variant):

```bash
cd experiments
./gaussian_comp_other_paper.R 0.01 NONDP
```

Outputs are written to:

```
experiments/results/<script_name>/
```

These include:
- stops_*.csv (stopping times and test statistics)
- *_power_curve.png (if later generated using plot.py)

## Supported Mechanisms

- Gaussian (gaussian*.R)
- Laplace (laplace*.R)
- DPSGD (dpsgd.R)

All mechanisms provide a function `Mechanism(x)` that returns a privatized output sample.

## Output Format

Each CSV contains one row per sequential trial:
- trial: trial index
- stopped_at_n: sample size when stopping occurred
- T1_: estimated type-I error
- T2_: estimated type-II error

## Requirements

- R version 4.0 or newer
- Packages: KernSmooth, parallel, stats, rmutil