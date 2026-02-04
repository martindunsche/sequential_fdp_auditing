# Sequential Auditing for f-Differential Privacy

This repository provides the artifact for the submission **"Sequential Auditing for f-Differential Privacy"**.

## Overview

We present **APT (Advanced Privacy Testing)**, a sequential auditing framework for f-Differential Privacy that:
- Detects violations across the ful tradeoff curve.
- Adaptively determines near-optimal sample sizes, avoiding the excessively large samples common in auditing studies
- Supports both **whitebox** and **blackbox** settings
- Can be executed in single-run frameworks for DP-SGD

## Repository Structure

```
sequential_fdp_auditing/
├── experiments_log/           
│   ├── gaussian.R             # Blackbox Gaussian mechanism experiments
│   ├── gaussian_parametric.R  # Whitebox Gaussian mechanism experiments
│   ├── laplace.R              # Blackbox Laplace mechanism experiments
│   ├── dpsgd.R                # Blackbox DP-SGD experiments
│   ├── real_world.R           # Real-world DP-SGD auditing (Figure 5)
│   ├── run_all.sh             # Main script to run all experiments
│   ├── plot.py                
│   └── results/               # Pre-computed experimental results
├── src_log/                   
│   ├── functions_functional.R # Main auditing functions (APT algorithm)
│   ├── classifiers.R          
│   ├── mechanisms.R           
│   └── KDE_estimator.R        
└── sgd_exp_code/              # Real-world DP-SGD auditing source code (Python)
    ├── scripts/               
    │   └── gen_scores_DP_whitebox.py
    ├── src/                   
    └── requirements.txt       
```

## Requirements

### R Dependencies (for main experiments)

The R scripts require:
- R ≥ 4.0
- `parallel` (built-in)
- `KernSmooth` (for kernel density estimation)

### Python Dependencies (for DP-SGD score generation)

See `sgd_exp_code/requirements.txt`. 

## Pre-computed Results

Pre-computed experimental results are available in `experiments_log/results/`:
- `gaussian/`: Blackbox Gaussian mechanism results
- `gaussian_parametric/`: Whitebox Gaussian mechanism results  
- `laplace/`: Laplace mechanism results
- `dpsgd/`: DP-SGD results
- `real_world/`: Real-world DP-SGD auditing results
- `gaussian_comp_other_paper/`: Comparison with prior work
- `laplace_comp_other_paper/`: Comparison with prior work


