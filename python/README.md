# SGD Experiment Code for Sequential DP Auditing

This directory contains the code for generating audit scores for sequential differential privacy (DP) auditing of DP-SGD, as described in **Section 4.4** of the paper.

## Overview

The main script `scripts/gen_scores_DP_whitebox.py` implements white-box auditing of DP-SGD using **Dirac canaries**. It trains a differentially private model while simultaneously generating the **in-scores** and **out-scores** needed for sequential auditing, where 

- **In-scores**: Scores computed when the canary is included in training (sampled with probability `q` at each step).
- **Out-scores**: Scores computed for canaries that were never included (sampled from a Gaussian null distribution).

## Requirements

The code is written in python 3.10. Install the required dependencies:

```bash
pip install -r requirements.txt
```

## Usage

Run the score generation script with the following (optional) parameters:

```bash
python scripts/gen_scores_DP_whitebox.py \
    --logical-batch-size 4096 \
    --max-physical-batch-size 128 \
    --aug-multiplicity 1 \
    --max-grad-norm 1.0 \
    --epsilon 8.0 \
    --delta 1e-5 \
    --epochs 200 \
    --lr 4.0 \
    --momentum 0.0 \
    --noise-multiplier 3.0 \
    --ckpt-interval 20 \
    --canary-count 10000 \
    --pkeep 0.5 \
    --database-seed <128-bit-integer> \
    --data-dir ./data \
    --log-dir ./logs \
    --num-workers 4
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--logical-batch-size` | 4096 | Logical batch size for gradient accumulation |
| `--max-physical-batch-size` | 128 | Maximum physical batch size (for memory efficiency) |
| `--max-grad-norm` | 1.0 | Per-sample gradient clipping threshold |
| `--noise-multiplier` | 3.0 | Noise multiplier σ for DP-SGD |
| `--epsilon` | 8.0 | Target privacy budget ε |
| `--delta` | 1e-5 | Privacy parameter δ |
| `--epochs` | 200 | Number of training epochs |
| `--canary-count` | 10000 | Number of Dirac canaries to use |
| `--pkeep` | 0.5 | Probability of including each canary in training |
| `--ckpt-interval` | 20 | Checkpoint and score save interval (epochs) |
| `--database-seed` | Random | 128-bit seed for reproducibility |

## Output

The script creates an experiment directory with the following structure:

```
data/mislabeled-canaries-<seed>-<canary_count>-<pkeep>-cifar10/
├── hparams.json                    # Hyperparameters used
├── ckpt/                           # Model checkpoints
│   └── ckpt_epoch_<N>.pt
├── in_scores_<epoch>.csv           # In-scores at each checkpoint
├── out_scores_<epoch>.csv          # Out-scores at each checkpoint
└── privacy_params_<epoch>.csv      # Privacy parameters (ε, δ) at each checkpoint
```

### Score Files

- **`in_scores_<epoch>.csv`**: Average scores for canaries that were included in training. Shape: `(canary_count,)`
- **`out_scores_<epoch>.csv`**: Scores for canaries that were never included (sampled from null distribution). Shape: `(canary_count,)`
- **`privacy_params_<epoch>.csv`**: Current privacy budget (ε, δ) at the checkpoint.

## Resume Training

The script automatically resumes from the latest checkpoint if one exists in the experiment directory. It restores:
- Model and optimizer state
- Privacy accountant history
- Accumulated scores

## Example

Generate scores with a specific seed for reproducibility:

```bash
python scripts/gen_scores_DP_whitebox.py \
    --database-seed 27198899012190525004019618245709479116 \
    --canary-count 10000 \
    --pkeep 0.5 \
    --epochs 200 \
    --data-dir ./data
```

This will create output in:
```
./data/mislabeled-canaries-27198899012190525004019618245709479116-10000-0.5-cifar10/
```

## Notes

- **GPU**: The script automatically uses CUDA if available.
- **Memory**: Use `--max-physical-batch-size` to control GPU memory usage with gradient accumulation.
- **Reproducibility**: Provide `--database-seed` for deterministic results.
