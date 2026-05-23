#!/usr/bin/env python3
"""Compute epsilon^* and plot the MMD surrogate margin curve for each mechanism.

For each mechanism, estimates the population MMD(P, Q) under the chosen RBF
kernel, sweeps the claimed epsilon grid, and plots:

  (a) MMD surrogate margin  Delta_MMD(eps) = MMD_hat(P,Q) - tau_MMD(eps, delta)
  (b) Decomposed components: MMD_hat (flat) vs tau_MMD(eps, delta) (curve)

The crossover epsilon^* is where Delta_MMD = 0; for eps > eps^* the empirical
MMD-surrogate margin is nonpositive.

Example:
  python experiments_log/mmd_surrogate_diagnostic/compute_mmd_eps_star.py \\
    --mechanisms NonDPGaussian1 NonDPLaplace1 \\
    --mechanism-epsilon 0.1 \\
    --eps-grid 0.01:3.0:80 \\
    --mmd-samples 5000 \\
    --outdir experiments_log/results/mmd_surrogate_diagnostic
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src_log"))
from sequential23_mmd import (
    eps_star_from_mmd,
    estimate_mmd,
    median_bandwidth,
    sample_mechanism_outputs,
    tau_mmd,
)


MARGIN_YLIM = (-1.1, 0.3)


def parse_eps_grid(spec):
    spec = spec.strip()
    if ":" in spec:
        parts = spec.split(":")
        if len(parts) != 3:
            raise ValueError("--eps-grid must be start:stop:num or comma-separated values.")
        return np.linspace(float(parts[0]), float(parts[1]), int(parts[2]))
    return np.array([float(x) for x in spec.split(",") if x.strip()], dtype=float)


def plot_mechanism(mechanism, mmd_hat, eps_star, eps_grid, tau_vals, margin_vals, outdir):
    fig, ax = plt.subplots(figsize=(5.0, 4.0))

    ax.axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    ax.plot(eps_grid, margin_vals, color="steelblue", linewidth=1.8)
    ax.fill_between(
        eps_grid, margin_vals, 0,
        where=np.array(margin_vals) < 0,
        alpha=0.15, color="tomato",
    )
    ax.fill_between(
        eps_grid, margin_vals, 0,
        where=np.array(margin_vals) >= 0,
        alpha=0.10, color="steelblue",
    )

    if np.isfinite(eps_star):
        ax.axvline(eps_star, color="tomato", linewidth=1.0, linestyle=":")
        ax.annotate(
            r"$\widehat{\varepsilon}_{\mathrm{MMD}}^\star"
            rf" \approx {eps_star:.2f}$",
            xy=(eps_star, 0.0),
            xytext=(eps_star + 0.05 * (eps_grid[-1] - eps_grid[0]), max(margin_vals) * 0.5),
            fontsize=8,
            color="tomato",
            arrowprops=dict(arrowstyle="->", color="tomato", lw=0.8),
        )

    ax.set_xlabel(r"Claimed $\varepsilon$")
    ax.set_ylabel(
        r"$\widehat{\Delta}_{\mathrm{MMD}}(\varepsilon)"
        r"= \widehat{\mathrm{MMD}}(P,Q) - \tau_{\mathrm{MMD}}(\varepsilon,\delta)$",
        fontsize=8,
    )
    ax.text(
        0.5, -0.18,
        r"Negative values: no positive estimated MMD-surrogate margin.",
        transform=ax.transAxes,
        ha="center", va="top", fontsize=7.5, color="dimgray", style="italic",
    )
    ax.set_ylim(*MARGIN_YLIM)

    fig.tight_layout()

    stem = outdir / f"{mechanism}_mmd_margin"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    return stem


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mechanisms", nargs="+",
                        default=["NonDPGaussian1", "NonDPLaplace1"])
    parser.add_argument("--mechanism-epsilon", type=float, default=0.1)
    parser.add_argument("--delta", type=float, default=1e-5)
    parser.add_argument("--eps-grid", type=str, default="0.01:3.0:80",
                        help="start:stop:num or comma-separated claimed-epsilon values.")
    parser.add_argument("--mmd-samples", type=int, default=5000)
    parser.add_argument("--bandwidth", type=str, default="median20")
    parser.add_argument("--block-size", type=int, default=1024)
    parser.add_argument("--seed", type=int, default=20260522)
    parser.add_argument("--threshold-version", choices=["paper", "notebook"], default="paper")
    parser.add_argument("--outdir", type=str,
                        default="experiments_log/results/mmd_surrogate_diagnostic")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    eps_grid = parse_eps_grid(args.eps_grid)

    summary_rows = []
    sweep_rows = []

    print(f"\n{'Mechanism':<22} {'bw':>10} {'MMD_hat':>10} {'eps_star':>10}")
    print("-" * 55)

    for mech_i, mechanism in enumerate(args.mechanisms):
        seed = args.seed + 100000 * mech_i
        X, Y = sample_mechanism_outputs(
            mechanism=mechanism,
            n_samples=args.mmd_samples,
            mechanism_epsilon=args.mechanism_epsilon,
            seed=seed,
        )

        if args.bandwidth.lower().startswith("median"):
            bw = median_bandwidth(X, Y, n0=20, fallback=1.0)
            bw_label = "median heuristic (n0=20)"
            X_mmd = X[20:] if len(X) > 20 else X
            Y_mmd = Y[20:] if len(Y) > 20 else Y
        else:
            bw = float(args.bandwidth)
            bw_label = f"fixed bw={bw:g}"
            X_mmd, Y_mmd = X, Y

        mmd_hat = estimate_mmd(X_mmd, Y_mmd, bw=bw, block_size=args.block_size)
        eps_star = eps_star_from_mmd(mmd_hat, args.delta)
        print(f"{mechanism:<22} {bw:>10.6g} {mmd_hat:>10.6g} {eps_star:>10.4g}")

        tau_vals = [tau_mmd(float(e), args.delta, args.threshold_version) for e in eps_grid]
        margin_vals = [mmd_hat - t for t in tau_vals]

        for e, t, m in zip(eps_grid, tau_vals, margin_vals):
            sweep_rows.append({
                "mechanism": mechanism,
                "eps_claim": float(e),
                "tau_mmd": t,
                "mmd_hat": mmd_hat,
                "delta_mmd": m,
            })

        summary_rows.append({
            "mechanism": mechanism,
            "mechanism_epsilon": args.mechanism_epsilon,
            "delta": args.delta,
            "bandwidth": bw,
            "bandwidth_rule": bw_label,
            "mmd_hat": mmd_hat,
            "mmd_samples_used": len(X_mmd),
            "eps_star": eps_star,
        })

        stem = plot_mechanism(mechanism, mmd_hat, eps_star, eps_grid,
                              tau_vals, margin_vals, outdir)
        print(f"  -> {stem}.{{pdf,png}}")

    print()

    summary_df = pd.DataFrame(summary_rows)
    sweep_df = pd.DataFrame(sweep_rows)
    summary_df.to_csv(outdir / "mmd_eps_star.csv", index=False)
    sweep_df.to_csv(outdir / "mmd_margin_sweep.csv", index=False)
    print(f"Saved CSVs to {outdir}")

    latex_cols = ["mechanism", "bandwidth", "mmd_hat", "eps_star"]
    (outdir / "mmd_eps_star.tex").write_text(
        summary_df[latex_cols].to_latex(index=False, float_format=lambda x: f"{x:.4g}"),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
