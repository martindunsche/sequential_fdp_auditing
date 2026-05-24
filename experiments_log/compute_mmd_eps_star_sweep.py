#!/usr/bin/env python3
"""Sweep epsilon for the MMD-surrogate barrier diagnostic.

Important: the epsilon grid is used BOTH inside the benchmark mechanism and in
its claimed privacy statement.  For each epsilon value, the script:

  1. instantiates the benchmark mechanism at that epsilon,
  2. draws fresh samples from the neighboring pair D=(0), D'=(0,1),
  3. computes an RBF bandwidth from the first bandwidth_n0 samples,
  4. excludes those bandwidth samples from the MMD estimate,
  5. computes MMD_hat, eps_star_MMD, and the margin at the same claimed epsilon.

Thus there is no separate --mechanism-epsilon argument.  A negative margin,
or equivalently epsilon > eps_star_MMD, means the MMD surrogate has no positive
population-level margin for detecting the violation at that claimed epsilon.

The script writes CSV and PNG files only; it does not run the sequential auditor.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import importlib.util
import math
from pathlib import Path
import sys
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Helper loading
# ---------------------------------------------------------------------------


def _load_helper():
    """Load sequential23_mmd.py from this folder or a nearby src_log folder."""
    here = Path(__file__).resolve().parent
    candidates = [here / "sequential23_mmd.py"]
    for parent in [here, *here.parents]:
        candidates.append(parent / "src_log" / "sequential23_mmd.py")
    for path in candidates:
        if path.is_file():
            spec = importlib.util.spec_from_file_location("sequential23_mmd", path)
            if spec is None or spec.loader is None:
                continue
            module = importlib.util.module_from_spec(spec)
            sys.modules["sequential23_mmd"] = module
            spec.loader.exec_module(module)
            return module
    raise RuntimeError("Could not find sequential23_mmd.py next to this script or in src_log/.")


seq23 = _load_helper()


# ---------------------------------------------------------------------------
# Parsing and computation
# ---------------------------------------------------------------------------


def parse_grid(spec: str) -> np.ndarray:
    """Parse 'start:stop:num' or comma-separated values."""
    spec = spec.strip()
    if ":" in spec:
        parts = spec.split(":")
        if len(parts) != 3:
            raise ValueError("Grid must be start:stop:num or comma-separated values.")
        start, stop, num = parts
        return np.linspace(float(start), float(stop), int(num))
    return np.array([float(x) for x in spec.split(",") if x.strip()], dtype=float)


def _task(job: dict[str, Any]) -> dict[str, Any]:
    mechanism = job["mechanism"]
    eps = float(job["epsilon"])
    delta = float(job["delta"])
    seed = int(job["seed"])

    X, Y = seq23.sample_mechanism_outputs(
        mechanism=mechanism,
        epsilon=eps,
        n_samples=int(job["mmd_samples"]),
        seed=seed,
        delta=delta,
    )

    n0 = int(job["bandwidth_n0"])
    bw = seq23.median_bandwidth(X, Y, n0=n0, fallback=1.0)
    X_mmd = X[n0:] if len(X) > n0 else X
    Y_mmd = Y[n0:] if len(Y) > n0 else Y

    mmd_hat = seq23.estimate_mmd(X_mmd, Y_mmd, bw=bw, block_size=int(job["block_size"]))
    eps_star = seq23.eps_star_from_fixed_mmd(mmd_hat, delta)
    tau_at_claim = seq23.tau_mmd(eps, delta)
    margin_at_claim = mmd_hat - tau_at_claim

    return {
        "mechanism": mechanism,
        "epsilon": eps,
        "delta": delta,
        "run": int(job["run"]),
        "seed": seed,
        "bandwidth": bw,
        "mmd_hat": mmd_hat,
        "eps_star": eps_star,
        "tau_mmd_at_claim": tau_at_claim,
        "margin_at_claim": margin_at_claim,
        "claimed_eps_above_eps_star": bool(np.isfinite(eps_star) and eps > eps_star),
        "positive_mmd_margin": bool(margin_at_claim > 0.0),
        "mmd_samples_used": int(len(X_mmd)),
    }


def aggregate(rows: list[dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    group_cols = ["mechanism", "epsilon", "delta"]
    pieces = []
    for keys, sub in df.groupby(group_cols, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_cols, keys))
        row["runs"] = len(sub)
        for col in ["bandwidth", "mmd_hat", "eps_star", "tau_mmd_at_claim", "margin_at_claim"]:
            vals = sub[col].to_numpy(dtype=float)
            row[f"{col}_mean"] = float(np.nanmean(vals))
            row[f"{col}_sd"] = float(np.nanstd(vals, ddof=1)) if len(vals) > 1 else 0.0
            row[f"{col}_q05"] = float(np.nanquantile(vals, 0.05))
            row[f"{col}_q10"] = float(np.nanquantile(vals, 0.10))
            row[f"{col}_q90"] = float(np.nanquantile(vals, 0.90))
            row[f"{col}_q95"] = float(np.nanquantile(vals, 0.95))
        row["positive_margin_rate"] = float(np.mean(sub["positive_mmd_margin"].to_numpy(dtype=bool)))
        row["claimed_eps_above_eps_star_rate"] = float(np.mean(sub["claimed_eps_above_eps_star"].to_numpy(dtype=bool)))
        pieces.append(row)
    return pd.DataFrame(pieces).sort_values(["mechanism", "epsilon"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def _nice_limits(y: np.ndarray) -> tuple[float, float]:
    finite = y[np.isfinite(y)]
    if finite.size == 0:
        return -1.0, 1.0
    ymin = min(float(np.nanmin(finite)), 0.0)
    ymax = max(float(np.nanmax(finite)), 0.0)
    pad = 0.07 * max(ymax - ymin, 1e-9)
    return ymin - pad, ymax + pad


def _zero_crossing(x: np.ndarray, y: np.ndarray) -> float | None:
    """Linear interpolation for the first crossing of y=0, if present."""
    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]
    for i in range(1, len(x)):
        if not (np.isfinite(y[i - 1]) and np.isfinite(y[i])):
            continue
        if y[i - 1] == 0:
            return float(x[i - 1])
        if y[i - 1] * y[i] < 0:
            t = -y[i - 1] / (y[i] - y[i - 1])
            return float(x[i - 1] + t * (x[i] - x[i - 1]))
    return None


def plot_margin_by_eps(agg: pd.DataFrame, mechanism: str, out_png: Path) -> None:
    sub = agg[agg["mechanism"] == mechanism].sort_values("epsilon")
    x = sub["epsilon"].to_numpy(dtype=float)
    y = sub["margin_at_claim_mean"].to_numpy(dtype=float)
    ylo = sub["margin_at_claim_q05"].to_numpy(dtype=float)
    yhi = sub["margin_at_claim_q95"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(5.6, 4.25))
    ax.fill_between(x, y, 0.0, where=(y >= 0.0), interpolate=True, alpha=0.13, color="#4C78A8")
    ax.fill_between(x, y, 0.0, where=(y < 0.0), interpolate=True, alpha=0.14, color="#FF8A80")
    if np.any(np.abs(yhi - ylo) > 1e-12):
        ax.fill_between(x, ylo, yhi, alpha=0.12, color="#3E7FB1", linewidth=0)
    ax.plot(x, y, linewidth=2.4, color="#3E7FB1")
    ax.axhline(0.0, color="black", linestyle="--", linewidth=1.05)

    ymin, ymax = _nice_limits(np.concatenate([y, ylo, yhi]))
    ax.set_xlim(float(x.min()), float(x.max()))
    ax.set_ylim(ymin, ymax)

    cross = _zero_crossing(x, y)
    if cross is not None and x.min() <= cross <= x.max():
        accent = "#FF5A4A"
        ax.axvline(cross, color=accent, linestyle=":", linewidth=1.25)
        label_x = min(cross + 0.10 * (x.max() - x.min()), x.max() - 0.36 * (x.max() - x.min()))
        label_y = ymin + 0.83 * (ymax - ymin)
        ax.annotate(
            rf"estimated crossover $\varepsilon\approx {cross:.2f}$",
            xy=(cross, 0.0),
            xytext=(label_x, label_y),
            arrowprops=dict(arrowstyle="->", color=accent, lw=1.1),
            color=accent,
            fontsize=9.5,
        )

    ax.set_xlabel(r"Mechanism/claimed $\varepsilon$", fontsize=12)
    ax.set_ylabel(
        r"$\widehat{\Delta}_{\rm MMD}(\varepsilon)="
        r"\widehat{\rm MMD}_{\varepsilon}(P,Q)-\tau_{\rm MMD}(\varepsilon,\delta)$",
        fontsize=10,
    )
    ax.tick_params(labelsize=10)
    for spine in ax.spines.values():
        spine.set_linewidth(1.05)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_eps_star_vs_claim(agg: pd.DataFrame, mechanism: str, out_png: Path) -> None:
    sub = agg[agg["mechanism"] == mechanism].sort_values("epsilon")
    x = sub["epsilon"].to_numpy(dtype=float)
    y = sub["eps_star_mean"].to_numpy(dtype=float)
    ylo = sub["eps_star_q05"].to_numpy(dtype=float)
    yhi = sub["eps_star_q95"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(5.6, 4.25))
    xmax = float(x.max())
    ymax = max(float(np.nanmax(yhi[np.isfinite(yhi)])), xmax)
    diag = np.linspace(float(x.min()), xmax, 200)
    ax.plot(diag, diag, color="black", linestyle="--", linewidth=1.05, label=r"$\varepsilon^\star=\varepsilon$")
    if np.any(np.abs(yhi - ylo) > 1e-12):
        ax.fill_between(x, ylo, yhi, alpha=0.12, color="#3E7FB1", linewidth=0)
    ax.plot(x, y, linewidth=2.4, color="#3E7FB1", label=r"mean $\widehat\varepsilon^\star_{\rm MMD}$")
    ax.fill_between(x, y, x, where=(x > y), interpolate=True, alpha=0.14, color="#FF8A80")

    ax.set_xlim(float(x.min()), xmax)
    ax.set_ylim(0.0, 1.07 * max(ymax, 1e-9))
    ax.set_xlabel(r"Mechanism/claimed $\varepsilon$", fontsize=12)
    ax.set_ylabel(r"$\widehat\varepsilon^\star_{\rm MMD}$", fontsize=12)
    ax.tick_params(labelsize=10)
    ax.legend(frameon=False, fontsize=10, loc="best")
    for spine in ax.spines.values():
        spine.set_linewidth(1.05)
    fig.text(
        0.5,
        0.02,
        r"Red region: claimed $\varepsilon$ exceeds $\widehat\varepsilon^\star_{\rm MMD}$, so the MMD surrogate has no positive margin.",
        ha="center",
        va="center",
        fontsize=8.6,
        style="italic",
        color="dimgray",
    )
    fig.tight_layout(rect=(0.0, 0.075, 1.0, 1.0))
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description="Sweep epsilon for MMD epsilon-star diagnostics; CSV + PNG only.")
    parser.add_argument("--mechanisms", nargs="+", default=["NonDPGaussian1", "NonDPLaplace1"])
    parser.add_argument("--eps-grid", default="0.01:3.0:80",
                        help="Epsilon values used both in the mechanism and in the claim: start:stop:num or comma list.")
    parser.add_argument("--delta-gaussian", type=float, default=1e-5)
    parser.add_argument("--delta-laplace", type=float, default=0.0)
    parser.add_argument("--mmd-samples", type=int, default=5000)
    parser.add_argument("--runs", type=int, default=1,
                        help="Independent repetitions per mechanism/epsilon. Parallelization is over runs.")
    parser.add_argument("--cores", type=int, default=1,
                        help="Worker processes used for the repetitions at each mechanism/epsilon.")
    parser.add_argument("--bandwidth-n0", type=int, default=20,
                        help="First n0 samples used for the median bandwidth heuristic and then excluded.")
    parser.add_argument("--block-size", type=int, default=1024)
    parser.add_argument("--seed", type=int, default=20260522)
    parser.add_argument("--outdir", default="mmd_eps_star_sweep")
    args = parser.parse_args()

    eps_grid = parse_grid(args.eps_grid)
    if np.any(eps_grid <= 0):
        raise ValueError("All epsilon values must be positive.")
    if args.runs < 1:
        raise ValueError("--runs must be positive.")
    if args.cores < 1:
        raise ValueError("--cores must be positive.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict[str, Any]] = []
    print("\nMMD epsilon-star sweep")
    print("epsilon grid drives both the mechanism parameter and the claimed epsilon")
    print(f"runs per point = {args.runs}, cores = {args.cores}, samples = {args.mmd_samples}\n")

    for mech_i, mechanism in enumerate(args.mechanisms):
        delta = seq23.default_delta_for_mechanism(
            mechanism,
            delta_gaussian=args.delta_gaussian,
            delta_laplace=args.delta_laplace,
        )
        print(f"Mechanism: {mechanism}  (delta={delta:g})")
        for eps_i, eps in enumerate(eps_grid):
            jobs = []
            for run in range(args.runs):
                seed = args.seed + 10_000_000 * mech_i + 100_000 * eps_i + run
                jobs.append({
                    "mechanism": mechanism,
                    "epsilon": float(eps),
                    "delta": float(delta),
                    "delta_gaussian": float(args.delta_gaussian),
                    "run": run,
                    "seed": seed,
                    "mmd_samples": args.mmd_samples,
                    "bandwidth_n0": args.bandwidth_n0,
                    "block_size": args.block_size,
                })
            if args.cores == 1 or args.runs == 1:
                rows = [_task(job) for job in jobs]
            else:
                with ProcessPoolExecutor(max_workers=min(args.cores, args.runs)) as ex:
                    rows = list(ex.map(_task, jobs))
            all_rows.extend(rows)
            mean_margin = float(np.mean([r["margin_at_claim"] for r in rows]))
            mean_star = float(np.mean([r["eps_star"] for r in rows]))
            print(f"  eps={eps:.4g}  mean margin={mean_margin:.4g}  mean eps_star={mean_star:.4g}")
        print()

    runs_df = pd.DataFrame(all_rows).sort_values(["mechanism", "epsilon", "run"]).reset_index(drop=True)
    agg_df = aggregate(all_rows)

    runs_path = outdir / "mmd_eps_star_runs.csv"
    agg_path = outdir / "mmd_eps_star_by_eps.csv"
    runs_df.to_csv(runs_path, index=False)
    agg_df.to_csv(agg_path, index=False)

    for mechanism in args.mechanisms:
        plot_margin_by_eps(agg_df, mechanism, outdir / f"{mechanism}_mmd_margin_by_eps.png")
        plot_eps_star_vs_claim(agg_df, mechanism, outdir / f"{mechanism}_eps_star_vs_claimed_eps.png")

    print("Saved:")
    print(f"  {runs_path}")
    print(f"  {agg_path}")
    for mechanism in args.mechanisms:
        print(f"  {outdir / (mechanism + '_mmd_margin_by_eps.png')}")
        print(f"  {outdir / (mechanism + '_eps_star_vs_claimed_eps.png')}")


if __name__ == "__main__":
    main()
