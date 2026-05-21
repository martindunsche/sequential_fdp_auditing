#!/usr/bin/env python3
"""
Plot sequential-audit empirical rejection-rate curves with Wilson binomial
confidence ribbons.

This script is designed for the CSV files produced by experiments such as

    stops_new_burn50_0.5.csv

with columns like

    trial, stopped_at_n, T1_, T2_, s1, s2, beta_allowed, reason

The important convention is that a row is counted as a rejection only when
`reason` is either "violation" or "violation at burn-in". Rows with
`reason == "max_iter reached without violation"` remain non-rejections for
all plotted sample sizes n.

Example usages
--------------
Plot one CSV:

    python plot_power_curve_wilson_ribbon.py stops_new_burn50_0.5.csv

Recursively scan the current directory for stops*.csv:

    python plot_power_curve_wilson_ribbon.py --root .

Use exact ECDF-style steps instead of light display smoothing:

    python plot_power_curve_wilson_ribbon.py stops_new_burn50_0.5.csv --exact-step

Outputs are named, for example:

    stops_new_burn50_0.5_rejection_curve_wilson95_ribbon.png
"""

from __future__ import annotations

import argparse
import fnmatch
import math
import os
import re
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path
from statistics import NormalDist
from typing import Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Reasons that mean the audit actually rejected the privacy claim.
REJECTION_REASONS = {
    "violation",
    "violation at burn-in",
}

# Fallback only, used if a CSV has no reason column.
# Prefer the reason column whenever it exists.
DEFAULT_CENSOR_VALUES = {10050.0, 10100.0}


@dataclass(frozen=True)
class CurveData:
    x: np.ndarray
    phat: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    n_trials: int
    n_rejections: int
    x_start: float
    x_right: float
    reason_counts: Optional[pd.Series]


def wilson_interval(successes: np.ndarray, n_trials: int, alpha: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return phat, lower, upper Wilson intervals for binomial proportions."""
    if n_trials <= 0:
        raise ValueError("n_trials must be positive.")
    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must lie in (0, 1).")

    x = np.asarray(successes, dtype=float)
    phat = x / float(n_trials)

    z = NormalDist().inv_cdf(1.0 - alpha / 2.0)
    z2 = z * z
    denom = 1.0 + z2 / n_trials

    center = (phat + z2 / (2.0 * n_trials)) / denom
    half = (z / denom) * np.sqrt(
        phat * (1.0 - phat) / n_trials + z2 / (4.0 * n_trials * n_trials)
    )

    lower = np.clip(center - half, 0.0, 1.0)
    upper = np.clip(center + half, 0.0, 1.0)
    return phat, lower, upper


def infer_burn_from_filename(path: os.PathLike[str] | str) -> Optional[float]:
    """Infer burn-in value from names such as stops_new_burn50_0.5.csv."""
    match = re.search(r"burn(\d+(?:\.\d+)?)", Path(path).name)
    if match is None:
        return None
    return float(match.group(1))


def moving_average_edge(y: np.ndarray, window: int) -> np.ndarray:
    """Centered moving average with edge padding, so endpoints are not zero-padded."""
    y = np.asarray(y, dtype=float)
    if y.size == 0:
        return y

    window = int(window)
    if window <= 1 or y.size < 3:
        return y.copy()

    window = max(3, window | 1)  # odd and at least 3
    window = min(window, y.size if y.size % 2 == 1 else y.size - 1)
    if window <= 1:
        return y.copy()

    pad = window // 2
    padded = np.pad(y, (pad, pad), mode="edge")
    kernel = np.ones(window, dtype=float) / float(window)
    return np.convolve(padded, kernel, mode="valid")


def smooth_curve_for_display(
    phat: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    window: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Light visual smoothing to resemble the existing Figure-5 style.

    The underlying binomial counts are still computed exactly before this
    display smoothing. Use --exact-step to disable smoothing entirely.
    """
    phat_s = moving_average_edge(phat, window)
    lower_s = moving_average_edge(lower, window)
    upper_s = moving_average_edge(upper, window)

    # The rejection-rate curve and Wilson endpoints should be nondecreasing
    # because the event "rejected by time n" is monotone in n. The moving
    # average can introduce tiny local decreases, so enforce monotonicity.
    phat_s = np.maximum.accumulate(np.clip(phat_s, 0.0, 1.0))
    lower_s = np.maximum.accumulate(np.clip(lower_s, 0.0, 1.0))
    upper_s = np.maximum.accumulate(np.clip(upper_s, 0.0, 1.0))

    # Ensure the ribbon contains the displayed center curve.
    lower_s = np.minimum(lower_s, phat_s)
    upper_s = np.maximum(upper_s, phat_s)
    return phat_s, lower_s, upper_s


def get_rejection_times(
    df: pd.DataFrame,
    csv_path: os.PathLike[str] | str,
    fallback_censor_values: Sequence[float] = tuple(DEFAULT_CENSOR_VALUES),
) -> Tuple[np.ndarray, int, Optional[pd.Series]]:
    """
    Return sorted rejection times, denominator N, and reason-count summary.

    If the `reason` column is present, only explicit violation rows count as
    rejections. If the column is missing, fall back to the old behavior of
    treating stopped_at_n values as rejection times except known censor values.
    """
    if "stopped_at_n" not in df.columns:
        raise ValueError(f"{csv_path} has no required column 'stopped_at_n'.")

    stopped = pd.to_numeric(df["stopped_at_n"], errors="coerce")
    n_trials = int(len(df))
    if n_trials == 0:
        raise ValueError(f"{csv_path} has zero rows.")

    reason_counts: Optional[pd.Series] = None

    if "reason" in df.columns:
        reason_raw = df["reason"].fillna("missing").astype(str).str.strip()
        reason_norm = reason_raw.str.lower()
        reason_counts = reason_raw.value_counts(dropna=False)
        is_rejection = reason_norm.isin(REJECTION_REASONS)
        usable = is_rejection & stopped.notna()

        n_bad_rejection_times = int((is_rejection & stopped.isna()).sum())
        if n_bad_rejection_times:
            print(
                f"[warning] {csv_path}: {n_bad_rejection_times} rejection rows have nonnumeric stopped_at_n; "
                "they cannot be placed on the curve and are treated as non-rejections.",
                file=sys.stderr,
            )

        rejection_times = stopped.loc[usable].to_numpy(dtype=float)
    else:
        print(
            f"[warning] {csv_path}: no 'reason' column found; assuming numeric stopped_at_n values are "
            f"rejections except censor values {sorted(fallback_censor_values)}.",
            file=sys.stderr,
        )
        numeric_stops = stopped.loc[stopped.notna()].to_numpy(dtype=float)
        censor = np.asarray(list(fallback_censor_values), dtype=float)
        if censor.size:
            is_censored = np.isin(numeric_stops, censor)
            rejection_times = numeric_stops[~is_censored]
        else:
            rejection_times = numeric_stops

    rejection_times = np.asarray(rejection_times, dtype=float)
    rejection_times = rejection_times[np.isfinite(rejection_times)]
    rejection_times.sort()
    return rejection_times, n_trials, reason_counts


def choose_x_limits(
    df: pd.DataFrame,
    csv_path: os.PathLike[str] | str,
    rejection_times: np.ndarray,
    x_start_override: Optional[float],
) -> Tuple[float, float]:
    stopped = pd.to_numeric(df["stopped_at_n"], errors="coerce")
    finite_stopped = stopped.loc[stopped.notna() & np.isfinite(stopped)].to_numpy(dtype=float)

    if x_start_override is not None:
        x_start = float(x_start_override)
    else:
        inferred = infer_burn_from_filename(csv_path)
        if inferred is not None:
            x_start = inferred
        elif finite_stopped.size:
            x_start = float(np.nanmin(finite_stopped))
        else:
            x_start = 0.0

    if finite_stopped.size:
        x_right = float(np.nanmax(finite_stopped))
    elif rejection_times.size:
        x_right = float(np.nanmax(rejection_times))
    else:
        x_right = x_start + 1.0

    # If every row stops at/before x_start, leave a small visible x-range.
    if x_right <= x_start:
        x_right = x_start + 1.0

    return x_start, x_right


def build_curve_data(
    csv_path: os.PathLike[str] | str,
    alpha: float,
    exact_step: bool,
    dense_points: int,
    smooth_window: int,
    x_start_override: Optional[float],
    fallback_censor_values: Sequence[float] = tuple(DEFAULT_CENSOR_VALUES),
) -> CurveData:
    df = pd.read_csv(csv_path)
    rejection_times, n_trials, reason_counts = get_rejection_times(
        df, csv_path, fallback_censor_values=fallback_censor_values
    )
    x_start, x_right = choose_x_limits(df, csv_path, rejection_times, x_start_override)

    if exact_step:
        # Evaluate only at the x locations needed for an exact ECDF-style curve.
        grid_parts = [np.asarray([x_start, x_right], dtype=float)]
        if rejection_times.size:
            grid_parts.append(rejection_times[rejection_times >= x_start])
        x = np.unique(np.concatenate(grid_parts))
        x.sort()
    else:
        # Dense grid for a visually smooth ribbon/curve.
        dense_points = max(50, int(dense_points))
        x = np.linspace(x_start, x_right, num=dense_points)

    successes = np.searchsorted(rejection_times, x, side="right")
    phat, lower, upper = wilson_interval(successes, n_trials, alpha=alpha)

    if not exact_step:
        phat, lower, upper = smooth_curve_for_display(phat, lower, upper, window=smooth_window)

    return CurveData(
        x=x,
        phat=phat,
        lower=lower,
        upper=upper,
        n_trials=n_trials,
        n_rejections=int(rejection_times.size),
        x_start=x_start,
        x_right=x_right,
        reason_counts=reason_counts,
    )


def make_output_path(csv_path: os.PathLike[str] | str, outdir: Optional[str], suffix: str) -> Path:
    csv_path = Path(csv_path)
    destination_dir = Path(outdir) if outdir is not None else csv_path.parent
    destination_dir.mkdir(parents=True, exist_ok=True)
    return destination_dir / f"{csv_path.stem}{suffix}.png"


def plot_one(
    csv_path: os.PathLike[str] | str,
    alpha: float,
    gamma: float,
    exact_step: bool,
    dense_points: int,
    smooth_window: int,
    x_start: Optional[float],
    outdir: Optional[str],
    suffix: str,
    dpi: int,
    legend: bool,
    title: Optional[str],
    fallback_censor_values: Sequence[float],
) -> Path:
    curve = build_curve_data(
        csv_path=csv_path,
        alpha=alpha,
        exact_step=exact_step,
        dense_points=dense_points,
        smooth_window=smooth_window,
        x_start_override=x_start,
        fallback_censor_values=fallback_censor_values,
    )

    out_path = make_output_path(csv_path, outdir=outdir, suffix=suffix)
    ci_level = 100.0 * (1.0 - alpha)

    plt.rcParams.update(
        {
            "figure.figsize": (10.5, 7.0),
            "font.size": 16,
            "axes.labelsize": 18,
            "xtick.labelsize": 15,
            "ytick.labelsize": 15,
            "axes.linewidth": 1.4,
        }
    )

    fig, ax = plt.subplots()

    curve_color = "#1f77b4"
    alpha_color = "red"

    if exact_step:
        ax.fill_between(
            curve.x,
            curve.lower,
            curve.upper,
            step="post",
            color=curve_color,
            alpha=0.18,
            linewidth=0,
            label=f"Pointwise {ci_level:.0f}% Wilson CI",
        )
        ax.step(
            curve.x,
            curve.phat,
            where="post",
            color=curve_color,
            linewidth=5.0,
            solid_capstyle="round",
            label="Empirical rejection rate",
        )
    else:
        ax.fill_between(
            curve.x,
            curve.lower,
            curve.upper,
            color=curve_color,
            alpha=0.18,
            linewidth=0,
            label=f"Pointwise {ci_level:.0f}% Wilson CI",
        )
        ax.plot(
            curve.x,
            curve.phat,
            color=curve_color,
            linewidth=5.0,
            solid_capstyle="round",
            label="Empirical rejection rate",
        )

    ax.axhline(
        gamma,
        color=alpha_color,
        linestyle="--",
        linewidth=5.0,
        alpha=0.95,
        label=rf"$\gamma={gamma:g}$",
    )

    if curve.x_right - curve.x_start < 1e-9:
        ax.set_xlim(curve.x_start, curve.x_start + 1.0)
    else:
        ax.set_xlim(curve.x_start, curve.x_right)

    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("")
    ax.set_ylabel("Empirical rejection rate")
    if title:
        ax.set_title(title)

    ax.grid(True, axis="y", linestyle="--", linewidth=1.2, alpha=0.25)
    ax.grid(False, axis="x")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if legend:
        ax.legend(frameon=False, loc="lower right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)

    final_rate = curve.n_rejections / float(curve.n_trials)
    print(
        f"[✓] {csv_path} -> {out_path} | "
        f"N={curve.n_trials}, rejections={curve.n_rejections}, final_rate={final_rate:.4f}"
    )
    if curve.reason_counts is not None:
        compact_counts = ", ".join(f"{idx}: {val}" for idx, val in curve.reason_counts.items())
        print(f"    reasons: {compact_counts}")

    return out_path


def iter_csvs_from_paths(paths: Sequence[str], pattern: str) -> Iterable[Path]:
    seen = set()
    for raw in paths:
        path = Path(raw)
        if path.is_dir():
            for root, dirs, files in os.walk(path):
                dirs[:] = [d for d in dirs if "_comp_" not in d and d != ".ipynb_checkpoints"]
                for name in files:
                    if not fnmatch.fnmatch(name, pattern):
                        continue
                    if name.endswith("-checkpoint.csv"):
                        continue
                    candidate = Path(root) / name
                    resolved = candidate.resolve()
                    if resolved not in seen:
                        seen.add(resolved)
                        yield candidate
        else:
            if path.name.endswith("-checkpoint.csv"):
                continue
            resolved = path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                yield path


def parse_censor_values(values: str) -> List[float]:
    if values.strip() == "":
        return []
    out: List[float] = []
    for item in values.split(","):
        item = item.strip()
        if not item:
            continue
        out.append(float(item))
    return out


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot empirical rejection-rate curves with shaded Wilson binomial confidence ribbons."
    )
    parser.add_argument(
        "paths",
        nargs="*",
        help="CSV files or directories. If omitted, recursively scans --root for stops*.csv.",
    )
    parser.add_argument("--root", default=".", help="Root directory to scan when no positional paths are given.")
    parser.add_argument("--pattern", default="stops*.csv", help="Filename pattern used when scanning directories.")
    parser.add_argument("--alpha", type=float, default=0.05, help="CI alpha; 0.05 gives pointwise 95%% intervals.")
    parser.add_argument("--gamma", type=float, default=0.05, help="Horizontal significance-level reference line.")
    parser.add_argument(
        "--exact-step",
        action="store_true",
        help="Disable display smoothing and plot the exact ECDF/Wilson ribbon as a step function.",
    )
    parser.add_argument(
        "--dense-points",
        type=int,
        default=2500,
        help="Number of x-grid points for the smoothed display curve.",
    )
    parser.add_argument(
        "--smooth-window",
        type=int,
        default=41,
        help="Moving-average window for display smoothing. Ignored with --exact-step.",
    )
    parser.add_argument(
        "--x-start",
        type=float,
        default=None,
        help="Override the left x-axis value. By default this is inferred from burnXX in the filename.",
    )
    parser.add_argument("--outdir", default=None, help="Optional output directory. Default: next to each CSV.")
    parser.add_argument(
        "--suffix",
        default=None,
        help="Output filename suffix before .png. Default: _rejection_curve_wilson95_ribbon for alpha=0.05.",
    )
    parser.add_argument("--dpi", type=int, default=300, help="PNG resolution.")
    parser.add_argument("--legend", action="store_true", help="Show a legend on each plot.")
    parser.add_argument("--title", default=None, help="Optional title for every plot.")
    parser.add_argument(
        "--fallback-censor-values",
        default="10050,10100",
        help="Comma-separated censor values used only when a CSV has no reason column.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    ci_level = int(round(100.0 * (1.0 - args.alpha)))
    suffix = args.suffix or f"_rejection_curve_wilson{ci_level}_ribbon"
    fallback_censor_values = parse_censor_values(args.fallback_censor_values)

    input_paths = args.paths if args.paths else [args.root]
    csvs = list(iter_csvs_from_paths(input_paths, pattern=args.pattern))

    if not csvs:
        print(f"[!] No CSV files found for pattern {args.pattern!r}.", file=sys.stderr)
        return 1

    failures = 0
    for csv_path in csvs:
        try:
            plot_one(
                csv_path=csv_path,
                alpha=args.alpha,
                gamma=args.gamma,
                exact_step=args.exact_step,
                dense_points=args.dense_points,
                smooth_window=args.smooth_window,
                x_start=args.x_start,
                outdir=args.outdir,
                suffix=suffix,
                dpi=args.dpi,
                legend=args.legend,
                title=args.title,
                fallback_censor_values=fallback_censor_values,
            )
        except Exception as exc:  # keep going when scanning many files
            failures += 1
            print("\n[✗] ERROR while processing CSV:", file=sys.stderr)
            print(f"    {csv_path}", file=sys.stderr)
            print(f"    {type(exc).__name__}: {exc}", file=sys.stderr)
            traceback.print_exc()

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
