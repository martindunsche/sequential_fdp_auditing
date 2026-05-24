#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


TICK_FONT_SIZE = 22
LABEL_FONT_SIZE = 22
REJECTION_COLOR = "#0072B2"
STOP_LINE_COLOR = "#E43D3D"
BOX_FACE_COLOR = "#F4A259"
BOX_EDGE_COLOR = "#A66A5F"


def resolve_path(path_value: str, summary_path: Path, script_dir: Path) -> Path:
    raw = Path(path_value)
    candidates = [
        raw,
        Path.cwd() / raw,
        script_dir / raw,
        summary_path.parent / raw,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return raw


def runtime_values(csv_path: Path) -> np.ndarray:
    if not csv_path.exists():
        return np.array([], dtype=float)

    df = pd.read_csv(csv_path)
    if "stopped_at_n" not in df.columns:
        return np.array([], dtype=float)

    stops = pd.to_numeric(df["stopped_at_n"], errors="coerce")
    return stops[stops.notna()].to_numpy(float)


def runtime_column(df: pd.DataFrame) -> str:
    for column in ["mean_runtime", "mean_stop_all", "mean_stop_reject"]:
        if column in df.columns:
            return column
    raise ValueError("Summary CSV has no runtime column: expected mean_runtime or mean_stop_all.")


def boxplot_whiskers(values: np.ndarray) -> np.ndarray:
    if values.size == 0:
        return np.array([], dtype=float)
    q1, q2, q3 = np.percentile(values, [25, 50, 75])
    iqr = q3 - q1
    lo = values[values >= q1 - 1.5 * iqr].min(initial=q1)
    hi = values[values <= q3 + 1.5 * iqr].max(initial=q3)
    return np.array([lo, q1, q2, q3, hi], dtype=float)


def plot_combined(summary_csv: Path, outfile: Path) -> None:
    script_dir = Path(__file__).resolve().parent
    df = pd.read_csv(summary_csv)
    df["eps"] = pd.to_numeric(df["eps"], errors="coerce")
    df["rejection_rate"] = pd.to_numeric(df["rejection_rate"], errors="coerce")
    runtime_col = runtime_column(df)
    df[runtime_col] = pd.to_numeric(df[runtime_col], errors="coerce")
    df = df.dropna(subset=["eps", "rejection_rate", runtime_col]).sort_values("eps")

    if df.empty:
        raise ValueError(f"No finite combined plot data in {summary_csv}")

    box_vals = [
        runtime_values(resolve_path(csv_path, summary_csv, script_dir))
        for csv_path in df["csv"].astype(str)
    ]
    has_box = [x.size > 0 for x in box_vals]

    whiskers = np.concatenate([boxplot_whiskers(x) for x in box_vals if x.size > 0])
    y_all = np.concatenate([df[runtime_col].to_numpy(float), whiskers])
    y_all = y_all[np.isfinite(y_all)]
    ymin, ymax = float(y_all.min()), float(y_all.max())
    pad = 0.06 * (ymax - ymin)
    if not np.isfinite(pad) or pad <= 0:
        pad = 1.0
    ylim = (max(0.0, ymin - pad), ymax + pad)

    eps = df["eps"].to_numpy(float)
    x_range = (float(eps.min()), float(eps.max()))
    x_pad = 0.03 * (x_range[1] - x_range[0])
    if not np.isfinite(x_pad) or x_pad <= 0:
        x_pad = 0.02

    unique_eps = np.unique(eps)
    if unique_eps.size > 1:
        box_width = min(0.018 * (x_range[1] - x_range[0]), 0.6 * np.diff(unique_eps).min())
    else:
        box_width = 0.02
    if not np.isfinite(box_width) or box_width <= 0:
        box_width = 0.02

    plt.rcParams.update({
        "figure.figsize": (10.5, 7.0),
        "font.size": TICK_FONT_SIZE,
        "axes.labelsize": LABEL_FONT_SIZE,
        "xtick.labelsize": TICK_FONT_SIZE,
        "ytick.labelsize": TICK_FONT_SIZE,
        "axes.linewidth": 1.4,
    })

    fig, ax = plt.subplots()
    ax.plot(
        eps,
        df["rejection_rate"].to_numpy(float),
        color=REJECTION_COLOR,
        linewidth=6.0,
        solid_capstyle="round",
        zorder=3,
    )
    ax.set_xlabel(r"$\epsilon$ for mechanism $M_\epsilon$ and claimed epsilon")
    ax.set_ylabel("Empirical rejection rate")
    ax.set_xlim(x_range[0] - x_pad, x_range[1] + x_pad)
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, axis="y", linestyle="--", linewidth=1.2, alpha=0.25)
    ax.grid(False, axis="x")
    ax.spines["top"].set_visible(False)

    ax_stop = ax.twinx()
    if any(has_box):
        ax_stop.boxplot(
            [x for x in box_vals if x.size > 0],
            positions=eps[has_box],
            widths=box_width,
            showfliers=False,
            patch_artist=True,
            boxprops={"facecolor": BOX_FACE_COLOR, "alpha": 0.45, "edgecolor": BOX_EDGE_COLOR, "linewidth": 1.3},
            whiskerprops={"color": BOX_EDGE_COLOR, "linewidth": 1.3},
            capprops={"color": BOX_EDGE_COLOR, "linewidth": 1.3},
            medianprops={"color": BOX_EDGE_COLOR, "linewidth": 2.0},
        )

    ax_stop.plot(
        eps,
        df[runtime_col].to_numpy(float),
        color=STOP_LINE_COLOR,
        linewidth=2.6,
        marker="o",
        markersize=5.0,
        solid_capstyle="round",
        zorder=4,
    )

    ax_stop.set_ylabel("Runtime (stopped at n)")
    ax_stop.set_ylim(*ylim)
    ax_stop.spines["top"].set_visible(False)

    tick_start = np.floor((x_range[0] + 1e-12) * 5) / 5
    tick_end = np.floor((x_range[1] + 1e-12) * 5) / 5
    xticks = np.arange(tick_start, tick_end + 0.001, 0.2)
    if xticks.size == 0:
        xticks = np.linspace(x_range[0], x_range[1], min(5, len(eps)))
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{x:g}" for x in xticks])

    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)
    print(f"[✓] Saved {outfile}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--summary",
        default="results/nondp_gaussian_auditor_sweep/nondp_gaussian_auditor_summary.csv",
        help="Path to nondp_gaussian_auditor_summary.csv",
    )
    parser.add_argument(
        "--out",
        default="results/nondp_gaussian_auditor_sweep/NonDPGaussian1_rejection_and_avg_stop_by_eps.png",
        help="Output PNG path",
    )
    args = parser.parse_args()

    plot_combined(Path(args.summary), Path(args.out))


if __name__ == "__main__":
    main()
