import os
import traceback
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = "."
ALPHA = 0.05
CENSOR_VALUES = {10050, 10100}  # plateaus that are NOT real rejections


def _smooth_monotone(x, y, num=2500, window=41):
    """
    Make a visually smooth, nondecreasing curve by:
      1) interpolating to a dense grid
      2) moving-average smoothing
      3) enforcing nondecreasing via cumulative max
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.size < 2:
        return x, y

    xg = np.linspace(float(x[0]), float(x[-1]), num=num)
    yg = np.interp(xg, x, y)

    window = max(5, int(window) | 1)  # ensure odd and >=5
    kernel = np.ones(window) / window
    ys = np.convolve(yg, kernel, mode="same")

    ys = np.clip(ys, 0.0, 1.0)
    ys = np.maximum.accumulate(ys)

    return xg, ys


def plot_one(csv_path: str) -> None:
    df = pd.read_csv(csv_path)

    stop_times = pd.to_numeric(df["stopped_at_n"], errors="coerce").to_numpy()
    stop_times = stop_times[~np.isnan(stop_times)]

    total_trials = len(df)  # denominator = rows in file (incl. censored/missing)

    # censor handling
    unique_vals = set(stop_times)
    censor_cutoff = None
    for v in sorted(unique_vals, reverse=True):
        if v in CENSOR_VALUES:
            censor_cutoff = v
            break

    if censor_cutoff is not None:
        true_stops = stop_times[stop_times < censor_cutoff]
    else:
        true_stops = stop_times

    true_stops = np.sort(true_stops)

    base = os.path.splitext(csv_path)[0]
    out_path = base + "_power_curve_smooth.png"

    plt.rcParams.update({
        "figure.figsize": (10.5, 7.0),
        "font.size": 16,
        "axes.labelsize": 18,
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
        "axes.linewidth": 1.4,
    })

    fig, ax = plt.subplots()

    if true_stops.size == 0:
        ax.plot([0, 1], [0, 0], linewidth=5.0, solid_capstyle="round")
        right = float(xs[-1])
        left = 100.0
        if right - left < 1e-9:
            ax.set_xlim(left, left + 1.0)
        else:
            ax.set_xlim(left, right)

    else:
        # ECDF at unique stop times
        x_unique, counts = np.unique(true_stops, return_counts=True)
        y = np.cumsum(counts) / float(total_trials)

        # --- add a visual anchor at (100, 0) ---
        x0 = 100.0
        if x_unique[0] > x0:
            x_unique = np.insert(x_unique, 0, x0)
            y = np.insert(y, 0, 0.0)
        else:
            # if there are already stops <= 100, force the curve to start at (100,0)
            # by adding an explicit point (100,0) before the rest
            x_unique = np.insert(x_unique, 0, x0)
            y = np.insert(y, 0, 0.0)

        # smooth for display (still monotone)
        xs, ys = _smooth_monotone(x_unique, y, num=2500, window=41)
        ax.plot(xs, ys, linewidth=6.0, solid_capstyle="round")

        # show from x=100 onward (looks like a "power vs n" plot)
        right = float(xs[-1])
        left = 100.0
        if right - left < 1e-9:
            ax.set_xlim(left, left + 1.0)
        else:
            ax.set_xlim(left, right)



    # Thick alpha line
    ax.axhline(ALPHA, color="red", linestyle="--", linewidth=6.0, alpha=0.95)
    

    ax.set_ylim(0, 1.0)
    ax.set_xlabel("Stopped at n")
    ax.set_ylabel("Empirical rejection rate")

    ax.grid(True, axis="y", linestyle="--", linewidth=1.2, alpha=0.25)
    ax.grid(False, axis="x")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)

    print(f"[✓] Saved {out_path}")


def main():
    for root, dirs, files in os.walk(ROOT):
        # skip any subfolder with "_comp_" in its name + Jupyter checkpoints
        dirs[:] = [d for d in dirs if "_comp_" not in d and d != ".ipynb_checkpoints"]

        for f in files:
            if not (f.startswith("stops") and f.endswith(".csv")):
                continue
            if f.endswith("-checkpoint.csv"):
                continue

            csv_path = os.path.join(root, f)
            try:
                plot_one(csv_path)
            except Exception as e:
                print("\n[✗] ERROR while processing CSV:")
                print(f"    {csv_path}")
                print(f"    {type(e).__name__}: {e}")
                traceback.print_exc()
                continue


if __name__ == "__main__":
    main()
