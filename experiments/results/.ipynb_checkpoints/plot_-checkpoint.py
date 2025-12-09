import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = "."
ALPHA = 0.05
TOTAL_TRIALS = 200
CENSOR_VALUES = {10050, 10100}


def load_curve(csv_path):
    df = pd.read_csv(csv_path)
    stop_times = df["stopped_at_n"].values

    # detect censoring
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

    if true_stops.size == 0:
        ks = np.arange(0, 101)
        ys = np.zeros_like(ks, float)
        return ks, ys

    max_display = max(true_stops.max(), 100)
    ks = np.arange(0, int(max_display) + 1)
    counts = np.searchsorted(true_stops, ks, side="right")
    ys = counts / TOTAL_TRIALS
    return ks, ys


def label_from_filename(filename):
    lbl = filename.replace("stops_new", "").replace(".csv", "").lstrip("_")
    return lbl


def plot_group(folder, files, suffix):
    if not files:
        return

    plt.figure(figsize=(8, 6))

    for f in sorted(files):
        ks, ys = load_curve(os.path.join(folder, f))
        plt.step(ks, ys, where="post", linewidth=2, label=label_from_filename(f))

    plt.axhline(ALPHA, color="red", linestyle="--", linewidth=2, label=f"α = {ALPHA}")
    plt.ylim(0, 1)
    plt.xlabel("Stopped at n")
    plt.ylabel("Fraction of trials stopped (Power)")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.title(f"{folder} – {suffix}")
    plt.legend(loc="upper left", fontsize=10)
    plt.tight_layout()

    out = os.path.join(folder, f"{folder}_{suffix}.png")
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"[✓] Saved {out}")


def main():
    for folder in os.listdir(ROOT):
        p = os.path.join(ROOT, folder)
        if not os.path.isdir(p) or folder == "__pycache__":
            continue

        csvs = [f for f in os.listdir(p)
                if f.startswith("stops_new") and f.endswith(".csv")]

        if not csvs:
            continue

        small = [f for f in csvs if "_100" not in f]
        large = [f for f in csvs if "_100" in f]

        plot_group(p, small, "normalN")
        plot_group(p, large, "largeN")


if __name__ == "__main__":
    main()
