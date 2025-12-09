import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = "."
ALPHA = 0.05
TOTAL_TRIALS = 200  # adjust if needed
CENSOR_VALUES = {10050, 10100}  # plateaus that are NOT real rejections

def plot_one(csv_path):
    df = pd.read_csv(csv_path)

    stop_times = df["stopped_at_n"].values

    # Determine censor cutoff dynamically
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
        rejection_rate = np.zeros_like(ks, dtype=float)
    else:
        max_display = max(true_stops.max(), 100)
        ks = np.arange(0, int(max_display) + 1)
        counts = np.searchsorted(true_stops, ks, side="right")
        rejection_rate = counts / TOTAL_TRIALS

    base = os.path.splitext(csv_path)[0]
    out_path = base + "_power_curve.png"

    plt.figure(figsize=(7, 5))
    plt.step(ks, rejection_rate, where="post", linewidth=2, label="Empirical rejection rate")
    plt.axhline(ALPHA, color="red", linestyle="--", linewidth=2)
    plt.ylim(0, 1)
    plt.xlabel("Stopped at n")
    plt.ylabel("Empirical rejection rate")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"[✓] Saved {out_path}")


def main():
    for root, dirs, files in os.walk(ROOT):
        for f in files:
            if f.startswith("stops_new") and f.endswith(".csv"):
                plot_one(os.path.join(root, f))

if __name__ == "__main__":
    main()
