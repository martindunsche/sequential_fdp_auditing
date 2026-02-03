import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

CSV = "mu_results.csv"
OUT = "mu_rejection_blue_1p3_to_1p1.png"

TICK_FONT_SIZE = 22
LABEL_FONT_SIZE = 22

def _smooth_monotone_decreasing(x, y, num=1000, window=31):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.size < 2:
        return x, y

    xg = np.linspace(x.min(), x.max(), num=num)
    yg = np.interp(xg, x, y)

    window = max(5, int(window) | 1)  # odd
    ys = np.convolve(yg, np.ones(window) / window, mode="same")
    ys = np.clip(ys, 0.0, 1.0)

    # enforce y nonincreasing as x increases:
    # reverse -> should be nondecreasing -> use maximum accumulate -> reverse back
    ys = np.maximum.accumulate(ys[::-1])[::-1]
    return xg, ys

df = pd.read_csv(CSV)
df["mu"] = pd.to_numeric(df["mu"], errors="coerce")
df["violation_rate"] = pd.to_numeric(df["violation_rate"], errors="coerce")
df = df.dropna(subset=["mu", "violation_rate"])

# ONLY plot from 1.3 down to 1.1
df = df[(df["mu"] >= 1.1) & (df["mu"] <= 1.3)].sort_values("mu")

x = df["mu"].to_numpy(float)               # increasing for smoothing
y = df["violation_rate"].to_numpy(float)

xs, ys = _smooth_monotone_decreasing(x, y)

plt.rcParams.update({
    "figure.figsize": (10.5, 7.0),
    "font.size": TICK_FONT_SIZE,
    "axes.labelsize": LABEL_FONT_SIZE,
    "xtick.labelsize": TICK_FONT_SIZE,
    "ytick.labelsize": TICK_FONT_SIZE,
    "axes.linewidth": 1.4,
})

fig, ax = plt.subplots()
ax.plot(xs, ys, linewidth=6.0, solid_capstyle="round")  # only the blue curve

ax.set_xlabel(r"$\mu$ (reverse order)")
ax.set_ylabel("Empirical rejection rate")
ax.set_ylim(0, 1.0)

ax.grid(True, axis="y", linestyle="--", linewidth=1.2, alpha=0.25)
ax.grid(False, axis="x")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ticks and reverse axis (left->right: 1.3 ... 1.1)
xt = np.sort(x)[::-1]
ax.set_xticks(xt)
ax.set_xticklabels([f"{v:g}" for v in xt])
ax.invert_xaxis()

fig.tight_layout()
fig.savefig(OUT, dpi=300)
plt.close(fig)
print(f"[✓] Saved {OUT}")
