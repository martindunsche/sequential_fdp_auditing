#!/usr/bin/env python3
"""MMD-surrogate utilities for the Sec. 4.3 diagnostic.

Only includes the benchmark mechanisms, median-bandwidth selection, MMD
estimation, and MMD threshold/epsilon-star formulas.  No sequential auditor.

Dynamic-sweep convention: at grid value epsilon, the benchmark mechanism is
instantiated with the same epsilon that is used as the claimed privacy level.
Thus P_epsilon and Q_epsilon change across the sweep.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Dict, Sequence, Tuple

import numpy as np
from scipy.spatial.distance import cdist, pdist

Array = np.ndarray


def as_2d(x: Array | Sequence[float] | float) -> Array:
    arr = np.asarray(x, dtype=float)
    if arr.ndim == 0:
        return arr.reshape(1, 1)
    if arr.ndim == 1:
        return arr.reshape(-1, 1)
    return arr.reshape(arr.shape[0], -1)


def rbf_kernel(x: Array, y: Array, bw: float) -> Array:
    bw = float(bw)
    if not np.isfinite(bw) or bw <= 0:
        raise ValueError(f"Bandwidth must be positive and finite, got {bw!r}.")
    X = as_2d(x)
    Y = as_2d(y)
    dist = cdist(X, Y, metric="euclidean")
    return np.exp(-(dist * dist) / (2.0 * bw * bw))


def median_bandwidth(x: Array, y: Array, n0: int = 20, fallback: float = 1.0) -> float:
    """Median heuristic from first n0 samples from each side."""
    X = as_2d(x)
    Y = as_2d(y)
    k = min(int(n0), len(X), len(Y))
    if k <= 1:
        return float(fallback)
    z = np.concatenate([X[:k], Y[:k]], axis=0)
    d = pdist(z, metric="euclidean")
    med = float(np.median(d)) if d.size else float(fallback)
    if not np.isfinite(med) or med <= 0:
        return float(fallback)
    return med


def _block_kernel_sum(x: Array, y: Array, bw: float, block_size: int = 1024) -> float:
    X = as_2d(x)
    Y = as_2d(y)
    bs = int(block_size)
    total = 0.0
    for i in range(0, len(X), bs):
        Xi = X[i:i + bs]
        for j in range(0, len(Y), bs):
            Yj = Y[j:j + bs]
            total += float(rbf_kernel(Xi, Yj, bw=bw).sum())
    return total


def mmd2_biased(x: Array, y: Array, bw: float, block_size: int = 1024) -> float:
    X = as_2d(x)
    Y = as_2d(y)
    n, m = len(X), len(Y)
    if n == 0 or m == 0:
        raise ValueError("Need at least one sample from each distribution.")
    kxx = _block_kernel_sum(X, X, bw=bw, block_size=block_size) / (n * n)
    kyy = _block_kernel_sum(Y, Y, bw=bw, block_size=block_size) / (m * m)
    kxy = _block_kernel_sum(X, Y, bw=bw, block_size=block_size) / (n * m)
    return float(kxx + kyy - 2.0 * kxy)


def estimate_mmd(x: Array, y: Array, bw: float, block_size: int = 1024) -> float:
    return float(math.sqrt(max(mmd2_biased(x, y, bw=bw, block_size=block_size), 0.0)))


def tau_mmd(epsilon: float, delta: float, threshold_version: str = "paper") -> float:
    eps = float(epsilon)
    d = float(delta)
    if threshold_version == "paper":
        return float(math.sqrt(2.0) * (1.0 - 2.0 * (1.0 - d) / (1.0 + math.exp(eps))))
    if threshold_version == "notebook":
        t = math.exp(eps)
        return float(2.0 * (t - 1.0 + 2.0 * d) / (t + 1.0))
    raise ValueError("threshold_version must be 'paper' or 'notebook'.")


def eps_star_from_fixed_mmd(mmd_value: float, delta: float) -> float:
    """eps_star solving tau_MMD(eps_star,delta)=MMD for one fixed P,Q."""
    m = float(mmd_value)
    d = float(delta)
    if m < 0:
        raise ValueError("MMD value must be nonnegative.")
    if m >= math.sqrt(2.0):
        return math.inf
    denom = 1.0 - m / math.sqrt(2.0)
    val = 2.0 * (1.0 - d) / denom - 1.0
    return float(math.log(val)) if val > 0 else float("nan")


@dataclass(frozen=True)
class MechanismSpec:
    name: str
    family: str   # 'laplace' or 'gaussian'
    variant: str  # 'DP', 'NonDP1', or 'NonDP2'


def parse_mechanism_name(name: str) -> MechanismSpec:
    mapping = {
        "DPLaplace": MechanismSpec("DPLaplace", "laplace", "DP"),
        "NonDPLaplace1": MechanismSpec("NonDPLaplace1", "laplace", "NonDP1"),
        "NonDPLaplace2": MechanismSpec("NonDPLaplace2", "laplace", "NonDP2"),
        "DPGaussian": MechanismSpec("DPGaussian", "gaussian", "DP"),
        "NonDPGaussian1": MechanismSpec("NonDPGaussian1", "gaussian", "NonDP1"),
        "NonDPGaussian2": MechanismSpec("NonDPGaussian2", "gaussian", "NonDP2"),
    }
    if name not in mapping:
        raise ValueError(f"Unknown mechanism {name!r}; choices are {sorted(mapping)}")
    return mapping[name]




def default_delta_for_mechanism(name: str, delta_gaussian: float = 1e-5, delta_laplace: float = 0.0) -> float:
    """Default claimed delta used in the comparison: Gaussian uses delta_gaussian, Laplace uses delta_laplace."""
    spec = parse_mechanism_name(name)
    if spec.family == "gaussian":
        return float(delta_gaussian)
    if spec.family == "laplace":
        return float(delta_laplace)
    raise ValueError(f"Unknown noise family {spec.family!r}.")


def _sample_one_dataset(
    data: Sequence[float],
    spec: MechanismSpec,
    eps: float,
    rng: np.random.Generator,
    n_samples: int,
    delta: float = 1e-5,
    min_n: float = 1e-12,
) -> Array:
    """Vectorized samples matching the uploaded R benchmark scripts."""
    eps = float(eps)
    if eps <= 0:
        raise ValueError("epsilon must be positive.")
    x = np.asarray(data, dtype=float)
    n = float(len(x))
    s = float(x.sum())
    size = int(n_samples)

    if spec.family == "laplace":
        tau = rng.laplace(loc=0.0, scale=2.0 / eps, size=size)
        n_tilde = np.maximum(float(min_n), n + tau)
        rho1 = rng.laplace(loc=0.0, scale=2.0 / (n_tilde * eps), size=size)
        rho2 = rng.laplace(loc=0.0, scale=2.0 / (n * eps), size=size)
        if spec.variant == "DP":
            return s / n_tilde + rho1
        if spec.variant == "NonDP1":
            return np.full(size, s / n) + rho2
        if spec.variant == "NonDP2":
            return np.full(size, s / n) + rho1

    if spec.family == "gaussian":
        if spec.variant == "DP":
            # Corrected Gaussian DP version in gaussian_comp_other_paper(1).R.
            eps_part = eps / 2.0
            delta_part = float(delta) / 2.0
            sigma = math.sqrt(2.0 * math.log(1.25 / delta_part)) / eps_part
            tau = rng.normal(loc=0.0, scale=sigma, size=size)
            n_tilde = np.maximum(float(min_n), n + tau)
            z_s = rng.normal(loc=0.0, scale=sigma, size=size)
            return (s + z_s) / n_tilde
        rho2 = rng.normal(loc=0.0, scale=2.0 / (n * eps), size=size)
        if spec.variant == "NonDP1":
            return np.full(size, s / n) + rho2
        if spec.variant == "NonDP2":
            tau = rng.normal(loc=0.0, scale=2.0 / eps, size=size)
            n_tilde = np.maximum(float(min_n), n + tau)
            rho1 = rng.normal(loc=0.0, scale=2.0 / (n_tilde * eps), size=size)
            return np.full(size, s / n) + rho1

    raise ValueError(f"Unsupported mechanism specification: {spec}")


def sample_mechanism_outputs(
    mechanism: str,
    epsilon: float,
    n_samples: int,
    seed: int,
    delta: float = 1e-5,
    data_p: Sequence[float] = (0.0,),
    data_q: Sequence[float] = (0.0, 1.0),
) -> Tuple[Array, Array]:
    spec = parse_mechanism_name(mechanism)
    rng = np.random.default_rng(int(seed))
    x = _sample_one_dataset(data_p, spec, epsilon, rng, n_samples, delta=delta)
    y = _sample_one_dataset(data_q, spec, epsilon, rng, n_samples, delta=delta)
    return x.reshape(-1, 1), y.reshape(-1, 1)


def mmd_dynamic_row(
    mechanism: str,
    epsilon: float,
    delta: float,
    n_samples: int,
    seed: int,
    rep: int = 0,
    bandwidth_n0: int = 20,
    block_size: int = 1024,
    threshold_version: str = "paper",
) -> Dict[str, Any]:
    """One independent dynamic-diagnostic row.

    The mechanism parameter and the claimed epsilon are both equal to `epsilon`.
    The bandwidth is recomputed from this row's own samples.
    """
    X, Y = sample_mechanism_outputs(
        mechanism=mechanism,
        epsilon=epsilon,
        n_samples=n_samples,
        seed=seed,
        delta=delta,
    )
    bw = median_bandwidth(X, Y, n0=bandwidth_n0, fallback=1.0)
    n0 = min(int(bandwidth_n0), len(X), len(Y))
    X_mmd = X[n0:] if len(X) > n0 else X
    Y_mmd = Y[n0:] if len(Y) > n0 else Y
    mmd_hat = estimate_mmd(X_mmd, Y_mmd, bw=bw, block_size=block_size)
    tau = tau_mmd(epsilon, delta=delta, threshold_version=threshold_version)
    eps_star = eps_star_from_fixed_mmd(mmd_hat, delta=delta)
    margin = float(mmd_hat - tau)
    return {
        "mechanism": mechanism,
        "epsilon": float(epsilon),
        "eps_claim": float(epsilon),
        "eps_mechanism": float(epsilon),
        "delta": float(delta),
        "rep": int(rep),
        "seed": int(seed),
        "bandwidth": float(bw),
        "bandwidth_n0": int(n0),
        "mmd_hat": float(mmd_hat),
        "tau_mmd": float(tau),
        "delta_mmd": margin,
        "eps_star_for_fixed_mechanism": float(eps_star),
        "claim_above_eps_star": bool(np.isfinite(eps_star) and float(epsilon) >= float(eps_star)),
        "positive_mmd_margin": bool(margin > 0.0),
        "mmd_samples_drawn": int(n_samples),
        "mmd_samples_used": int(len(X_mmd)),
    }


def crossing_points(eps: Sequence[float], margin: Sequence[float]) -> list[float]:
    x = np.asarray(eps, dtype=float)
    y = np.asarray(margin, dtype=float)
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    out: list[float] = []
    for i in range(len(x) - 1):
        y0, y1 = y[i], y[i + 1]
        if not (np.isfinite(y0) and np.isfinite(y1)):
            continue
        if y0 == 0.0:
            out.append(float(x[i]))
        if y0 * y1 < 0.0:
            t = -y0 / (y1 - y0)
            out.append(float(x[i] + t * (x[i + 1] - x[i])))
    if len(x) and np.isfinite(y[-1]) and y[-1] == 0.0:
        out.append(float(x[-1]))
    dedup: list[float] = []
    for z in out:
        if not dedup or abs(z - dedup[-1]) > 1e-10:
            dedup.append(z)
    return dedup
