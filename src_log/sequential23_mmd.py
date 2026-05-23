#!/usr/bin/env python3
"""MMD surrogate diagnostics and a paper-correct sequential auditor for [23].

This module is intended to be dropped into the repository

    martindunsche/sequential_fdp_auditing

for the comparison discussed in Section 4.3.

It contains two conceptually separate components:

1. A large-sample fixed-batch MMD estimator for estimating the population
   signal m = MMD(P, Q) induced by a chosen kernel/bandwidth.  This is used for
   the MMD-surrogate barrier diagnostic.

2. A lightly adapted implementation of the one-sided sequential two-sample test
   from Gonzalez, Rubio, Ramdas, and Ribero, "Sequentially Auditing
   Differential Privacy" ([23] in the paper draft).  The OnlineNewtonStep,
   OnlineGradientAscent, and OneSidedTwoSampleSequentialTest classes are based
   on the public Google Research notebook

   https://github.com/google-research/google-research/blob/master/dp_sequential_test/numerical_experiments.ipynb

   which is released under the Apache License 2.0.  The only substantive change
   here is that the MMD threshold tau is set to the formula stated in the paper,

       tau(eps, delta) = sqrt(2) * (1 - 2(1 - delta)/(1 + exp(eps))).

   The public notebook uses a leading constant 2 in the corresponding expression;
   unless there is an implicit normalization elsewhere, that differs from the
   paper by a factor sqrt(2).  This module defaults to the paper formula, while
   still allowing threshold_version='notebook' for debugging.

Dependencies: numpy, scipy.  matplotlib/pandas are only needed by the runner.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple
import math
import numpy as np
from scipy.spatial.distance import cdist, pdist

Array = np.ndarray


# ---------------------------------------------------------------------------
# Kernel and MMD utilities
# ---------------------------------------------------------------------------

def _as_2d(x: Array | Sequence[float] | float) -> Array:
    """Convert scalar/vector/list input to shape (n, d)."""
    arr = np.asarray(x, dtype=float)
    if arr.ndim == 0:
        return arr.reshape(1, 1)
    if arr.ndim == 1:
        return arr.reshape(-1, 1) if arr.size != 1 else arr.reshape(1, 1)
    return arr.reshape(arr.shape[0], -1)


def ensure_list(input_data):
    """Compatibility helper mirroring the public [23] notebook."""
    if not isinstance(input_data, list):
        return [input_data]
    return input_data


def rbf_kernel(samples_x, samples_y=None, bw: float = 1.0) -> Array:
    """Squared exponential/RBF kernel matrix.

    This follows the public [23] notebook convention:
        K(x, y) = exp(-||x-y||^2/(2*bw^2)).
    """
    if bw <= 0 or not np.isfinite(bw):
        raise ValueError(f"Bandwidth must be positive and finite, got {bw!r}.")
    samples_y = samples_x if samples_y is None else samples_y
    samples_x = ensure_list(samples_x)
    samples_y = ensure_list(samples_y)
    X = np.vstack(samples_x).astype(float)
    Y = np.vstack(samples_y).astype(float)
    distances = cdist(X, Y, "euclidean")
    return np.exp(-(distances * distances) / (2.0 * bw * bw))


def median_bandwidth(x: Array, y: Array, n0: int = 20, fallback: float = 1.0) -> float:
    """Median heuristic bandwidth using the first n0 samples from each side.

    The [23] paper states that the first 20 samples are used to set the bandwidth
    and then excluded from testing in the additive-noise experiments.  This
    function implements that rule.  If all pairwise distances are zero or the
    median is non-finite, it falls back to `fallback`.
    """
    X = _as_2d(x)
    Y = _as_2d(y)
    k = min(n0, len(X), len(Y))
    if k <= 1:
        return float(fallback)
    Z = np.concatenate([X[:k], Y[:k]], axis=0)
    distances = pdist(Z, metric="euclidean")
    med = float(np.median(distances)) if distances.size else float(fallback)
    if not np.isfinite(med) or med <= 0:
        return float(fallback)
    return med


def _block_kernel_sum(X: Array, Y: Array, bw: float, block_size: int = 1024) -> float:
    """Compute sum_{i,j} K(X_i, Y_j) without materializing huge matrices."""
    X = _as_2d(X)
    Y = _as_2d(Y)
    total = 0.0
    for i in range(0, len(X), block_size):
        Xi = X[i : i + block_size]
        for j in range(0, len(Y), block_size):
            Yj = Y[j : j + block_size]
            total += float(rbf_kernel(Xi, Yj, bw).sum())
    return total


def mmd2_biased(x: Array, y: Array, bw: float, block_size: int = 1024) -> float:
    """Biased fixed-batch estimator of MMD^2.

    This is nonnegative up to numerical precision and is stable for diagnostic
    reporting.  It estimates the population object in the MMD surrogate null,
    not the online witness score used inside the betting process.
    """
    X = _as_2d(x)
    Y = _as_2d(y)
    n = len(X)
    m = len(Y)
    if n == 0 or m == 0:
        raise ValueError("Need at least one sample from each distribution.")
    kxx = _block_kernel_sum(X, X, bw, block_size=block_size) / (n * n)
    kyy = _block_kernel_sum(Y, Y, bw, block_size=block_size) / (m * m)
    kxy = _block_kernel_sum(X, Y, bw, block_size=block_size) / (n * m)
    return float(kxx + kyy - 2.0 * kxy)


def estimate_mmd(x: Array, y: Array, bw: float, block_size: int = 1024) -> float:
    return float(math.sqrt(max(mmd2_biased(x, y, bw, block_size=block_size), 0.0)))


def tau_mmd(epsilon: float, delta: float, threshold_version: str = "paper") -> float:
    """MMD threshold tau(eps, delta).

    threshold_version='paper' uses Theorem 3.1 / Definition 3.2 in [23]:
        sqrt(2) * (1 - 2(1-delta)/(1+exp(eps))).

    threshold_version='notebook' uses the public notebook's expression, which is
    equal to 2 * (exp(eps)-1+2*delta)/(exp(eps)+1).  This is provided only for
    debugging/reproducing the notebook behavior.
    """
    eps = float(epsilon)
    d = float(delta)
    if threshold_version == "paper":
        return float(math.sqrt(2.0) * (1.0 - 2.0 * (1.0 - d) / (1.0 + math.exp(eps))))
    if threshold_version == "notebook":
        t = math.exp(eps)
        return float(2.0 * (t - 1.0 + 2.0 * d) / (t + 1.0))
    raise ValueError("threshold_version must be 'paper' or 'notebook'.")


def eps_star_from_mmd(mmd_value: float, delta: float) -> float:
    """Crossover eps where tau_mmd(eps, delta) equals mmd_value.

    Uses the paper threshold.  For m >= sqrt(2), the crossover is infinite,
    because the bounded-kernel MMD cannot force saturation before the maximum.
    """
    m = float(mmd_value)
    d = float(delta)
    if m < 0:
        raise ValueError("MMD value must be nonnegative.")
    if m >= math.sqrt(2.0):
        return math.inf
    denom = 1.0 - m / math.sqrt(2.0)
    val = 2.0 * (1.0 - d) / denom - 1.0
    if val <= 0:
        return float("nan")
    return float(math.log(val))


# ---------------------------------------------------------------------------
# Sequential MMD auditor from [23], with paper-correct tau by default.
# ---------------------------------------------------------------------------

class OnlineNewtonStep:
    """Online Newton Step betting strategy from the public [23] notebook."""

    def __init__(self, tau: float) -> None:
        self.scaling_constant = 2.0 / (2.0 - math.log(3.0))
        self.sum_grads_squared = 1.0
        self.previous_lambda = 0.0
        self.tau = float(tau)

    def next_bet(self, payoff_history: Sequence[float]) -> float:
        if not payoff_history:
            return 0.0
        previous_payoff = float(payoff_history[-1])
        denom = 1.0 + self.previous_lambda * previous_payoff
        # Numerical guard: denom should remain positive by construction.
        denom = max(denom, 1e-12)
        z = -previous_payoff / denom
        self.sum_grads_squared += z ** 2
        lower_limit = 0.0
        upper_limit = 1.0 / (8.0 + 4.0 * self.tau)
        lambd = max(
            min(
                self.previous_lambda - self.scaling_constant * z / self.sum_grads_squared,
                upper_limit,
            ),
            lower_limit,
        )
        self.previous_lambda = float(lambd)
        return float(lambd)


class OnlineGradientAscent:
    """Online Gradient Ascent for learning the MMD witness function.

    The returned value v_t is the current learned-witness score.  It is the
    quantity used by the betting process, but it is not the population MMD.
    """

    def __init__(self) -> None:
        self.m_t = 0.0
        self.gradient_second_moments: List[float] = []
        self.history_products: List[float] = []
        self.history_auxiliaryterm: List[float] = []

    def step(self, x_hist: List[Array], y_hist: List[Array], bw: float) -> float:
        x = x_hist[-1]
        y = y_hist[-1]
        steps = len(x_hist)

        increment = (
            rbf_kernel(x, x, bw)[0, 0]
            + rbf_kernel(y, y, bw)[0, 0]
            - 2.0 * rbf_kernel(x, y, bw)[0, 0]
        )
        increment = float(max(increment, 0.0))
        self.m_t += increment
        self.m_t = max(float(self.m_t), 1e-12)
        self.gradient_second_moments.append(self.m_t)

        if len(x_hist) == 1:
            v_t = 0.0
        else:
            kernel_matrix = (
                rbf_kernel(x_hist[:-1], x, bw)
                - rbf_kernel(x_hist[:-1], y, bw)
                + rbf_kernel(y_hist[:-1], y, bw)
                - rbf_kernel(y_hist[:-1], x, bw)
            ).flatten()
            denom = 2.0 * np.sqrt(np.asarray(self.gradient_second_moments[:-1], dtype=float))
            v_t = float(np.sum(kernel_matrix * np.asarray(self.history_products, dtype=float) / denom))

        aux_term = v_t / math.sqrt(self.m_t) + increment / (4.0 * self.m_t)
        self.history_auxiliaryterm.append(float(aux_term))

        s_t = sum(
            self.history_auxiliaryterm[i] * self.history_products[i] ** 2
            for i in range(steps - 1)
        )
        norm_sq_over_4 = max(s_t + aux_term, 1e-12)
        gamma_t = min(1.0, 1.0 / (2.0 * math.sqrt(norm_sq_over_4)))
        self.history_products.append(1.0)
        self.history_products = [float(z * gamma_t) for z in self.history_products]
        return float(v_t)


class OneSidedTwoSampleSequentialTest:
    """One-sided sequential MMD two-sample test using betting."""

    def __init__(
        self,
        epsilon: float,
        delta: float,
        bw: float,
        threshold_version: str = "paper",
    ) -> None:
        self.wealth = 1.0
        self.wealth_hist: List[float] = [1.0]
        self.tau = tau_mmd(epsilon, delta, threshold_version=threshold_version)
        self.online_newton_step = OnlineNewtonStep(self.tau)
        self.online_gradient_ascent = OnlineGradientAscent()
        self.bandwidth = float(bw)
        self.lambd = 0.0
        self.payoff_history: List[float] = []
        self.v_history: List[float] = []
        self.x_hist: List[Array] = []
        self.y_hist: List[Array] = []

    def step(self, x, y) -> float:
        self.x_hist.append(np.asarray(x, dtype=float))
        self.y_hist.append(np.asarray(y, dtype=float))
        v_t = self.online_gradient_ascent.step(self.x_hist, self.y_hist, self.bandwidth)
        payoff = float(v_t - self.tau)
        multiplier = 1.0 + self.lambd * payoff
        # Numerically, the construction should keep this nonnegative.  We guard
        # to avoid a floating point underflow/negative wealth due to roundoff.
        multiplier = max(multiplier, 0.0)
        self.wealth *= multiplier
        self.wealth_hist.append(float(self.wealth))
        self.v_history.append(float(v_t))
        self.payoff_history.append(payoff)
        self.lambd = self.online_newton_step.next_bet(self.payoff_history)
        return float(self.wealth)


def run_sequential23_test(
    x: Array,
    y: Array,
    epsilon: float,
    delta: float,
    bw: float,
    alpha: float = 0.05,
    max_steps: Optional[int] = None,
    threshold_version: str = "paper",
) -> Tuple[bool, int, float]:
    """Run one paired stream through the [23] sequential auditor.

    Returns: (rejected, stopping_time, final_wealth).  If not rejected,
    stopping_time is max_steps.
    """
    X = _as_2d(x)
    Y = _as_2d(y)
    n = min(len(X), len(Y)) if max_steps is None else min(len(X), len(Y), int(max_steps))
    tester = OneSidedTwoSampleSequentialTest(
        epsilon=epsilon,
        delta=delta,
        bw=bw,
        threshold_version=threshold_version,
    )
    boundary = 1.0 / float(alpha)
    for t in range(n):
        wealth = tester.step(X[t], Y[t])
        if wealth >= boundary:
            return True, t + 1, float(wealth)
    return False, n, float(tester.wealth)


# ---------------------------------------------------------------------------
# Paper-style additive-noise mechanisms for mean estimation.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MechanismSpec:
    name: str
    noise_family: str  # 'gaussian' or 'laplace'
    variant: str       # 'DP', 'NonDP1', 'NonDP2'


def laplace(rng: np.random.Generator, scale: float, size=None) -> Array:
    return rng.laplace(loc=0.0, scale=float(scale), size=size)


def gaussian(rng: np.random.Generator, scale: float, size=None) -> Array:
    return rng.normal(loc=0.0, scale=float(scale), size=size)


def _noise(rng: np.random.Generator, family: str, scale: float, size=None) -> Array:
    if family == "laplace":
        return laplace(rng, scale, size=size)
    if family == "gaussian":
        return gaussian(rng, scale, size=size)
    raise ValueError("family must be 'laplace' or 'gaussian'.")


def parse_mechanism_name(name: str) -> MechanismSpec:
    name = name.strip()
    mapping = {
        "DPLaplace": MechanismSpec(name, "laplace", "DP"),
        "NonDPLaplace1": MechanismSpec(name, "laplace", "NonDP1"),
        "NonDPLaplace2": MechanismSpec(name, "laplace", "NonDP2"),
        "DPGaussian": MechanismSpec(name, "gaussian", "DP"),
        "NonDPGaussian1": MechanismSpec(name, "gaussian", "NonDP1"),
        "NonDPGaussian2": MechanismSpec(name, "gaussian", "NonDP2"),
    }
    if name not in mapping:
        raise ValueError(f"Unknown mechanism {name!r}; choices are {sorted(mapping)}")
    return mapping[name]


def evaluate_mean_mechanism(
    data: Sequence[float],
    spec: MechanismSpec,
    mechanism_epsilon: float,
    rng: np.random.Generator,
    min_n: float = 1e-12,
    privatized_count_scale_factor: float = 2.0,
    output_noise_scale_factor: float = 2.0,
) -> float:
    """Evaluate the additive-noise mean mechanism from [23]-style experiments.

    For Laplace, this follows the formulas printed in the [23] paper:
      n_tilde = max(1e-12, n + tau), tau ~ Laplace(0, 2/epsilon)
      rho1 ~ Laplace(0, 2/(n_tilde*epsilon))
      rho2 ~ Laplace(0, 2/(n*epsilon))

    For Gaussian variants, the paper says the mechanisms are analogous but use
    additive Gaussian noise.  We implement that direct analogy by replacing each
    Laplace draw by a Gaussian draw with the same scale parameter interpreted as
    a standard deviation.  This is configurable through the *_scale_factor args.
    """
    eps = float(mechanism_epsilon)
    if eps <= 0:
        raise ValueError("mechanism_epsilon must be positive for additive-noise mechanisms.")
    x = np.asarray(data, dtype=float)
    n = float(len(x))
    if n <= 0:
        raise ValueError("Dataset must be nonempty.")
    tau_count = float(_noise(rng, spec.noise_family, privatized_count_scale_factor / eps))
    n_tilde = max(float(min_n), n + tau_count)

    rho1_scale = output_noise_scale_factor / (n_tilde * eps)
    rho2_scale = output_noise_scale_factor / (n * eps)
    rho1 = float(_noise(rng, spec.noise_family, rho1_scale))
    rho2 = float(_noise(rng, spec.noise_family, rho2_scale))

    if spec.variant == "DP":
        return float(x.sum() / n_tilde + rho1)
    if spec.variant == "NonDP1":
        return float(x.sum() / n + rho2)
    if spec.variant == "NonDP2":
        return float(x.sum() / n + rho1)
    raise ValueError(f"Unknown variant {spec.variant!r}.")


def sample_mechanism_outputs(
    mechanism: str,
    n_samples: int,
    mechanism_epsilon: float,
    seed: int,
    data_p: Sequence[float] = (0.0,),
    data_q: Sequence[float] = (0.0, 1.0),
    **kwargs,
) -> Tuple[Array, Array]:
    """Draw paired samples from P=A(data_p) and Q=A(data_q)."""
    spec = parse_mechanism_name(mechanism)
    rng = np.random.default_rng(seed)
    X = np.empty(int(n_samples), dtype=float)
    Y = np.empty(int(n_samples), dtype=float)
    for i in range(int(n_samples)):
        X[i] = evaluate_mean_mechanism(data_p, spec, mechanism_epsilon, rng, **kwargs)
        Y[i] = evaluate_mean_mechanism(data_q, spec, mechanism_epsilon, rng, **kwargs)
    return X.reshape(-1, 1), Y.reshape(-1, 1)
