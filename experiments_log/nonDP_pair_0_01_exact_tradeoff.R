# Exact DP discrepancy and tradeoff curves for the NonDP Gaussian mechanism
# Neighboring pair: x = (0), x' = (0, 1)
# Mechanism: M_eps(x) = mean(x) + N(0, (2 / (length(x) * eps_alg))^2)

f_eps_delta <- function(alpha, eps, delta) {
  pmax(0, 1 - delta - exp(eps) * alpha, exp(-eps) * (1 - delta - alpha))
}

pair_params_nonDP_0_01 <- function(eps_alg) {
  list(
    P = list(mu = 0,   sd = 2 / eps_alg),  # x=(0), n=1
    Q = list(mu = 0.5, sd = 1 / eps_alg)   # x'=(0,1), n=2
  )
}

# Intervals where A*x^2 + B*x + C > 0.
.quad_gt_intervals <- function(A, B, C, tol = 1e-14) {
  if (abs(A) < tol) {
    if (abs(B) < tol) {
      if (C > 0) return(list(c(-Inf, Inf))) else return(list())
    }
    r <- -C / B
    if (B > 0) return(list(c(r, Inf))) else return(list(c(-Inf, r)))
  }

  disc <- B^2 - 4 * A * C
  if (disc <= 0) {
    if (A > 0) return(list(c(-Inf, Inf))) else return(list())
  }

  r1 <- (-B - sqrt(disc)) / (2 * A)
  r2 <- (-B + sqrt(disc)) / (2 * A)
  lo <- min(r1, r2)
  hi <- max(r1, r2)

  if (A > 0) list(c(-Inf, lo), c(hi, Inf)) else list(c(lo, hi))
}

.normal_prob_intervals <- function(mu, sd, intervals) {
  if (length(intervals) == 0) return(0)
  total <- 0
  for (iv in intervals) {
    total <- total + pnorm(iv[2], mean = mu, sd = sd) - pnorm(iv[1], mean = mu, sd = sd)
  }
  total
}

# Probability under N(under_mu, under_sd^2) of {log N(mu_num,sd_num^2) / N(mu_den,sd_den^2) > threshold}.
.prob_logratio_gt <- function(mu_num, sd_num, mu_den, sd_den, threshold, under_mu, under_sd) {
  A <- -1 / (2 * sd_num^2) + 1 / (2 * sd_den^2)
  B <-  mu_num / sd_num^2 - mu_den / sd_den^2
  C <- -mu_num^2 / (2 * sd_num^2) + mu_den^2 / (2 * sd_den^2) + log(sd_den / sd_num) - threshold
  intervals <- .quad_gt_intervals(A, B, C)
  .normal_prob_intervals(under_mu, under_sd, intervals)
}

# Hockey-stick divergence D_eps(A || B) = sup_S [A(S) - exp(eps) B(S)]
# for A=N(mu_a,sd_a^2), B=N(mu_b,sd_b^2).
normal_hockey <- function(mu_a, sd_a, mu_b, sd_b, eps_claim) {
  pa <- .prob_logratio_gt(mu_a, sd_a, mu_b, sd_b, eps_claim, mu_a, sd_a)
  pb <- .prob_logratio_gt(mu_a, sd_a, mu_b, sd_b, eps_claim, mu_b, sd_b)
  max(0, pa - exp(eps_claim) * pb)
}

# Exact Neyman-Pearson tradeoff T(P,Q)(alpha).
# Reject P on the region where log(q/p) is large.
normal_tradeoff_exact <- function(alpha, mu_p, sd_p, mu_q, sd_q) {
  prob_p <- function(t) .prob_logratio_gt(mu_q, sd_q, mu_p, sd_p, t, mu_p, sd_p)
  prob_q <- function(t) .prob_logratio_gt(mu_q, sd_q, mu_p, sd_p, t, mu_q, sd_q)

  out <- numeric(length(alpha))
  for (i in seq_along(alpha)) {
    a0 <- alpha[i]
    if (a0 <= 0) { out[i] <- 1; next }
    if (a0 >= 1) { out[i] <- 0; next }

    lower <- -1
    upper <- 1
    iter <- 0
    while (prob_p(lower) < a0 && iter < 200) {
      lower <- 2 * lower - 1
      iter <- iter + 1
    }
    iter <- 0
    while (prob_p(upper) > a0 && iter < 200) {
      upper <- 2 * upper + 1
      iter <- iter + 1
    }

    root <- uniroot(function(t) prob_p(t) - a0, lower = lower, upper = upper,
                    tol = 1e-10)$root
    out[i] <- 1 - prob_q(root)
  }
  out
}

# ---------------- Example ----------------
# Use eps_alg <- 1 if "NonDPGaussian1" means the mechanism itself is fixed at eps=1.
# Use eps_alg <- eps_claim if each plotted x-axis value also changes the mechanism noise.

eps_alg <- 0.6
eps_claim <- 0.6
delta <- 1e-5

params <- pair_params_nonDP_0_01(eps_alg)
P <- params$P
Q <- params$Q

D_PQ <- normal_hockey(P$mu, P$sd, Q$mu, Q$sd, eps_claim)
D_QP <- normal_hockey(Q$mu, Q$sd, P$mu, P$sd, eps_claim)
cat("D_eps(P || Q) =", D_PQ, "\n")
cat("D_eps(Q || P) =", D_QP, "\n")
cat("max bidirectional discrepancy =", max(D_PQ, D_QP), "\n")
cat("delta =", delta, "\n")

alpha <- seq(1e-5, 1 - 1e-5, length.out = 1000)
f <- f_eps_delta(alpha, eps_claim, delta)
T_PQ <- normal_tradeoff_exact(alpha, P$mu, P$sd, Q$mu, Q$sd)
T_QP <- normal_tradeoff_exact(alpha, Q$mu, Q$sd, P$mu, P$sd)

plot(alpha, f, type = "l", lwd = 2, ylim = c(0, 1),
     xlab = expression(alpha), ylab = "tradeoff value / type-II error",
     main = "x=(0), x'=(0,1): f curve and exact tradeoff curves")
lines(alpha, T_PQ, lwd = 2, lty = 2)
lines(alpha, T_QP, lwd = 2, lty = 3)
legend("topright", lwd = 2, lty = c(1, 2, 3), bty = "n",
       legend = c("f_{eps,delta}", "T(P,Q)", "T(Q,P)"))

# Positive vertical gaps against the fixed f_{eps,delta} curve.
plot(alpha, pmax(f - T_PQ, 0), type = "l", lwd = 2,
     xlab = expression(alpha), ylab = "positive vertical gap",
     main = "Positive violations of the tradeoff lower bound")
lines(alpha, pmax(f - T_QP, 0), lwd = 2, lty = 2)
legend("topright", lwd = 2, lty = c(1, 2), bty = "n",
       legend = c("max(f - T(P,Q), 0)", "max(f - T(Q,P), 0)"))
