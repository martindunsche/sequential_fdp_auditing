#!/usr/bin/env Rscript

source("../../src_log/classifiers.R")
source("../../src_log/KDE_estimator.R")
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

mu_txt <- get_arg("--mu")
if (is.na(mu_txt)) stop('Usage: ./gaussian_non_seq.R --mu "0.5"')
mu <- as.numeric(mu_txt)
if (!is.finite(mu) || mu <= 0) stop("--mu must be a positive finite number")

mu_tag <- function(s) gsub("[^0-9A-Za-z._-]+", "_", sprintf("%g", s))

x1 <- c(1, rep(0, 9))
x2 <- rep(0, 10)
Mechanism <- function(x) sum(x) + rnorm(1, mean = 0, sd = 1)

claimed_curve <- function(alpha) pnorm(qnorm(1 - alpha) - mu)

ALPHA <- 0.05
N_values <- c(100, 300, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000)
n_runs <- 250L
mc_cores <- 250

# save under per-mu folder
base_dir <- file.path("nonseq_results", paste0("mu_", mu_tag(mu)))
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

h <- 0.1
eta_search_max <- 15
classifier <- make_kde_classifier()

`%||%` <- function(x, y) if (!is.null(x) && length(x) == 1L) x else y

one_shot_nonseq <- function(N) {
  P <- replicate(N, Mechanism(x1))
  Q <- replicate(N, Mechanism(x2))

  cfg <- list(
    Mechanism = Mechanism,
    x1 = x1,
    x2 = x2,
    h = h,
    eta_search_max = eta_search_max,
    claimed_curve = claimed_curve,
    M_burn = 0L
  )

  atk <- classifier$init(P = P, Q = Q, cfg = cfg)
  if (is.null(atk)) stop("classifier init returned NULL")

  e1 <- classifier$classify(P, atk)
  e2 <- classifier$classify(Q, atk)
  if (anyNA(e1) || anyNA(e2)) stop("classify() returned NA values")

  T1 <- mean(e1)
  T2 <- 1 - mean(e2)

  s1 <- max(sqrt(T1 * (1 - T1)), 1 / (2 * sqrt(N)))
  s2 <- max(sqrt(T2 * (1 - T2)), 1 / (2 * sqrt(N)))
# --- allowed point from (T1, T2): 45°-line intersection with claimed curve ---
f45 <- function(r) claimed_curve(r) - (T2 + (r - T1))

eps  <- 1e-12
grid <- seq(eps, 1 - eps, length.out = 2001)
vals <- vapply(grid, f45, numeric(1))
idx  <- which(vals[-1] * vals[-length(vals)] <= 0)

if (!length(idx)) {
  # fallback: use midpoint-eta from the tradeoff table
  j_mid <- as.integer(ceiling(nrow(td) / 2))
  alpha_allowed <- as.numeric(td$alpha[j_mid])
  beta_allowed  <- claimed_curve(alpha_allowed)
} else {
  lo <- grid[idx[1]]
  hi <- grid[idx[1] + 1]
  alpha_allowed <- tryCatch(stats::uniroot(f45, c(lo, hi))$root,
                            error = function(e) {
                              j_mid <- as.integer(ceiling(nrow(td) / 2))
                              as.numeric(td$alpha[j_mid])
                            })
  beta_allowed <- claimed_curve(alpha_allowed)
}



  crit <- qnorm(ALPHA / 2)
  violation1 <- as.integer(sqrt(N) * (T1 - alpha_allowed) / s1 <= crit)
  violation2 <- as.integer(sqrt(N) * (T2 - beta_allowed)  / s2 <= crit)
  violation  <- as.integer(violation1 == 1L & violation2 == 1L)
  data.frame(
    violation = violation,
    n = N,
    T1 = T1,
    T2 = T2,
    s1 = s1,
    s2 = s2,
    alpha_allowed = alpha_allowed,
    beta_allowed = beta_allowed,
    eta_star = (atk$eta_star %||% NA_real_),
    used_fallback_eta = (atk$used_fallback_eta %||% NA),
    error = NA_character_,              # <-- ADD THIS
    stringsAsFactors = FALSE
  )

}

message("=== mu = ", mu, " | output dir: ", base_dir, " ===")

for (N in N_values) {
  n_dir <- file.path(base_dir, sprintf("N_%04d", N))
  dir.create(n_dir, showWarnings = FALSE, recursive = TRUE)
  message("Starting N = ", N)

  res_list <- mclapply(
    seq_len(n_runs),
    function(i) {
      df <- tryCatch(one_shot_nonseq(N),
                     error = function(e) data.frame(
                       violation = NA_integer_, n = N,
                       T1 = NA_real_, T2 = NA_real_, s1 = NA_real_, s2 = NA_real_,
                       alpha_allowed = NA_real_, beta_allowed = NA_real_,
                       eta_star = NA_real_, used_fallback_eta = NA,
                       error = conditionMessage(e),
                       stringsAsFactors = FALSE
                     ))
      write.csv(df, file.path(n_dir, sprintf("nonseq_run_%04d.csv", i)), row.names = FALSE)
      df
    },
    mc.cores = mc_cores
  )

  big <- do.call(rbind, res_list)
  write.csv(big, file.path(n_dir, "nonseq_all_runs.csv"), row.names = FALSE)

  message("Finished N = ", N, " | violation_rate=", mean(big$violation == 1, na.rm = TRUE))
}
