#!/usr/bin/env Rscript
# profile_single_trial.R
# ---------------------------------------------------------------
# Profiles one complete sequential audit trial using Rprof.
# Run from the repo root:
#   Rscript scripts/profile_single_trial.R
# ---------------------------------------------------------------

source("R/kde_estimator.R")
source("R/classifiers.R")
source("R/mechanisms.R")
source("R/audit_engine.R")

set.seed(42)
mu_val  <- 0.5
x1      <- rep(0, 10)
x2      <- c(1, rep(0, 9))
M_burn  <- 100
trials  <- 1

Mechanism <- make_sum_laplace(sigma = 1)

claimed_curve <- function(alpha) {
  ifelse(
    alpha < exp(-mu_val) / 2,
    1 - exp(mu_val) * alpha,
    ifelse(
      exp(-mu_val) / 2 <= alpha & alpha <= 1 / 2,
      exp(-mu_val) / (4 * alpha),
      ifelse(alpha > 1 / 2, exp(-mu_val) * (1 - alpha), 0)
    )
  )
}

classifier <- make_kde_classifier()

# --- Pre-compute q_alpha (not part of per-trial cost) ---
q_alpha <- simulate_gaussian_sup_quantile(
  M = M_burn, alpha = 0.025, sims = 500, k_max = 500
)
cat(sprintf("q_alpha = %.4f\n", q_alpha))

# --- Profile one trial ---
prof_file <- tempfile(fileext = ".out")
Rprof(prof_file, interval = 0.01, line.profiling = FALSE)

result <- sequential_audit_simple(
  Mechanism     = Mechanism,
  x1            = x1,
  x2            = x2,
  M_burn        = M_burn,
  h             = 0.1,
  eta_search_max = 15,
  claimed_curve = claimed_curve,
  q_alpha       = q_alpha,
  classifier    = classifier,
  max_iter      = 2000L,
  eval_step     = 10L
)

Rprof(NULL)

cat(sprintf("\nTrial stopped at n = %g  (reason: %s)\n\n",
            result$stopped_at_n, result$reason %||% "?"))

cat("=== Profiling summary ===\n")
summary_df <- summaryRprof(prof_file)$by.self
summary_df <- summary_df[order(-summary_df$self.time), ]
print(head(summary_df, 20))

cat("\n=== By total ===\n")
by_total <- summaryRprof(prof_file)$by.total
by_total <- by_total[order(-by_total$total.time), ]
print(head(by_total, 15))

unlink(prof_file)
