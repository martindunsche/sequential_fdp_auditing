#!/usr/bin/env Rscript

# dpsgd.R
# ------------------------------------------------------------
# Runner script in the SAME form as gaussian.R, but using the
# noisy SGD mechanism and its claimed_curve.
#
# Usage:
#   ./noisy_sgd.R <tau> [--cores=250] [--trials=250] [--burn=100]
# ------------------------------------------------------------

get_script_name <- function(fallback = "script") {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(fallback)
  script_path <- sub("^--file=", "", file_arg)
  tools::file_path_sans_ext(basename(script_path))
}

get_flag_value <- function(args, flag, default = NULL) {
  key <- paste0(flag, "=")
  hit <- grep(paste0("^", key), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key), "", hit[1])
}

`%||%` <- function(x, y) if (!is.null(x) && length(x) == 1L) x else y

# ---- sources (ORDER MATTERS) ----
source("../src_log/KDE_estimator.R")         # alpha_value / beta_value
source("../src_log/classifiers.R")           # make_kde_classifier()
source("../src_log/functions_functional.R")  # run_experiment()
source("../src_log/mechanisms.R")            # (optional; we define SGD below)

# ---- fixed dataset ----
x2 <- c(1, rep(0, 9))
x1 <- rep(0, 10)

# ---- noisy SGD params (keep as in your script) ----
sigma <- 0.2
theta_0 <- 0
m <- 5
eta_learn <- 0.2
T_ <- 10

make_noisy_sgd <- function(theta_0, eta_learn, sigma, T_, m) {
  function(x) {
    theta <- theta_0
    for (i in seq_len(T_)) {
      xs <- sample(x, size = m)
      theta <- theta - eta_learn * (mean(theta - xs) + stats::rnorm(1, 0, sigma))
    }
    theta
  }
}

Mechanism <- make_noisy_sgd(theta_0, eta_learn, sigma, T_, m)

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)

positional <- args[!grepl("^--", args)]
if (length(positional) < 1) {
  stop("Usage: ./noisy_sgd.R <tau> [--cores=250] [--trials=250] [--burn=100]")
}

tau_arg <- positional[1]
tau <- as.integer(tau_arg)
if (is.na(tau) || tau < 1) stop("<tau> must be a positive integer")

cores_arg  <- get_flag_value(args, "--cores",  "250")
trials_arg <- get_flag_value(args, "--trials", "250")
burn_arg   <- get_flag_value(args, "--burn",   "100")

mc_cores <- as.integer(cores_arg)
trials   <- as.integer(trials_arg)
M_burn   <- as.integer(burn_arg)

if (is.na(mc_cores) || mc_cores < 1) stop("--cores must be a positive integer")
if (is.na(trials)   || trials   < 1) stop("--trials must be a positive integer")
if (is.na(M_burn)   || M_burn   < 1) stop("--burn must be a positive integer")

avail <- parallel::detectCores(logical = TRUE)
mc_cores <- min(mc_cores, avail)

cat(
  "Running with trials =", trials,
  "burn =", M_burn,
  "mc_cores =", mc_cores,
  "(avail:", avail, ")\n"
)

# ---- claimed curve (same math as your script, just structured cleanly) ----

calc_mu <- function(x, m, eta_learn, tau) {
  sum(eta_learn * (1 - eta_learn)^(tau - x) / m)
}

# Build mu_vector
mu_vector <- 0
for (k in seq_len(tau)) {
  combos <- utils::combn(tau, k)
  mu_vals <- apply(combos, 2, calc_mu, m = m, eta_learn = eta_learn, tau = tau)
  mu_vector <- c(mu_vector, mu_vals)
}

sigma_tilde <- eta_learn * sigma *
  sqrt((1 - (1 - eta_learn)^(2 * tau)) / (1 - (1 - eta_learn)^2))

claimed_curve <- function(alpha) {
  mean(stats::pnorm(stats::qnorm(1 - alpha) - mu_vector / sigma_tilde))
}

# ---- choose classifier (KDE) ----
classifier <- make_kde_classifier()

# ---- run ----
res <- run_experiment(
  Mechanism     = Mechanism,
  x1            = x1,
  x2            = x2,
  claimed_curve = claimed_curve,
  classifier    = classifier,
  M_burn        = M_burn,
  h             = 0.1,
  trials        = trials,
  mc_cores      = mc_cores
)

# ---- diagnostics ----
reasons <- vapply(res$raw, function(z) z$reason %||% "OK", character(1))

cat("\nReason summary:\n")
print(sort(table(reasons), decreasing = TRUE), quote = FALSE)

# ---- save ----
script_name <- get_script_name("noisy_sgd")
outdir <- file.path("results", script_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(outdir, paste0("stops_new_burn", M_burn, "_tau", tau_arg, ".csv"))

n_trials <- length(res$stops)

out_df <- data.frame(
  trial        = seq_len(n_trials),
  stopped_at_n = res$stops,
  T1_          = res$last_T1,
  T2_          = res$last_T2,
  s1           = res$last_s1,
  s2           = res$last_s2,
  beta_allowed = res$beta_allowed,
  reason       = reasons,
  row.names    = NULL
)

write.csv(out_df, outfile, row.names = FALSE)
cat("Saved stopping times to:", outfile, "\n")
