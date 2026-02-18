#!/usr/bin/env Rscript

# laplace.R
# ------------------------------------------------------------
# Runner script in the SAME form as gaussian.R:
# - positional <mu>
# - flags: --cores, --trials, --burn
# - functional engine + pluggable classifier
# - prints Reason summary
# - writes CSV with burn in filename + reason column
#
# Usage:
#   ./laplace.R <mu> [--cores=250] [--trials=250] [--burn=100]
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
source("../R/kde_estimator.R")    # alpha_value / beta_value
source("../R/classifiers.R")     # make_kde_classifier()
source("../R/mechanisms.R")      # make_sum_laplace()
source("../R/audit_engine.R")    # run_experiment()

# ---- inputs ----
x2 <- c(1, rep(0, 9))
x1 <- rep(0, 10)

Mechanism <- make_sum_laplace(sigma = 1)

args <- commandArgs(trailingOnly = TRUE)

# positional arg: mu
positional <- args[!grepl("^--", args)]
if (length(positional) < 1) {
  stop("Usage: ./laplace.R <mu> [--cores=250] [--trials=250] [--burn=100]")
}

mu_arg <- positional[1]
mu_val <- as.numeric(mu_arg)
if (is.na(mu_val)) stop("<mu> must be numeric")

# flags
cores_arg       <- get_flag_value(args, "--cores",       "250")
trials_arg      <- get_flag_value(args, "--trials",      "250")
burn_arg        <- get_flag_value(args, "--burn",        "100")
qalpha_sims_arg <- get_flag_value(args, "--qalpha_sims", "10000")
qalpha_kmax_arg <- get_flag_value(args, "--qalpha_kmax", "10000")

mc_cores    <- as.integer(cores_arg)
trials      <- as.integer(trials_arg)
M_burn      <- as.integer(burn_arg)
qalpha_sims <- as.integer(qalpha_sims_arg)
qalpha_kmax <- as.integer(qalpha_kmax_arg)

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

# ---- claimed curve (unchanged logic) ----
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
  mc_cores      = mc_cores,
  qalpha_sims   = qalpha_sims,
  qalpha_kmax   = qalpha_kmax
)

# ---- diagnostics ----
reasons <- vapply(res$raw, function(z) z$reason %||% "OK", character(1))

cat("\nReason summary:\n")
print(sort(table(reasons), decreasing = TRUE), quote = FALSE)

# ---- save ----
script_name <- get_script_name("laplace")
outdir <- file.path("../results", script_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(outdir, paste0("stops_new_burn", M_burn, "_", mu_arg, ".csv"))

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
