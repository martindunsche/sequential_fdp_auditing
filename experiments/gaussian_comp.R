#!/usr/bin/env Rscript

# -------------------------
# Helper: robust script name
# -------------------------
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
# -------------------------
# Load source files
# -------------------------
source("../R/kde_estimator.R")
source("../R/classifiers.R")
source("../R/mechanisms.R")
source("../R/audit_engine.R")

# -------------------------
# Problem setup
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: <script>.R <mu> <DP|NonDP>", call. = FALSE)
}
mu_val <- as.numeric(args[1])
if (is.na(mu_val)) {
  stop("<mu> must be numeric.", call. = FALSE)
}

x1 <- c(0)
x2 <- c(0, 1)

# Non-DP Gaussian: uses true sample size n (NonDPGaussian1)
make_nonDP_gaussian <- function(mu_val ) {
  function(x) {
    n <- length(x)
    # Private mean uses true n
    mean_nonpriv <- sum(x) / n
    rho <- stats::rnorm(1, mean = 0, sd = 2/(n * mu_val))
    return(mean_nonpriv + rho)
  }
}

# DP Gaussian: uses privatized sample size ñ (DPGaussian1)
make_DP_gaussian <- function(mu_val) {
  function(x) {
    n <- length(x)
    tau <- stats::rnorm(1, mean = 0, sd = 2/mu_val)
    n_tilde <- max(1e-12, n + tau)

    mean_priv <- sum(x) / n_tilde
    rho <- stats::rnorm(1, mean = 0, sd = 2/(n_tilde * mu_val))

    return(mean_priv + rho)
  }
}


# -------------------------
# Parse command-line args
# -------------------------

mech_opt <- toupper(args[2])
if (!mech_opt %in% c("DP", "NONDP")) {
  stop("<mechanism> must be one of: DP, NonDP", call. = FALSE)
}

Mechanism <- switch(
  mech_opt,
  "DP"    = make_DP_gaussian(mu_val ),
  "NONDP" = make_nonDP_gaussian(mu_val)
)

cat("Running mechanism:", mech_opt, "\n")

# Optional flags
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

# Cap to available cores
avail <- parallel::detectCores(logical = TRUE)
mc_cores <- min(mc_cores, avail)

cat(
  "Running with trials =", trials,
  "burn =", M_burn,
  "mc_cores =", mc_cores,
  "(avail:", avail, ")\n"
)

# -------------------------
# Claimed privacy curve
# -------------------------
claimed_curve <- function(alpha) {
  stats::pnorm(stats::qnorm(1 - alpha) - mu_val)
}
classifier <- make_kde_classifier()

# -------------------------
# Run experiment
# -------------------------
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

# -------------------------
# Save results
# -------------------------
script_name <- get_script_name("gauss_eps")
outdir <- file.path("../results", script_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(outdir, paste0("stops_burn", M_burn, "_", mu_val, mech_opt,".csv"))

pad <- function(x, n, fill = NA) {
  if (is.null(x) || length(x) == 0) return(rep(fill, n))
  if (length(x) == 1) return(rep(x, n))
  if (length(x) == n) return(x)
  stop("Unexpected length: got ", length(x), " expected 0,1,", n)
}

n_trials <- length(res$stops)

out_df <- data.frame(
  trial        = seq_len(n_trials),
  stopped_at_n = res$stops,
  T1_          = pad(res$last_T1, n_trials),
  T2_          = pad(res$last_T2, n_trials),
  s1           = pad(res$last_s1, n_trials),
  s2           = pad(res$last_s2, n_trials),
  beta_allowed = pad(res$beta_allowed, n_trials),
  reason       = pad(res$reasons, n_trials, fill = NA_character_),
  row.names    = NULL
)


write.csv(out_df, outfile, row.names = FALSE)
cat("Saved stopping times to:", outfile, "\n")