#!/usr/bin/env Rscript

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
source("../src_log/classifiers.R")
source("../src_log/functions_functional.R")
source("../src_log/mechanisms.R")
source("../src_log/KDE_estimator.R")

# -------------------------
# Problem setup
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: <script>.R <eps> <DP|NonDP>", call. = FALSE)
}
eps <- as.numeric(args[1])
if (is.na(eps)) {
  stop("<eps> must be numeric.", call. = FALSE)
}
eps_claim <- eps
delta_claim <- 1e-10

x1 <- c(0)
x2 <- c(0, 1)

# Non-DP Gaussian: uses true sample size n (NonDPGaussian1),
# matching the additive-noise benchmark convention in the Laplace script.
make_nonDP_gaussian <- function(eps) {
  function(x) {
    n <- length(x) + 1000000
    mean_nonpriv <- sum(x) / n
    rho <- stats::rnorm(1, mean = 0, sd = 2/(n * eps))  # uses true n
    return(mean_nonpriv + rho)
  }
}

# DP Gaussian: uses privatized sample size n_tilde (DPGaussian1),
# matching the additive-noise benchmark convention in the Laplace script.
make_DP_gaussian <- function(eps) {
  function(x) {
    n <- length(x)

    tau <- stats::rnorm(1, mean = 0, sd = 2/eps)
    n_tilde <- max(1e-12, n + tau)

    mean_priv <- sum(x) / n_tilde

    rho <- stats::rnorm(1, mean = 0, sd = 2/(n_tilde * eps))

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
  "DP"    = make_DP_gaussian(eps),
  "NONDP" = make_nonDP_gaussian(eps)
)

cat("Running mechanism:", mech_opt, "\n")

# Optional flags
cores_arg  <- get_flag_value(args, "--cores",  "250")
trials_arg <- get_flag_value(args, "--trials", "250")
burn_arg   <- get_flag_value(args, "--burn",   "100")

mc_cores <- as.integer(cores_arg)
trials   <- as.integer(trials_arg)
M_burn   <- as.integer(burn_arg)

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
cat("Claimed approximate-DP parameters: eps =", eps_claim,
    "delta =", delta_claim, "\n")

# -------------------------
# Claimed privacy curve
# -------------------------
# f-DP lower bound equivalent to an (eps_claim, delta_claim)-DP claim:
#   f_{eps,delta}(alpha)
#   = max{0, 1 - delta - exp(eps) * alpha,
#            exp(-eps) * (1 - delta - alpha)}.
# This replaces the GDP/Gaussian-shift curve, because the benchmark
# mechanisms here are the approximate-DP Gaussian analogues of the
# Laplace mechanisms rather than equal-variance Gaussian shift pairs.
claimed_curve <- function(alpha) {
  pmax(
    0,
    1 - delta_claim - exp(eps_claim) * alpha,
    exp(-eps_claim) * (1 - delta_claim - alpha)
  )
}

classifier <- make_kde_classifier()

# -------------------------
# Run experiment
# -------------------------
start_time <- Sys.time()
cat("Started at:", format(start_time, "%Y-%m-%d %H:%M:%S %Z"), "\n")

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
  show_progress = TRUE
)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")
cat("Finished at:", format(end_time, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Elapsed time:", round(as.numeric(elapsed), 2), "minutes\n")

# -------------------------
# Save results
# -------------------------
script_name <- get_script_name("gaussian_eps")

outdir <- file.path("results", script_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(outdir, paste0("stops_burn", M_burn, "_", eps, mech_opt, ".csv"))

pad <- function(x, n, fill = NA) {
  if (is.null(x) || length(x) == 0) return(rep(fill, n))
  if (length(x) == 1) return(rep(x, n))
  if (length(x) == n) return(x)
  stop("Unexpected length: got ", length(x), " expected 0,1,", n)
}

n_trials <- length(res$stops)

# ---- robust reason extraction: use res$reasons if meaningful, else fall back to res$raw ----
get_reason1 <- function(z) {
  # try a few plausible locations (covers different run_experiment versions)
  cand <- list(
    z$reason,
    z$stop_reason,
    z$reasons,
    z$meta$reason,
    z$info$reason
  )

  for (r in cand) {
    if (is.null(r) || length(r) == 0) next
    r <- as.character(r[1])
    if (!is.na(r) && nzchar(r)) return(r)
  }
  "OK"
}

reasons <- NULL

if (!is.null(res$reasons) && length(res$reasons) > 0) {
  reasons <- as.character(res$reasons)
  reasons[is.na(reasons) | reasons == ""] <- "OK"
  reasons <- pad(reasons, n_trials, fill = "OK")
} else if (!is.null(res$raw) && length(res$raw) > 0) {
  reasons <- vapply(res$raw, get_reason1, character(1))
} else {
  reasons <- rep("OK", n_trials)
}

# Optional console summary (useful to confirm it's not all NA anymore)
cat("\nReason summary:\n")
print(sort(table(reasons), decreasing = TRUE))

out_df <- data.frame(
  trial        = seq_len(n_trials),
  stopped_at_n = res$stops,
  T1_          = pad(res$last_T1, n_trials),
  T2_          = pad(res$last_T2, n_trials),
  s1           = pad(res$last_s1, n_trials),
  s2           = pad(res$last_s2, n_trials),
  beta_allowed = pad(res$beta_allowed, n_trials),
  reason       = pad(reasons, n_trials, fill = "OK"),
  row.names    = NULL
)

write.csv(out_df, outfile, row.names = FALSE)
cat("Saved stopping times to:", outfile, "\n")
