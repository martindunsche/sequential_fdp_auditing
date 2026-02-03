#!/usr/bin/env Rscript

# real_world.R
# ------------------------------------------------------------
# Runner script that uses PRE-GENERATED mechanism outputs from TWO CSVs
# (P stream + Q stream), without changing any src scripts.
#
# Usage:
#   ./real_world.R <mu> --dataP=IN.csv --dataQ=OUT.csv
#                 [--cores=8] [--trials=250] [--burn=100]
#                 [--alpha=0.05] [--max_iter=10000] [--eval_step=10]
#                 [--refit_threshold=0.10] [--boundary_tweaks=0]
#                 [--qalpha_sims=10000] [--qalpha_kmax=10000]
#
# Each of dataP/dataQ CSV supports either:
#   (1) long: columns trial,value
#   (2) wide: columns trial,P   (for --dataP)   OR columns trial,Q (for --dataQ)
#
# IMPORTANT: This runner assumes each mechanism call returns a SCALAR numeric.
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
source("../src_log/KDE_estimator.R")         # alpha_value / beta_value (may be unused here)
source("../src_log/classifiers.R")           # classifier factories (incl parametric)
source("../src_log/functions_functional.R")  # audit engine
source("../src_log/mechanisms.R")            # mechanisms (may be unused here)

# ---- setup ----
x1 <- rep(0, 10)
x2 <- c(1, rep(0, 9))

args <- commandArgs(trailingOnly = TRUE)

# Required positional arg: mu
positional <- args[!grepl("^--", args)]
if (length(positional) < 1) {
  stop("Usage: ./real_world.R <mu> --dataP=IN.csv --dataQ=OUT.csv [--cores=8] [--trials=250] [--burn=100] ...")
}

mu_arg <- positional[1]
mu_val <- as.numeric(mu_arg)
if (is.na(mu_val)) stop("<mu> must be numeric")

# Required flags for separate CSVs
dataP_arg <- get_flag_value(args, "--dataP", NULL)
dataQ_arg <- get_flag_value(args, "--dataQ", NULL)
if (is.null(dataP_arg) || is.null(dataQ_arg)) {
  stop("Missing required flags: --dataP=PATH.csv and --dataQ=PATH.csv")
}

# Optional flags
cores_arg  <- get_flag_value(args, "--cores",  "8")
trials_arg <- get_flag_value(args, "--trials", "250")
burn_arg   <- get_flag_value(args, "--burn",   "100")

alpha_arg           <- get_flag_value(args, "--alpha", "0.05")
max_iter_arg        <- get_flag_value(args, "--max_iter", "10000")
eval_step_arg       <- get_flag_value(args, "--eval_step", "10")
refit_thr_arg       <- get_flag_value(args, "--refit_threshold", "0.10")
boundary_tweaks_arg <- get_flag_value(args, "--boundary_tweaks", "0")

qalpha_sims_arg <- get_flag_value(args, "--qalpha_sims", "10000")
qalpha_kmax_arg <- get_flag_value(args, "--qalpha_kmax", "10000")

mc_cores <- as.integer(cores_arg)
trials   <- as.integer(trials_arg)
M_burn   <- as.integer(burn_arg)

alpha_level     <- as.numeric(alpha_arg)
max_iter        <- as.integer(max_iter_arg)
eval_step       <- as.integer(eval_step_arg)
refit_threshold <- as.numeric(refit_thr_arg)
boundary_tweaks <- as.logical(as.integer(boundary_tweaks_arg))

qalpha_sims <- as.integer(qalpha_sims_arg)
qalpha_kmax <- as.integer(qalpha_kmax_arg)

if (is.na(mc_cores) || mc_cores < 1) stop("--cores must be a positive integer")
if (is.na(trials)   || trials   < 1) stop("--trials must be a positive integer")
if (is.na(M_burn)   || M_burn   < 1) stop("--burn must be a positive integer")
if (is.na(alpha_level) || alpha_level <= 0 || alpha_level >= 1) stop("--alpha must be in (0,1)")
if (is.na(max_iter) || max_iter < 1) stop("--max_iter must be a positive integer")
if (is.na(eval_step) || eval_step < 1) stop("--eval_step must be a positive integer")
if (is.na(refit_threshold) || refit_threshold < 0) stop("--refit_threshold must be >= 0")
if (is.na(qalpha_sims) || qalpha_sims < 100) stop("--qalpha_sims seems too small")
if (is.na(qalpha_kmax) || qalpha_kmax < 10) stop("--qalpha_kmax seems too small")

# Cap to available cores
avail <- parallel::detectCores(logical = TRUE)
mc_cores <- min(mc_cores, avail)

cat(
  "Running with trials =", trials,
  "burn =", M_burn,
  "mc_cores =", mc_cores,
  "(avail:", avail, ")\n",
  "Reading CSV P:", dataP_arg, "\n",
  "Reading CSV Q:", dataQ_arg, "\n",
  "alpha =", alpha_level,
  "max_iter =", max_iter,
  "eval_step =", eval_step,
  "refit_threshold =", refit_threshold,
  "boundary_tweaks =", boundary_tweaks, "\n"
)

# ---- claimed curve ----
claimed_curve <- function(alpha) {
  stats::pnorm(stats::qnorm(1 - alpha) - mu_val)
}

# ---- choose classifier (PARAMETRIC GAUSSIAN) ----
classifier <- make_kde_classifier()

# ---- load precomputed samples from two CSVs (one-column) and resample per trial ----

datP <- utils::read.csv(dataP_arg, header = TRUE, stringsAsFactors = FALSE)
datQ <- utils::read.csv(dataQ_arg, header = TRUE, stringsAsFactors = FALSE)

extract_numeric_column <- function(dat, which_stream = c("P","Q")) {
  which_stream <- match.arg(which_stream)
  if (ncol(dat) < 1) stop(which_stream, " CSV has no columns")

  # If exactly one column, use it; otherwise prefer 'value' if present
  v <- if (ncol(dat) == 1) dat[[1]] else if ("value" %in% names(dat)) dat[["value"]] else dat[[1]]

  v <- as.numeric(v)
  if (anyNA(v)) stop(which_stream, " CSV: could not parse numeric values (NA introduced).")
  v
}

P_all <- extract_numeric_column(datP, "P")
Q_all <- extract_numeric_column(datQ, "Q")

need_per_trial <- M_burn + max_iter

if (length(P_all) < 1) stop("P CSV has no usable rows")
if (length(Q_all) < 1) stop("Q CSV has no usable rows")

# For each trial, draw a fresh stream of length need_per_trial from the big pool.
# Use replace=TRUE to allow many trials even if the pool isn't huge.
P_list <- replicate(trials, sample(P_all, size = need_per_trial, replace = TRUE), simplify = FALSE)
Q_list <- replicate(trials, sample(Q_all, size = need_per_trial, replace = TRUE), simplify = FALSE)

cat("Loaded P/Q pools and resampled per trial.\n")
cat("pool sizes: |P_all| =", length(P_all), " |Q_all| =", length(Q_all), "\n")
cat("samples per trial =", need_per_trial, " (= burn + max_iter)\n")
cat("Example lengths (trial 1): P =", length(P_list[[1]]), "Q =", length(Q_list[[1]]), "\n")


# ---- per-trial streaming mechanism (no src changes) ----
make_precomputed_mechanism <- function(P_stream, Q_stream, x1, x2) {
  iP <- 0L
  iQ <- 0L
  function(x) {
    is_P <- isTRUE(all.equal(x, x1))
    if (is_P) {
      iP <<- iP + 1L
      if (iP > length(P_stream)) stop("Ran out of precomputed P samples")
      P_stream[iP]
    } else {
      iQ <<- iQ + 1L
      if (iQ > length(Q_stream)) stop("Ran out of precomputed Q samples")
      Q_stream[iQ]
    }
  }
}

# ---- compute q_alpha (same as run_experiment does) ----
q_alpha <- simulate_gaussian_sup_quantile(
  M = M_burn,
  alpha = alpha_level / 2,
  sims = qalpha_sims,
  k_max = qalpha_kmax
)

# ---- run trials (without run_experiment) ----
results_list <- parallel::mclapply(
  seq_len(trials),
  function(i) {
    Mechanism_i <- make_precomputed_mechanism(P_list[[i]], Q_list[[i]], x1, x2)

    tryCatch(
      sequential_audit_simple(
        Mechanism = Mechanism_i,
        x1 = x1,
        x2 = x2,
        M_burn = M_burn,
        h = 0.1,
        eta_search_max = 15,
        claimed_curve = claimed_curve,
        q_alpha = q_alpha,
        classifier = classifier,
        max_iter = max_iter,
        eval_step = eval_step,
        refit_threshold = refit_threshold,
        boundary_tweaks = boundary_tweaks
      ),
      error = function(e) {
        message(sprintf("[Trial %d] Error: %s", i, conditionMessage(e)))
        make_empty_result(conditionMessage(e))
      }
    )
  },
  mc.cores = mc_cores
)

get1 <- function(z, nm) z[[nm]] %||% NA_real_

violations <- vapply(results_list, get1, numeric(1), nm = "violation")
stops      <- vapply(results_list, get1, numeric(1), nm = "stopped_at_n")
T1_last    <- vapply(results_list, get1, numeric(1), nm = "T1_")
T2_last    <- vapply(results_list, get1, numeric(1), nm = "T2_")
s1_last    <- vapply(results_list, get1, numeric(1), nm = "s1")
s2_last    <- vapply(results_list, get1, numeric(1), nm = "s2")
beta_last  <- vapply(results_list, get1, numeric(1), nm = "beta_allowed")

res <- list(
  violation_rate = mean(violations == 1, na.rm = TRUE),
  avg_stop_n     = mean(stops, na.rm = TRUE),
  last_T1        = T1_last,
  last_T2        = T2_last,
  last_s1        = s1_last,
  last_s2        = s2_last,
  beta_allowed   = beta_last,
  stops          = stops,
  raw            = results_list
)

cat("Violation rate:", res$violation_rate, "\n")
cat("Average stopping n:", res$avg_stop_n, "\n")

# ---- diagnostics ----
reasons <- vapply(res$raw, function(z) z$reason %||% "OK", character(1))

cat("\nReason summary:\n")
print(sort(table(reasons), decreasing = TRUE), quote = FALSE)

# ---- save ----
script_name <- get_script_name("real_world")
outdir <- file.path("results", script_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tagP <- tools::file_path_sans_ext(basename(dataP_arg))
tagQ <- tools::file_path_sans_ext(basename(dataQ_arg))
data_tag <- paste0("P", tagP, "_Q", tagQ)
data_tag <- gsub("[^A-Za-z0-9._-]+", "_", data_tag)

outfile <- file.path(outdir, paste0("stops_burn", M_burn, "_mu", mu_arg, "_", data_tag, ".csv"))

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
