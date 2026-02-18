#!/usr/bin/env Rscript

# ------------------------------------------------------------
# helpers
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

# ------------------------------------------------------------
# sources
# ------------------------------------------------------------

source("../R/kde_estimator.R")
source("../R/classifiers.R")
source("../R/mechanisms.R")
source("../R/audit_engine.R")

# ------------------------------------------------------------
# setup
# ------------------------------------------------------------

x1 <- rep(0, 10)
x2 <- c(1, rep(0, 9))

Mechanism <- make_sum_gauss(sigma = 1)

args <- commandArgs(trailingOnly = TRUE)

# positional arg: mu
positional <- args[!grepl("^--", args)]
if (length(positional) < 1) {
  stop("Usage: ./gaussian.R <mu> [--cores=250] [--trials=250] [--burn=100]")
}

mu_arg <- positional[1]
mu_val <- as.numeric(mu_arg)
if (is.na(mu_val)) stop("<mu> must be numeric")

# optional flags
cores_arg  <- get_flag_value(args, "--cores",  "250")
trials_arg <- get_flag_value(args, "--trials", "250")
burn_arg   <- get_flag_value(args, "--burn",   "100")

mc_cores <- as.integer(cores_arg)
trials   <- as.integer(trials_arg)
M_burn   <- as.integer(burn_arg)

if (is.na(mc_cores) || mc_cores < 1) stop("--cores must be positive")
if (is.na(trials)   || trials   < 1) stop("--trials must be positive")
if (is.na(M_burn)   || M_burn   < 1) stop("--burn must be positive")

avail <- parallel::detectCores(logical = TRUE)
mc_cores <- min(mc_cores, avail)

cat(
  "Running with trials =", trials,
  "burn =", M_burn,
  "mc_cores =", mc_cores,
  "(avail:", avail, ")\n"
)

claimed_curve <- function(alpha) {
  stats::pnorm(stats::qnorm(1 - alpha) - mu_val)
}

classifier <- make_kde_classifier()

# ------------------------------------------------------------
# RUN (this must come BEFORE any NA diagnostics)
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# NA DIAGNOSTICS  (THIS IS THE CORRECT PLACE)
# ------------------------------------------------------------

reasons <- vapply(
  res$raw,
  function(z) z$reason %||% "OK",
  character(1)
)

# optional console summary
cat("\nReason summary:\n")
print(sort(table(reasons), decreasing = TRUE))

# ------------------------------------------------------------
# SAVE RESULTS
# ------------------------------------------------------------

script_name <- get_script_name("gaussian")
outdir <- file.path("../results", script_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(outdir, paste0("stops_burn", M_burn, "_", mu_arg, ".csv"))

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
