#!/usr/bin/env Rscript

get_script_name <- function(fallback = "script") {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(fallback)
  script_path <- sub("^--file=", "", file_arg)
  tools::file_path_sans_ext(basename(script_path))
}

source("../src_cheap/functions.R")
source("../src_cheap/KDE_estimator.R")
library(KernSmooth)

# -------------------------
# Command line argument
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: run_plrt_from_kl.R <mu>")
mu_arg <- args[1]
mu_val <- as.numeric(mu_arg)
if (is.na(mu_val)) stop("<mu> must be numeric")

claimed_curve <- function(alpha) {
  stats::pnorm(stats::qnorm(1 - alpha) - mu_val)
}

# -------------------------
# Load KL samples
# -------------------------
kl_x0_raw <- as.numeric(scan("plrt_KL_x0.csv", what = double(), sep = ",", quiet = TRUE))
kl_x1_raw <- as.numeric(scan("plrt_KL_x1.csv", what = double(), sep = ",", quiet = TRUE))

cat("Loaded", length(kl_x0_raw), "x0 KL samples\n")
cat("Loaded", length(kl_x1_raw), "x1 KL samples\n")

make_sequential_mechanism <- function(P_samples, Q_samples) {
  iP <- 0; iQ <- 0
  function(label) {
    if (identical(label, "x0")) {
      iP <<- iP + 1
      if (iP > length(P_samples)) stop("Ran out of x0 samples")
      return(P_samples[iP])
    }
    if (identical(label, "x1")) {
      iQ <<- iQ + 1
      if (iQ > length(Q_samples)) stop("Ran out of x1 samples")
      return(Q_samples[iQ])
    }
    stop("Unknown label")
  }
}

x1 <- "x0"
x2 <- "x1"

# -------------------------
# Perform TRIALS independent runs
# each with a fresh shuffle
# -------------------------
n_trials <- 200
stops <- numeric(n_trials)

for (t in seq_len(n_trials)) {

  # shuffle KL losses BEFORE this run
  kl_x0 <- sample(kl_x0_raw)
  kl_x1 <- sample(kl_x1_raw)

  # fresh mechanism with fresh counters
  Mechanism <- make_sequential_mechanism(kl_x0, kl_x1)

  # run ONE independent sequential audit
  res_t <- run_experiment(
    Mechanism = Mechanism,
    x1 = x1,
    x2 = x2,
    M_burn = 100,
    h = 0.1,
    trials = 1,     # DO NOT run 200 inside here
    refresh_every = 100
  )

  stops[t] <- res_t$stops
}

cat("Violation rate:", mean(!is.na(stops)), "\n")
cat("Average stopping n:", mean(stops, na.rm = TRUE), "\n")

# -------------------------
# Save results
# -------------------------
script_name <- get_script_name("run_plrt_from_kl")
outdir <- file.path("results", script_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(outdir, paste0("stops_new", mu_arg, ".csv"))

write.csv(
  data.frame(
    trial = seq_len(n_trials),
    stopped_at_n = stops
  ),
  outfile,
  row.names = FALSE
)

cat("Saved stopping times to:", outfile, "\n")
