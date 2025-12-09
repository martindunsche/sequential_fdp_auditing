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

# -------------------------
# Load source files
# -------------------------
source("../src_cheap/functions.R")
source("../src_cheap/mechanisms.R")
source("../src_cheap/KDE_estimator.R")

# -------------------------
# Problem setup
# -------------------------
eps <- 0.01
x1 <- c(0)
x2 <- c(0, 1)

# Non-DP Gaussian: uses true sample size n (NonDPGaussian1)
make_nonDP_gaussian <- function(eps) {
  function(x) {
    n <- length(x)
    # Private mean uses true n
    mean_nonpriv <- sum(x) / n
    rho <- stats::rnorm(1, mean = 0, sd = 2/(n * eps))
    return(mean_nonpriv + rho)
  }
}

# DP Gaussian: uses privatized sample size ñ (DPGaussian1)
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
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: <script>.R <mu> <DP|NonDP>", call. = FALSE)
}

mu_val <- as.numeric(args[1])
if (is.na(mu_val)) {
  stop("<mu> must be numeric.", call. = FALSE)
}

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

# -------------------------
# Claimed privacy curve
# -------------------------
claimed_curve <- function(alpha) {
  stats::pnorm(stats::qnorm(1 - alpha) - mu_val)
}

# -------------------------
# Run experiment
# -------------------------
res <- run_experiment(
  Mechanism = Mechanism,
  x1 = x1,
  x2 = x2,
  h = 0.1
)

# -------------------------
# Save results
# -------------------------
script_name <- get_script_name("gauss_eps")
outdir <- file.path("results", script_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(outdir, sprintf("stops_new_100_mu%s_%s.csv",
                                     args[1], mech_opt))

write.csv(
  data.frame(
    trial = seq_along(res$stops),
    stopped_at_n = res$stops,
      T1_          = res$last_T1,
    T2_          = res$last_T2
  ),
  outfile,
  row.names = FALSE
)

cat("Saved stopping times to:", outfile, "\n")
