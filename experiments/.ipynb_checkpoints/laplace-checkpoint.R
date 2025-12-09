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
# Load code
# -------------------------
source("../src_cheap/functions.R")
source("../src_cheap/mechanisms.R")
source("../src_cheap/KDE_estimator.R")

# -------------------------
# Inputs
# -------------------------
x2 <- c(1, rep(0, 9))
x1 <- rep(0, 10)

Mechanism <- make_sum_laplace(sigma = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: ./laplace.R <mu>")
mu_arg <- args[1]
mu_val <- as.numeric(mu_arg)
if (is.na(mu_val)) stop("<mu> must be numeric")

# -------------------------
# Claimed curve
# -------------------------
claimed_curve <- function(alpha) {
  ifelse(
    alpha < exp(-mu_val) / 2,
    1 - exp(mu_val) * alpha,
    ifelse(
      exp(-mu_val) / 2 <= alpha & alpha <= 1/2,
      exp(-mu_val) / (4 * alpha),
      ifelse(alpha > 1/2, exp(-mu_val) * (1 - alpha), 0)
    )
  )
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
script_name <- get_script_name("laplace")
outdir <- file.path("results", script_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(outdir, paste0("stops_new_100", mu_arg, ".csv"))

write.csv(
  data.frame(trial = seq_along(res$stops), stopped_at_n = res$stops),
  outfile,
  row.names = FALSE
)

cat("Saved stopping times to:", outfile, "\n")
