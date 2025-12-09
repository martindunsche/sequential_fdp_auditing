#!/usr/bin/env Rscript

get_script_name <- function(fallback = "script") {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(fallback)
  script_path <- sub("^--file=", "", file_arg)
  tools::file_path_sans_ext(basename(script_path))
}

library(rmutil)
library(KernSmooth)

source("../src_cheap/functions.R")
source("../src_cheap/mechanisms.R")
source("../src_cheap/KDE_estimator.R")

x2 <- c(1, rep(0, 9))
x1 <- rep(0, 10)

n <- 10
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
      theta <- theta - eta_learn * (mean(theta - xs) + rnorm(1, 0, sigma))
    }
    theta
  }
}

Mechanism <- make_noisy_sgd(theta_0, eta_learn, sigma, T_, m)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: noisy_sgd.R <tau>")
tau_arg <- args[1]
tau <- as.numeric(tau_arg)
if (is.na(tau)) stop("<tau> must be numeric")

calc_mu <- function(x, m, eta_learn, tau) {
  sum(eta_learn * (1 - eta_learn)^(tau - x) / m)
}

mu_vector <- 0
for (k in seq_len(tau)) {
  combos <- combn(tau, k)
  mu_vals <- apply(combos, 2, calc_mu, m = m, eta_learn = eta_learn, tau = tau)
  mu_vector <- c(mu_vector, mu_vals)
}

sigma_tilde <- eta_learn * sigma *
  sqrt((1 - (1 - eta_learn)^(2 * tau)) / (1 - (1 - eta_learn)^2))

claimed_curve <- function(alpha) {
  sum(stats::pnorm(stats::qnorm(1 - alpha) - mu_vector / sigma_tilde) / 2^tau)
}

res <- run_experiment(
  Mechanism = Mechanism,
  x1 = x1,
  x2 = x2,
  h = 0.1
)

script_name <- get_script_name("noisy_sgd")
outdir <- file.path("results", script_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(outdir, paste0("stops_new_100", tau_arg, ".csv"))

write.csv(
  data.frame(trial = seq_along(res$stops), stopped_at_n = res$stops),
  outfile,
  row.names = FALSE
)

cat("Saved stopping times to:", outfile, "\n")
