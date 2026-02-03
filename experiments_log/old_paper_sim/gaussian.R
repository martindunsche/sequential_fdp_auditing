#!/usr/bin/env Rscript
# --------------------------------------------------
# Batch KDE simulation by N (1000 CSVs per N) - base R version
# --------------------------------------------------

library(KernSmooth)
library(parallel)

# --------------------------------------------------
# Functions
# --------------------------------------------------

inner_int <- function(x, eta, grid_width, p, p_q_ratio) {
  p_vector <- p
  p_vector[p_q_ratio <= eta + x] <- 0
  sum(p_vector * grid_width)
}

alpha_value <- function(eta, x_vector, h, grid_width, hat_p, hat_p_q_ratio) {
  1 - sum(
    sapply(
      x_vector,
      inner_int,
      eta = eta,
      grid_width = grid_width,
      p = hat_p,
      p_q_ratio = hat_p_q_ratio
    ) / h * (x_vector[2] - x_vector[1])
  )
}

beta_value <- function(eta, x_vector, h, grid_width, hat_q, hat_p_q_ratio) {
  sum(
    sapply(
      x_vector,
      inner_int,
      eta = eta,
      grid_width = grid_width,
      p = hat_q,
      p_q_ratio = hat_p_q_ratio
    ) / h * (x_vector[2] - x_vector[1])
  )
}

KDE_Estimator <- function(eta_max, Mechanism, x1, x2, N, h) {
  p_sample <- sapply(rep(list(x1), N), Mechanism)
  q_sample <- sapply(rep(list(x2), N), Mechanism)

  t_min <- min(c(p_sample, q_sample))
  t_max <- max(c(p_sample, q_sample))

  bw_p <- dpik(
    p_sample,
    kernel = "normal",
    gridsize = 1000L,
    range.x = c(t_min, t_max)
  )

  bw_q <- dpik(
    q_sample,
    kernel = "normal",
    gridsize = 1000L,
    range.x = c(t_min, t_max)
  )

  KDE_p <- bkde(
    p_sample,
    kernel = "normal",
    bandwidth = bw_p,
    gridsize = 1000L,
    range.x = c(t_min, t_max)
  )

  KDE_q <- bkde(
    q_sample,
    kernel = "normal",
    bandwidth = bw_q,
    gridsize = 1000L,
    range.x = c(t_min, t_max)
  )

  hat_p <- KDE_p$y
  hat_q <- KDE_q$y
  grid_width <- KDE_p$x[2] - KDE_p$x[1]

  hat_p_q_ratio <- hat_p / hat_q
  hat_p_q_ratio[!is.finite(hat_p_q_ratio)] <- 0

  x_vector <- seq(from = -h / 2, to = h / 2, length.out = 1000)
  eta_vector <- seq(from = 0, to = eta_max, length.out = 1000)

  alpha <- sapply(
    eta_vector,
    alpha_value,
    x_vector = x_vector,
    h = h,
    grid_width = grid_width,
    hat_p = hat_p,
    hat_p_q_ratio = hat_p_q_ratio
  )

  beta <- sapply(
    eta_vector,
    beta_value,
    x_vector = x_vector,
    h = h,
    grid_width = grid_width,
    hat_q = hat_q,
    hat_p_q_ratio = hat_p_q_ratio
  )

  data.frame(
    eta = eta_vector,
    alpha = alpha,
    beta = beta,
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------
# Mechanism & parameters
# --------------------------------------------------

sigma <- 1

x1 <- c(1, rep(0, 9))
x2 <- rep(0, 10)

sum_gauss <- function(x) {
  sum(x) + rnorm(1, mean = 0, sd = sigma)
}

N_values <- c(100, 300, 500, 1000,1500,2000,2500,3000,3500,4000)
h <- 0.1
eta_max <- 15
n_runs <- 250

# --------------------------------------------------
# Output layout
# --------------------------------------------------

base_dir <- "kde_results"
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------
# Parallel settings
# --------------------------------------------------

n_cores <- max(1L, parallel::detectCores() - 1L)

# --------------------------------------------------
# Run: 1000 CSVs per N, saved into separate folders
# --------------------------------------------------

for (N in N_values) {
  n_dir <- file.path(base_dir, sprintf("N_%04d", N))
  dir.create(n_dir, showWarnings = FALSE, recursive = TRUE)

  message("Starting N = ", N, " (saving to: ", n_dir, ")")

  ok <- parallel::mclapply(
    X = seq_len(n_runs),
    FUN = function(i) {
      df <- KDE_Estimator(
        eta_max = eta_max,
        Mechanism = sum_gauss,
        x1 = x1,
        x2 = x2,
        N = N,
        h = h
      )

      out_file <- file.path(n_dir, sprintf("kde_run_%04d.csv", i))
      utils::write.csv(df, out_file, row.names = FALSE)

      TRUE
    },
    mc.cores = n_cores
  )

  message("Finished N = ", N, " | wrote ", sum(unlist(ok)), " files")
}
