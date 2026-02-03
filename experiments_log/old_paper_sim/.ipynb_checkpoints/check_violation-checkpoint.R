#!/usr/bin/env Rscript

# ---------- args ----------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
mu_txt <- get_arg("--mu", NA_character_)
if (is.na(mu_txt)) stop("Usage: ./check_violation.R --mu <number>")
mu <- as.numeric(mu_txt)
if (!is.finite(mu)) stop("--mu must be a finite number")

# ---------- helpers ----------
claimed_f <- function(alpha) pnorm(qnorm(1 - alpha) - mu)

get_n_from_path <- function(path) {
  m <- gregexpr("(?i)(^|[^A-Za-z0-9])N_?([0-9]+)", path, perl = TRUE)
  hits <- regmatches(path, m)[[1]]
  if (length(hits) == 0) stop("Couldn't infer n: no 'N<digits>' or 'N_<digits>' found in path: ", path)
  ns <- as.integer(sub(".*(?i)N_?([0-9]+).*", "\\1", hits, perl = TRUE))
  ns[length(ns)]
}

# ---------- main ----------
base_dir <- "kde_results"
files <- list.files(base_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)

audit_one <- function(f) {
  n <- get_n_from_path(f)
  omega <- sqrt(log(4/0.05) / (2*n))
  message(sprintf("[check_violation] %s : n=%d omega=%g mu=%g", f, n, omega, mu))

  d <- read.csv(f)

  # avoid NaNs: need 0 < alpha and alpha+omega < 1
  a <- d$alpha
  keep <- is.finite(a) & is.finite(d$beta) & (a > 0) & (a + omega < 1)
  d <- d[keep, , drop = FALSE]
  if (nrow(d) == 0) return(NA_integer_)

 a <- d$alpha
b <- d$beta

ok <- which(is.finite(a) & is.finite(b) & (a > 0) & (a + omega < 1))
if (!length(ok)) return(NA_integer_)

# distance to intersection of claimed curve with 45-degree line through (a0,b0)
dist45 <- rep(NA_real_, length(ok))

for (k in seq_along(ok)) {
  j <- ok[k]
  a0 <- a[j]
  b0 <- b[j]

  f45 <- function(a_) claimed_f(a_) - (b0 + (a_ - a0))

  # restrict root search to domain where claimed_f is valid
  lo <- 1e-12
  hi <- 1 - omega - 1e-12

  r <- tryCatch(stats::uniroot(f45, c(lo, hi))$root,
                error = function(e) NA_real_)

  if (!is.na(r) && is.finite(r)) {
    b_ <- claimed_f(r)
    dist45[k] <- sqrt((r - a0)^2 + (b_ - b0)^2)
  }
}

# pick the point with maximum 45-degree distance
i <- if (all(is.na(dist45))) NA_integer_ else ok[which.max(dist45)]

# fallback like you want
if (length(i) == 0 || is.na(i)) i <- 500L
i <- max(1L, min(as.integer(i), nrow(d)))  # clamp to valid row


  as.integer(d$beta[i] + omega < claimed_f(d$alpha[i] + omega))
}

out <- data.frame(
  file = files,
  n = vapply(files, get_n_from_path, integer(1)),
  violation = vapply(files, audit_one, integer(1)),
  stringsAsFactors = FALSE
)

mu_tag <- gsub("[^0-9A-Za-z._-]+", "_", sprintf("%g", mu))
out_name <- sprintf("audit_results_mu%s.csv", mu_tag)
write.csv(out, out_name, row.names = FALSE)
message(sprintf("[check_violation] wrote %s", out_name))
