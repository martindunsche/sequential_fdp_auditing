# kde_estimator.R
# ------------------------------------------------------------
# Integral helpers used by the KDE-based classifier.
# IMPORTANT: these functions assume `log_lr` is log(p) - log(q)
# ------------------------------------------------------------

inner_int <- function(x, eta, grid_width, p, log_lr) {
  thresh <- eta + x
  if (thresh <= 0) return(sum(p * grid_width))

  log_thresh <- log(thresh)

  p_vec <- p
  p_vec[log_lr <= log_thresh] <- 0

  sum(p_vec * grid_width)
}

alpha_value <- function(eta, x_vector, h, grid_width, hat_p, log_lr) {
  dx <- x_vector[2] - x_vector[1]

  val <- sum(
    vapply(
      x_vector,
      inner_int,
      numeric(1),
      eta = eta,
      grid_width = grid_width,
      p = hat_p,
      log_lr = log_lr
    ) * dx / h
  )

  1 - val
}

beta_value <- function(eta, x_vector, h, grid_width, hat_q, log_lr) {
  dx <- x_vector[2] - x_vector[1]

  sum(
    vapply(
      x_vector,
      inner_int,
      numeric(1),
      eta = eta,
      grid_width = grid_width,
      p = hat_q,
      log_lr = log_lr
    ) * dx / h
  )
}

# ------------------------------------------------------------
# Vectorized replacement for the sapply(etas, alpha_value, ...)
# and sapply(etas, beta_value, ...) double loops.
#
# Instead of 1000 etas x 1000 xv points = 1M inner_int() calls,
# this precomputes a cumulative sum of (hat_p * gw) sorted by
# log_lr, then uses a single vectorized findInterval() call over
# all (eta, x) threshold combinations.
#
# Returns list(a = alpha_vector, b = beta_vector) for all etas.
# ------------------------------------------------------------
alpha_beta_vectorized <- function(etas, xv, h, gw, hat_p, hat_q, log_lr) {
  dx <- xv[2] - xv[1]

  # Sort grid points by log-likelihood ratio (ascending)
  ord          <- order(log_lr)
  log_lr_s     <- log_lr[ord]

  # Cumulative sums from the right in sorted order:
  #   cs[i] = sum of (density * gw) for all grid points j with log_lr[j] >= log_lr_s[i]
  # Appending a trailing 0 handles the "above all thresholds" edge case.
  csP <- c(rev(cumsum(rev((hat_p * gw)[ord]))), 0)
  csQ <- c(rev(cumsum(rev((hat_q * gw)[ord]))), 0)

  # Build all (eta, x) threshold combinations in one outer product (1000x1000 matrix)
  thresh_mat <- outer(etas, xv, "+")   # nEta x nX

  # log-threshold: -Inf where thresh <= 0 (whole mass is included).
  # Compute log() only on positive entries to avoid NaN warnings.
  log_thresh_mat <- thresh_mat          # reuse storage
  pos <- thresh_mat > 0
  log_thresh_mat[ pos] <- log(thresh_mat[pos])
  log_thresh_mat[!pos] <- -Inf

  # For each threshold t, we want sum(density * gw where log_lr > t).
  # findInterval(t, log_lr_s) returns the largest index i s.t. log_lr_s[i] <= t,
  # so the desired sum is csP[idx + 1]  (elements at positions idx+1 .. n).
  # idx = 0  → t below all values → include everything → csP[1] = total mass
  # idx = n  → t above all values → include nothing  → csP[n+1] = 0 (appended)
  idx <- findInterval(as.vector(log_thresh_mat), log_lr_s)  # length nEta*nX

  Sp_mat <- matrix(csP[idx + 1L], nrow = length(etas))  # nEta x nX
  Sq_mat <- matrix(csQ[idx + 1L], nrow = length(etas))  # nEta x nX

  # alpha(eta) = 1 - (dx/h) * sum_x S_p(log(eta+x))
  # beta(eta)  =     (dx/h) * sum_x S_q(log(eta+x))
  a <- 1 - rowSums(Sp_mat) * (dx / h)
  b <-     rowSums(Sq_mat) * (dx / h)

  list(a = a, b = b)
}
