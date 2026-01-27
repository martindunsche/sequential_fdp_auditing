# KDE_estimator.R
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
