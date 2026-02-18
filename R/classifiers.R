# classifiers.R
# -------------------------------------------------------------------------
# Pluggable classifier strategies for the sequential audit engine.
#
# Each classifier is a list with functions:
#   - init(P, Q, cfg)    -> atk (must be non-NULL)
#   - refit(P, Q, cfg, atk) -> atk (may return NULL to mean "keep old atk")
#   - classify(z, atk)   -> logical vector
#   - tradeoff(atk)      -> data.frame (optional but recommended)
#
# IMPORTANT:
# - This file assumes `alpha_value()` and `beta_value()` exist
#   (source your KDE_estimator.R before running).
# -------------------------------------------------------------------------

make_kde_classifier <- function(
  kde_gridsize = NULL,     # NULL = auto based on sample size
  eps = 1e-12
) {
  list(
    init = function(P, Q, cfg) {
      kde_refit(
        P = P,
        Q = Q,
        h = cfg$h,
        eta_search_max = cfg$eta_search_max,
        claimed_curve = cfg$claimed_curve,
        kde_gridsize = kde_gridsize,
        eps = eps
      )
    },

    refit = function(P, Q, cfg, atk) {
      kde_refit(
        P = P,
        Q = Q,
        h = cfg$h,
        eta_search_max = cfg$eta_search_max,
        claimed_curve = cfg$claimed_curve,
        kde_gridsize = kde_gridsize,
        eps = eps
      )
    },

    classify = function(z, atk) atk$classify(z),

    tradeoff = function(atk) atk$tradeoff
  )
}

kde_refit <- function(
  P, Q, h, eta_search_max, claimed_curve,
  kde_gridsize = NULL,
  eps = 1e-12
) {
  if (!requireNamespace("KernSmooth", quietly = TRUE)) {
    stop("Need KernSmooth. Install with install.packages('KernSmooth').")
  }

  dpik <- KernSmooth::dpik
  bkde <- KernSmooth::bkde

  tmin <- min(P, Q)
  tmax <- max(P, Q)

  gs <- if (is.null(kde_gridsize)) {
    max(50L, min(1000L, length(P) %/% 10L))
  } else {
    as.integer(kde_gridsize)
  }

  bwP <- dpik(P, kernel = "normal", gridsize = gs, range.x = c(tmin, tmax))
  bwQ <- dpik(Q, kernel = "normal", gridsize = gs, range.x = c(tmin, tmax))

  KP <- bkde(P, kernel = "normal", bandwidth = bwP,
             gridsize = gs, range.x = c(tmin, tmax))
  KQ <- bkde(Q, kernel = "normal", bandwidth = bwQ,
             gridsize = gs, range.x = c(tmin, tmax))

  hat_p <- KP$y
  hat_q <- KQ$y
  gw <- KP$x[2] - KP$x[1]

  # alpha_value/beta_value expect log-likelihood ratio:
  log_lr <- log(hat_p + eps) - log(hat_q + eps)

  xv   <- seq(-h / 2, h / 2, length.out = 1000)
  etas <- seq(0, eta_search_max, length.out = 1000)

  ab <- alpha_beta_vectorized(etas, xv, h, gw, hat_p, hat_q, log_lr)
  a  <- ab$a
  b  <- ab$b

  ok <- which(a > 0 & b <= 1)
  eta_fallback <- etas[ceiling(length(etas) / 2)]

  if (!length(ok)) {
    eta_star <- eta_fallback
    used_fallback <- TRUE
  } else {
    used_fallback <- FALSE
    eta_star <- {
      d <- rep(NA_real_, length(ok))
      for (i in seq_along(ok)) {
        a0 <- a[ok][i]
        b0 <- b[ok][i]
        f45 <- function(a_) claimed_curve(a_) - (b0 + (a_ - a0))
        r <- tryCatch(stats::uniroot(f45, c(0, 1))$root,
                      error = function(e) NA_real_)
        if (!is.na(r)) {
          b_ <- claimed_curve(r)
          d[i] <- sqrt((r - a0)^2 + (b_ - b0)^2)
        }
      }
      if (all(is.na(d))) NA_real_ else etas[ok][which.max(d)]
    }

    if (is.na(eta_star) || !is.finite(eta_star) || eta_star <= 0) {
      eta_star <- eta_fallback
      used_fallback <- TRUE
    }
  }

  if (is.na(eta_star) || !is.finite(eta_star) || eta_star <= 0) {
    eta_star <- eta_fallback
    used_fallback <- TRUE
  }

  fP_raw <- stats::approxfun(KP$x, KP$y, rule = 2)
  fQ_raw <- stats::approxfun(KQ$x, KQ$y, rule = 2)

  fP <- function(z) pmax(fP_raw(z), 0)
  fQ <- function(z) pmax(fQ_raw(z), 0)

  list(
    tradeoff = data.frame(alpha = a, beta = b, eta = etas),
    eta_star = eta_star,
    used_fallback_eta = used_fallback,
    classify = function(z) {
      (log(fQ(z) + eps) - log(fP(z) + eps)) > log(eta_star)
    }
  )
}

make_parametric_gaussian_classifier <- function(
  sigma = 1,
  grid_size = 1000L,
  eps = 1e-12
) {
  list(
    init = function(P, Q, cfg) {
      parametric_gaussian_refit(
        P = P,
        Q = Q,
        h = cfg$h,
        eta_search_max = cfg$eta_search_max,
        claimed_curve = cfg$claimed_curve,
        sigma = sigma,
        grid_size = grid_size,
        eps = eps
      )
    },

    refit = function(P, Q, cfg, atk) {
      parametric_gaussian_refit(
        P = P,
        Q = Q,
        h = cfg$h,
        eta_search_max = cfg$eta_search_max,
        claimed_curve = cfg$claimed_curve,
        sigma = sigma,
        grid_size = grid_size,
        eps = eps
      )
    },

    classify = function(z, atk) atk$classify(z),

    tradeoff = function(atk) atk$tradeoff
  )
}

parametric_gaussian_refit <- function(
  P, Q, h, eta_search_max, claimed_curve,
  sigma = 1,
  grid_size = 1000L,
  eps = 1e-12
) {
  mu_P <- mean(P)
  mu_Q <- mean(Q)

  fP <- function(z) stats::dnorm(z, mean = mu_P, sd = sigma)
  fQ <- function(z) stats::dnorm(z, mean = mu_Q, sd = sigma)

  tmin <- min(P, Q)
  tmax <- max(P, Q)

  xgrid <- seq(tmin, tmax, length.out = as.integer(grid_size))
  hat_p <- fP(xgrid)
  hat_q <- fQ(xgrid)
  gw <- xgrid[2] - xgrid[1]

  # alpha_value/beta_value expect log-likelihood ratio:
  log_lr <- log(hat_p + eps) - log(hat_q + eps)

  xv   <- seq(-h / 2, h / 2, length.out = 1000)
  etas <- seq(0, eta_search_max, length.out = 1000)

  ab <- alpha_beta_vectorized(etas, xv, h, gw, hat_p, hat_q, log_lr)
  a  <- ab$a
  b  <- ab$b

  ok <- which(a > 0 & b <= 1)
  eta_fallback <- etas[ceiling(length(etas) / 2)]

  if (!length(ok)) {
    eta_star <- eta_fallback
    used_fallback <- TRUE
  } else {
    used_fallback <- FALSE
    eta_star <- {
      d <- rep(NA_real_, length(ok))
      for (i in seq_along(ok)) {
        a0 <- a[ok][i]
        b0 <- b[ok][i]
        f45 <- function(a_) claimed_curve(a_) - (b0 + (a_ - a0))
        r <- tryCatch(stats::uniroot(f45, c(0, 1))$root,
                      error = function(e) NA_real_)
        if (!is.na(r)) {
          b_ <- claimed_curve(r)
          d[i] <- sqrt((r - a0)^2 + (b_ - b0)^2)
        }
      }
      if (all(is.na(d))) NA_real_ else etas[ok][which.max(d)]
    }

    if (is.na(eta_star) || !is.finite(eta_star) || eta_star <= 0) {
      eta_star <- eta_fallback
      used_fallback <- TRUE
    }
  }

  if (is.na(eta_star) || !is.finite(eta_star) || eta_star <= 0) {
    eta_star <- eta_fallback
    used_fallback <- TRUE
  }

  list(
    tradeoff = data.frame(alpha = a, beta = b, eta = etas),
    eta_star = eta_star,
    used_fallback_eta = used_fallback,
    classify = function(z) {
      (log(fQ(z) + eps) - log(fP(z) + eps)) > log(eta_star)
    }
  )
}
