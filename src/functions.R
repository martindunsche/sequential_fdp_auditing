simulate_gaussian_sup_quantile <- function(M,
                                           alpha,
                                           sims  = 10000,
                                           k_max = 10000,
                                           eval_step = 10) {

  # Generate Gaussian increments
  z_mat <- matrix(rnorm(sims * (M + k_max)), nrow = sims)
  cumsums <- t(apply(z_mat, 1L, cumsum))

  # Only look at evaluation points
  eval_indices <- seq(M + eval_step, M + k_max, by = eval_step)

  cumsums_k <- cumsums[, eval_indices, drop = FALSE]

  denom <- sqrt(eval_indices) *
    sqrt(log(20 + eval_indices / M))

  # Normalize
  w_mat <- sweep(cumsums_k, 2L, denom, "/")

  row_maxes <- apply(w_mat, 1L, max)

  unname(stats::quantile(row_maxes, 1 - alpha, names = FALSE))
}




## --- Core: sequential audit ---
sequential_audit_simple <- function(
    Mechanism, x1, x2, M_burn, h, eta_search_max, claimed_curve, q_alpha,
    max_iter = 10000L, eval_step = 10L, refresh_every = 100L
) {
  if (!requireNamespace("KernSmooth", quietly = TRUE)) {
    stop("Package 'KernSmooth' is required. Install with install.packages('KernSmooth').")
  }
  
  ## Fit classifier from accumulated samples (no resampling)
  refit <- function(P, Q) {
    tmin <- min(P, Q)
    tmax <- max(P, Q)
    gs   <- min(1000L, length(P) %/% 10L)

    bwP <- KernSmooth::dpik(P, kernel = "normal", gridsize = gs, range.x = c(tmin, tmax))
    bwQ <- KernSmooth::dpik(Q, kernel = "normal", gridsize = gs, range.x = c(tmin, tmax))

    KP  <- KernSmooth::bkde(P, kernel = "normal", bandwidth = bwP, gridsize = gs, range.x = c(tmin, tmax))
    KQ  <- KernSmooth::bkde(Q, kernel = "normal", bandwidth = bwQ, gridsize = gs, range.x = c(tmin, tmax))
    
    fP  <- stats::approxfun(KP$x, KP$y, rule = 2, ties = mean)
    fQ  <- stats::approxfun(KQ$x, KQ$y, rule = 2, ties = mean)
    
    # grid for alpha/beta via your alpha_value/beta_value
    hat_p <- KP$y; hat_q <- KQ$y; gw <- KP$x[2] - KP$x[1]
    ratio <- hat_p / hat_q; ratio[!is.finite(ratio)] <- 0
    
    xv   <- seq(-h/2, h/2, length.out = 1000)
    etas <- seq(0, eta_search_max, length.out = 1000)
    
    a <- sapply(etas, alpha_value, xv, h, gw, hat_p, ratio)
    b <- sapply(etas, beta_value,  xv, h, gw, hat_q, ratio)
    
    ok <- which(a > 0 & b <= 1)
    if (!length(ok)) return(NULL)
    
    eta_star <- {
      dists <- numeric(length(ok))
      for (i in seq_along(ok)) {
        a_hat <- a[ok][i]
        b_hat <- b[ok][i]
        f45 <- function(a_prime) claimed_curve(a_prime) - (b_hat + (a_prime - a_hat))
        root <- tryCatch(uniroot(f45, c(0, 1))$root, error = function(e) NA)
        if (!is.na(root)) {
          b_prime <- claimed_curve(root)
          dists[i] <- sqrt((root - a_hat)^2 + (b_prime - b_hat)^2)
        } else {
          dists[i] <- NA
        }
      }
      if (all(is.na(dists))) NA else etas[ok][which.max(dists)]
    }

    list(
      bwP = bwP, bwQ = bwQ,      # <<< CHANGED: return bandwidths
      classify = function(z) {
        fPz <- fP(z); fQz <- fQ(z)
        as.numeric((fQz + 1e-12) / (pmax(fPz, 1e-12)) > eta_star)
      }
    )
  }
  
  ## Full-history evaluation
  eval_full <- function(P, Q, classify, q_alpha, M_burn, n) {
    e1 <- classify(P); e2 <- classify(Q)
    T1 <- mean(e1)
    T2 <- 1 - mean(e2)
    
    s1 <- sqrt(T1*(1-T1))
    s2 <- sqrt(T2*(1-T2))
      
    bump <- sqrt(log(20 + n / M_burn)) / sqrt(n)
    T1_ <- min(1, max(0, T1 + q_alpha * s1 * bump))
    T2_ <- min(1, max(0, T2 + q_alpha * s2 * bump))
    
    if (T1_ <= 0.01) T2_ <- max(0.99, T2_)
    if (T1_ >= 0.99) T2_ <- max(T2_,0.01)
    
    beta_allowed <- claimed_curve(T1_)
    c(T1_ = as.numeric(T1_), T2_ = as.numeric(T2_), beta_allowed = as.numeric(beta_allowed))
  }
  
  ## Burn-in
  P <- replicate(M_burn, Mechanism(x1))
  Q <- replicate(M_burn, Mechanism(x2))
  n <- M_burn
  
  atk <- refit(P, Q)
  bwP <- atk$bwP                         # <<< CHANGED
  bwQ <- atk$bwQ                         # <<< CHANGED

    if (is.null(atk)) {
    return(list(
      violation    = NA_real_,
      stopped_at_n = as.numeric(n),
      reason       = "could not build initial classifier (no alpha>0)",
      T1_          = NA_real_,
      T2_          = NA_real_,
      beta_allowed = NA_real_
    ))
  }
  
  
  comp <- eval_full(P, Q, atk$classify, q_alpha, M_burn, n)

   if (comp["T2_"] < comp["beta_allowed"]) {
    return(list(
      violation    = 1.0,
      stopped_at_n = as.numeric(n),
      T1_          = comp["T1_"],
      T2_          = comp["T2_"],
      beta_allowed = comp["beta_allowed"]
    ))
  }
  ## Sequential loop
  total_max_n <- n + max_iter
 n_last_kde <- n   # <<< initialize after burn-in

while (n < total_max_n) {

    k <- min(eval_step, total_max_n - n)

    # append new samples
    P <- c(P, replicate(k, Mechanism(x1)))
    Q <- c(Q, replicate(k, Mechanism(x2)))

    n <- n + k     # update current sample size

    # --- shrink rule based on last KDE update ---
    shrink <- 1 - (n_last_kde / n)^(1/5)
    trigger <- (shrink > 0.10)

    if (trigger) {
        tmp <- refit(P, Q)
        if (!is.null(tmp)) {
            atk <- tmp
            bwP <- atk$bwP
            bwQ <- atk$bwQ
            n_last_kde <- n    # <<< update marker
        }
    }

    comp <- eval_full(P, Q, atk$classify, q_alpha, M_burn, n)

    if (comp["T2_"] < comp["beta_allowed"]) {
        return(list(
            violation    = 1.0,
            stopped_at_n = as.numeric(n),
            T1_          = comp["T1_"],
            T2_          = comp["T2_"],
            beta_allowed = comp["beta_allowed"]
        ))
    }
}

  list(
    violation    = 0.0,
    stopped_at_n = as.numeric(n),
    T1_          = comp["T1_"],
    T2_          = comp["T2_"],
    beta_allowed = comp["beta_allowed"],
    reason       = "max_iter reached without violation"
  )
}



`%||%` <- function(x, y) if (length(x)) x else y


run_many_sequential_audits <- function(
    n_trials, Mechanism, x1, x2, M_burn, h, eta_search_max, claimed_curve,
    q_alpha, max_iter = 10000L, eval_step = 100L, refresh_every = 1000L,
    mc_cores = 125
) {
  
  # robust worker: never returns NULL/atomic
  results_list <- parallel::mclapply(
    seq_len(n_trials),
    function(i) {
      tryCatch(
        sequential_audit_simple(
          Mechanism, x1, x2, M_burn, h, eta_search_max, claimed_curve, q_alpha,
          max_iter, eval_step, refresh_every
        ),
        error = function(e) {
          message(sprintf("[Trial %d] Error: %s", i, conditionMessage(e)))
          list(
            violation    = NA_real_,
            stopped_at_n = NA_real_,
            T1_          = NA_real_,
            T2_          = NA_real_,
            beta_allowed = NA_real_,
            reason       = conditionMessage(e)
          )
        }
      )
    },
    mc.cores = mc_cores
  )
  
  violations <- vapply(results_list, function(z) z$violation    %||% NA_real_, numeric(1))
  stops      <- vapply(results_list, function(z) z$stopped_at_n %||% NA_real_, numeric(1))
  T1_last    <- vapply(results_list, function(z) z$T1_          %||% NA_real_, numeric(1))
  T2_last    <- vapply(results_list, function(z) z$T2_          %||% NA_real_, numeric(1))
  beta_last  <- vapply(results_list, function(z) z$beta_allowed %||% NA_real_, numeric(1))
  
  violated_indices <- which(violations == 1)
  violated_runs    <- if (length(violated_indices)) results_list[violated_indices] else list()
  
  # compact table of just the violated trials
  violated_table <- if (length(violated_indices)) {
    data.frame(
      trial        = violated_indices,
      stopped_at_n = stops[violated_indices],
      T1_          = T1_last[violated_indices],
      T2_          = T2_last[violated_indices],
      beta_allowed = beta_last[violated_indices],
      row.names    = NULL
    )
  } else {
    data.frame(
      trial        = integer(0),
      stopped_at_n = numeric(0),
      T1_          = numeric(0),
      T2_          = numeric(0),
      beta_allowed = numeric(0)
    )
  }
  
  list(
    violation_rate  = mean(violations == 1, na.rm = TRUE),
    avg_stop_n      = mean(stops, na.rm = TRUE),
    last_T1         = T1_last,
    last_T2         = T2_last,
    stops           = stops,
    raw             = results_list,
    violated_indices = violated_indices,
    violated_runs    = violated_runs,
    violated_table   = violated_table
  )
}

run_experiment <- function(
  Mechanism, x1, x2,
  M_burn = 100,
  h = 0.1,
  alpha = 0.05,
  trials = 250,
  eta_search_max = 15,
  max_iter = 10000,
  eval_step = 10,
  refresh_every = 100
){
  
  q_alpha <- simulate_gaussian_sup_quantile(
    M = M_burn,
    alpha = alpha / 2,
    sims = 10000,
    k_max = 10000
  )

  res <- run_many_sequential_audits(
    n_trials = trials,
    Mechanism = Mechanism,
    x1 = x1,
    x2 = x2,
    M_burn = M_burn,
    h = h,
    eta_search_max = eta_search_max,
    claimed_curve = claimed_curve,
    q_alpha = q_alpha,
    max_iter = max_iter,
    eval_step = eval_step,
    refresh_every = refresh_every
  )

  cat("Violation rate:", res$violation_rate, "\n")
  cat("Average stopping n:", res$avg_stop_n, "\n")

  return(res)
}

