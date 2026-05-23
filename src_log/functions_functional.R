# functions_functional.R
# -------------------------------------------------------------------------
# Functional sequential audit engine with pluggable classifiers.
#
# This file contains ONLY the audit engine and helpers.
# Classifiers live in classifiers.R.
#
# Required from outside:
# - claimed_curve: function(alpha) -> beta_allowed
# - classifier: list(init/refit/classify/tradeoff)
#
# IMPORTANT: alpha_value/beta_value are not used here directly,
# but your classifiers (KDE/parametric) will use them.
# -------------------------------------------------------------------------

simulate_gaussian_sup_quantile <- function(
  M,
  alpha,
  sims = 10000,
  k_max = 10000,
  eval_step = 10
) {
  z_mat <- matrix(stats::rnorm(sims * (M + k_max)), nrow = sims)
  cumsums <- t(apply(z_mat, 1L, cumsum))

  eval_indices <- seq(M + eval_step, M + k_max, by = eval_step)
  cumsums_k <- cumsums[, eval_indices, drop = FALSE]

  denom <- sqrt(eval_indices) * sqrt(log(20 + eval_indices / M))
  w_mat <- sweep(cumsums_k, 2L, denom, "/")

  row_maxes <- apply(w_mat, 1L, max)
  unname(stats::quantile(row_maxes, 1 - alpha, names = FALSE))
}

`%||%` <- function(x, y) if (!is.null(x) && length(x) == 1L) x else y

make_empty_result <- function(reason = "error") {
  list(
    violation    = NA_real_,
    stopped_at_n = NA_real_,
    T1_          = NA_real_,
    T2_          = NA_real_,
    s1           = NA_real_,
    s2           = NA_real_,
    beta_allowed = NA_real_,
    comp         = NA,
    tradeoff     = NA,
    reason       = reason
  )
}

pack_out <- function(violation, stopped_at_n, comp, tradeoff, reason = NULL) {
  out <- list(
    violation    = as.numeric(violation),
    stopped_at_n = as.numeric(stopped_at_n),
    T1_          = unname(comp["T1_"]),
    T2_          = unname(comp["T2_"]),
    s1           = unname(comp["s1"]),
    s2           = unname(comp["s2"]),
    beta_allowed = unname(comp["beta_allowed"]),
    comp         = comp,
    tradeoff     = tradeoff
  )
  if (!is.null(reason)) out$reason <- reason
  out
}

validate_classifier <- function(classifier) {
  needed <- c("init", "refit", "classify")
  missing <- setdiff(needed, names(classifier))
  if (length(missing) > 0) {
    stop("`classifier` is missing: ", paste(missing, collapse = ", "))
  }
  invisible(classifier)
}

audit_eval_full <- function(
  P, Q, classify, q_alpha, M_burn, n, claimed_curve,
  boundary_tweaks = FALSE
) {
  e1 <- classify(P)
  e2 <- classify(Q)

  # If classify ever returns NA, fail fast (prevents silent NA propagation)
  if (anyNA(e1) || anyNA(e2)) {
    stop("classify() returned NA values")
  }

  T1 <- mean(e1)
  T2 <- 1 - mean(e2)

  # Variance floors so bounds don't collapse when T1 or T2 hit {0,1}
  # (prevents pathological early decisions when classifier saturates)
  s1 <- max(sqrt(T1 * (1 - T1)), 1 / (2 * sqrt(n)))
  s2 <- max(sqrt(T2 * (1 - T2)), 1 / (2 * sqrt(n)))

  bump <- sqrt(log(20 + n / M_burn)) / sqrt(n)

  T1_ <- min(1, max(0, T1 + q_alpha * s1 * bump))
  T2_ <- min(1, max(0, T2 + q_alpha * s2 * bump))

  if (isTRUE(boundary_tweaks)) {
    if (T1_ <= 0.01) T2_ <- max(0.99, T2_)
    if (T1_ >= 0.99) T2_ <- max(T2_, 0.01)
  }

  c(
    T1 = as.numeric(T1),
    T2 = as.numeric(T2),
    s1 = as.numeric(s1),
    s2 = as.numeric(s2),
    T1_ = as.numeric(T1_),
    T2_ = as.numeric(T2_),
    beta_allowed = as.numeric(claimed_curve(T1_))
  )
}

audit_violation <- function(comp) {
  isTRUE(comp["T2_"] < comp["beta_allowed"])
}

audit_append_samples <- function(P, Q, Mechanism, x1, x2, k) {
  list(
    P = c(P, replicate(k, Mechanism(x1))),
    Q = c(Q, replicate(k, Mechanism(x2)))
  )
}

audit_should_refit <- function(n_last_refit, n, threshold = 0.10) {
  shrink <- 1 - (n_last_refit / n)^(1 / 5)
  list(trigger = isTRUE(shrink > threshold), shrink = shrink)
}

audit_step <- function(state, cfg) {
  k <- min(cfg$eval_step, cfg$total_max_n - state$n)

  samp <- audit_append_samples(
    P = state$P,
    Q = state$Q,
    Mechanism = cfg$Mechanism,
    x1 = cfg$x1,
    x2 = cfg$x2,
    k = k
  )

  n_new <- state$n + k

  refit_dec <- audit_should_refit(
    n_last_refit = state$n_last_refit,
    n = n_new,
    threshold = cfg$refit_threshold
  )

  atk_new <- state$atk
  n_last_refit_new <- state$n_last_refit

  if (refit_dec$trigger) {
    atk_try <- cfg$classifier$refit(
      P = samp$P,
      Q = samp$Q,
      cfg = cfg,
      atk = atk_new
    )
    if (!is.null(atk_try)) {
      atk_new <- atk_try
      n_last_refit_new <- n_new
    }
  }

  classify_fun <- function(z) cfg$classifier$classify(z, atk_new)

  comp_new <- audit_eval_full(
    P = samp$P,
    Q = samp$Q,
    classify = classify_fun,
    q_alpha = cfg$q_alpha,
    M_burn = cfg$M_burn,
    n = n_new,
    claimed_curve = cfg$claimed_curve,
    boundary_tweaks = cfg$boundary_tweaks
  )

  list(
    P = samp$P,
    Q = samp$Q,
    n = n_new,
    atk = atk_new,
    n_last_refit = n_last_refit_new,
    comp = comp_new,
    shrink = refit_dec$shrink
  )
}

sequential_audit_simple <- function(
  Mechanism, x1, x2,
  M_burn, h, eta_search_max,
  claimed_curve, q_alpha,
  classifier,
  max_iter = 10000L,
  eval_step = 10L,
  refit_threshold = 0.10,
  boundary_tweaks = FALSE
) {
  validate_classifier(classifier)

  P0 <- replicate(M_burn, Mechanism(x1))
  Q0 <- replicate(M_burn, Mechanism(x2))
  n0 <- M_burn

  cfg <- list(
    Mechanism = Mechanism,
    x1 = x1,
    x2 = x2,
    M_burn = M_burn,
    h = h,
    eta_search_max = eta_search_max,
    claimed_curve = claimed_curve,
    q_alpha = q_alpha,
    eval_step = eval_step,
    total_max_n = n0 + max_iter,
    refit_threshold = refit_threshold,
    classifier = classifier,
    boundary_tweaks = boundary_tweaks
  )

  atk0 <- classifier$init(P = P0, Q = Q0, cfg = cfg)
  if (is.null(atk0)) {
    return(make_empty_result("classifier init returned NULL at burn-in"))
  }

  comp0 <- audit_eval_full(
    P = P0,
    Q = Q0,
    classify = function(z) classifier$classify(z, atk0),
    q_alpha = q_alpha,
    M_burn = M_burn,
    n = n0,
    claimed_curve = claimed_curve,
    boundary_tweaks = boundary_tweaks
  )

  tradeoff0 <- if (!is.null(classifier$tradeoff)) classifier$tradeoff(atk0) else NA

  if (audit_violation(comp0)) {
    return(pack_out(1, n0, comp0, tradeoff0, reason = "violation at burn-in"))
  }

  state <- list(
    P = P0,
    Q = Q0,
    n = n0,
    atk = atk0,
    n_last_refit = n0,
    comp = comp0
  )

  while (state$n < cfg$total_max_n) {
    state <- audit_step(state, cfg)
    if (audit_violation(state$comp)) {
      tradeoff <- if (!is.null(classifier$tradeoff)) classifier$tradeoff(state$atk) else NA
      return(pack_out(1, state$n, state$comp, tradeoff, reason = "violation"))
    }
  }

  tradeoff <- if (!is.null(classifier$tradeoff)) classifier$tradeoff(state$atk) else NA
  pack_out(0, state$n, state$comp, tradeoff, reason = "max_iter reached without violation")
}

run_many_sequential_audits <- function(
  n_trials,
  Mechanism, x1, x2,
  M_burn, h, eta_search_max,
  claimed_curve, q_alpha,
  classifier,
  max_iter = 10000L,
  eval_step = 100L,
  mc_cores = 8L,
  refit_threshold = 0.10,
  boundary_tweaks = FALSE,
  show_progress = FALSE
) {
  validate_classifier(classifier)

  run_one_trial <- function(i) {
    tryCatch(
      sequential_audit_simple(
        Mechanism = Mechanism,
        x1 = x1,
        x2 = x2,
        M_burn = M_burn,
        h = h,
        eta_search_max = eta_search_max,
        claimed_curve = claimed_curve,
        q_alpha = q_alpha,
        classifier = classifier,
        max_iter = max_iter,
        eval_step = eval_step,
        refit_threshold = refit_threshold,
        boundary_tweaks = boundary_tweaks
      ),
      error = function(e) {
        message(sprintf("[Trial %d] Error: %s", i, conditionMessage(e)))
        make_empty_result(conditionMessage(e))
      }
    )
  }

  if (isTRUE(show_progress)) {
    if (!requireNamespace("pbmcapply", quietly = TRUE)) {
      stop("show_progress=TRUE requires the pbmcapply package.")
    }
    results_list <- pbmcapply::pbmclapply(
      seq_len(n_trials),
      run_one_trial,
      mc.cores = mc_cores,
      ignore.interactive = TRUE
    )
  } else {
    results_list <- parallel::mclapply(
      seq_len(n_trials),
      run_one_trial,
      mc.cores = mc_cores
    )
  }

  get1 <- function(z, nm) z[[nm]] %||% NA_real_

  violations <- vapply(results_list, get1, numeric(1), nm = "violation")
  stops      <- vapply(results_list, get1, numeric(1), nm = "stopped_at_n")
  T1_last    <- vapply(results_list, get1, numeric(1), nm = "T1_")
  T2_last    <- vapply(results_list, get1, numeric(1), nm = "T2_")
  s1_last    <- vapply(results_list, get1, numeric(1), nm = "s1")
  s2_last    <- vapply(results_list, get1, numeric(1), nm = "s2")
  beta_last  <- vapply(results_list, get1, numeric(1), nm = "beta_allowed")

  violated_indices <- which(violations == 1)

  violated_table <- if (length(violated_indices)) {
    data.frame(
      trial        = violated_indices,
      stopped_at_n = stops[violated_indices],
      T1_          = T1_last[violated_indices],
      T2_          = T2_last[violated_indices],
      s1           = s1_last[violated_indices],
      s2           = s2_last[violated_indices],
      beta_allowed = beta_last[violated_indices],
      row.names    = NULL
    )
  } else {
    data.frame(
      trial        = integer(0),
      stopped_at_n = numeric(0),
      T1_          = numeric(0),
      T2_          = numeric(0),
      s1           = numeric(0),
      s2           = numeric(0),
      beta_allowed = numeric(0)
    )
  }

  list(
    violation_rate   = mean(violations == 1, na.rm = TRUE),
    avg_stop_n       = mean(stops, na.rm = TRUE),
    last_T1          = T1_last,
    last_T2          = T2_last,
    last_s1          = s1_last,
    last_s2          = s2_last,
    beta_allowed     = beta_last,
    stops            = stops,
    raw              = results_list,
    violated_indices = violated_indices,
    violated_table   = violated_table
  )
}

run_experiment <- function(
  Mechanism, x1, x2,
  claimed_curve,
  classifier,
  M_burn = 100,
  h = 0.1,
  alpha = 0.05,
  trials = 250,
  eta_search_max = 15,
  max_iter = 10000,
  eval_step = 10,
  mc_cores = 8L,
  refit_threshold = 0.10,
  boundary_tweaks = FALSE,
  show_progress = FALSE,
  qalpha_sims = 10000,
  qalpha_kmax = 10000
) {
  validate_classifier(classifier)

  q_alpha <- simulate_gaussian_sup_quantile(
    M = M_burn,
    alpha = alpha / 2,
    sims = qalpha_sims,
    k_max = qalpha_kmax
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
    classifier = classifier,
    max_iter = max_iter,
    eval_step = eval_step,
    mc_cores = mc_cores,
    refit_threshold = refit_threshold,
    boundary_tweaks = boundary_tweaks,
    show_progress = show_progress
  )

  cat("Violation rate:", res$violation_rate, "\n")
  cat("Average stopping n:", res$avg_stop_n, "\n")
  res
}
