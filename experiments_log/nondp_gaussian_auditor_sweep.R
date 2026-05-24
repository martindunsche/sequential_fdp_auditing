#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# NonDPGaussian1 sequential f-DP auditor sweep, both orientations
# -----------------------------------------------------------------------------
# Purpose:
#   For each claimed epsilon, instantiate the corresponding NonDPGaussian1
#   benchmark mechanism, run the sequential f-DP auditor in both neighboring
#   orientations D -> D' and D' -> D, save one per-epsilon
#   stopping-time CSV with resume support, then summarize and draw two figures:
#     1) average stopping time to rejection vs epsilon,
#     2) empirical rejection rate vs epsilon.
#
# Privacy object audited:
#   The f-DP lower bound induced by the claimed (epsilon, delta)-DP statement,
#   with delta = 1e-5 by default:
#     f_{eps,delta}(alpha) = max{0,
#       1 - delta - exp(eps) alpha,
#       exp(-eps) (1 - delta - alpha)}.
#
# The mechanism parameter and the claimed epsilon are tied together.  That is,
# at epsilon = e, this script audits M_e against the claim (e, delta)-DP.
# Since DP is directional, the script audits both orientations and reports the
# combined result: rejection if either orientation rejects.  By default the
# sequential critical value uses a Bonferroni split over the two orientations.
# -----------------------------------------------------------------------------

`%||%` <- function(x, y) if (!is.null(x) && length(x) == 1L) x else y

get_flag_value <- function(args, flag, default = NULL) {
  key <- paste0(flag, "=")
  hit <- grep(paste0("^", key), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^", key), "", hit[1L])
}

get_flag_bool <- function(args, flag, default = FALSE) {
  val <- get_flag_value(args, flag, NULL)
  if (is.null(val)) return(default)
  tolower(val) %in% c("1", "true", "t", "yes", "y")
}

parse_eps_grid <- function(s) {
  if (is.null(s) || !nzchar(s)) {
    return(c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0))
  }
  s <- gsub("\\s+", "", s)
  if (grepl(",", s, fixed = TRUE)) {
    out <- as.numeric(strsplit(s, ",", fixed = TRUE)[[1L]])
  } else if (grepl(":", s, fixed = TRUE)) {
    parts <- as.numeric(strsplit(s, ":", fixed = TRUE)[[1L]])
    if (length(parts) != 3L || anyNA(parts)) {
      stop("--eps-grid with ':' must have form start:end:length, e.g. 0.01:3:80")
    }
    out <- seq(parts[1L], parts[2L], length.out = as.integer(parts[3L]))
  } else {
    out <- as.numeric(s)
  }
  if (anyNA(out) || any(out <= 0)) stop("All epsilon values must be positive numerics")
  unique(out)
}

fmt_eps <- function(eps) {
  # Stable, readable tags such as 0.01, 0.1, 1, 2.5.
  tag <- formatC(eps, format = "fg", digits = 10)
  tag <- sub("\\.$", "", tag)
  tag
}

safe_num <- function(x) {
  as.numeric(x %||% NA_real_)
}

find_source_dir <- function(src_dir_arg) {
  candidates <- unique(c(
    src_dir_arg,
    ".",
    "./src_log",
    "../src_log",
    dirname(normalizePath(sys.frame(1)$ofile %||% ".", mustWork = FALSE))
  ))
  needed <- c("KDE_estimator.R", "classifiers.R", "functions_functional.R")
  for (d in candidates) {
    if (is.null(d) || !nzchar(d)) next
    paths <- file.path(d, needed)
    if (all(file.exists(paths))) return(normalizePath(d, mustWork = TRUE))
  }
  stop(
    "Could not find source files. Use --src-dir=/path/to/dir containing ",
    paste(needed, collapse = ", ")
  )
}

source_audit_files <- function(src_dir) {
  # Order matters: classifiers.R expects alpha_value/beta_value at runtime.
  source(file.path(src_dir, "KDE_estimator.R"))
  source(file.path(src_dir, "classifiers.R"))
  source(file.path(src_dir, "functions_functional.R"))
  mech_path <- file.path(src_dir, "mechanisms.R")
  if (file.exists(mech_path)) source(mech_path)
}

make_nonDP_gaussian <- function(eps) {
  # NonDPGaussian1 from the comparison benchmark: true sample size n is used
  # in both the mean and the noise scale.
  function(x) {
    n <- length(x)
    mean_nonpriv <- sum(x) / n
    rho <- stats::rnorm(1L, mean = 0, sd = 2 / (n * eps))
    mean_nonpriv + rho
  }
}

make_approxdp_tradeoff <- function(eps, delta = 1e-5) {
  force(eps)
  force(delta)
  function(alpha) {
    pmax(
      0,
      1 - delta - exp(eps) * alpha,
      exp(-eps) * (1 - delta - alpha)
    )
  }
}

get_reason1 <- function(z) {
  cand <- list(
    z$reason,
    z$stop_reason,
    z$reasons,
    z$meta$reason,
    z$info$reason
  )
  for (r in cand) {
    if (is.null(r) || length(r) == 0L) next
    r <- as.character(r[1L])
    if (!is.na(r) && nzchar(r)) return(r)
  }
  "OK"
}

normalize_violation <- function(z) {
  v <- safe_num(z$violation)
  if (is.na(v)) {
    reason <- get_reason1(z)
    v <- as.numeric(reason %in% c("violation", "violation at burn-in"))
  }
  v
}

result_margin <- function(z) {
  safe_num(z$beta_allowed) - safe_num(z$T2_)
}

prefixed_result_df <- function(z, prefix) {
  reason <- get_reason1(z)
  violation <- normalize_violation(z)
  stopped <- safe_num(z$stopped_at_n)
  T1_ <- safe_num(z$T1_)
  T2_ <- safe_num(z$T2_)
  beta_allowed <- safe_num(z$beta_allowed)

  out <- data.frame(
    violation    = as.numeric(violation),
    stopped_at_n = stopped,
    T1_          = T1_,
    T2_          = T2_,
    s1           = safe_num(z$s1),
    s2           = safe_num(z$s2),
    beta_allowed = beta_allowed,
    margin       = beta_allowed - T2_,
    reason       = as.character(reason),
    row.names    = NULL
  )
  names(out) <- paste0(prefix, "_", names(out))
  out
}

select_direction_result <- function(forward, reverse) {
  vf <- normalize_violation(forward)
  vr <- normalize_violation(reverse)
  sf <- safe_num(forward$stopped_at_n)
  sr <- safe_num(reverse$stopped_at_n)
  mf <- result_margin(forward)
  mr <- result_margin(reverse)

  # Combined decision: reject if either orientation rejects.  If both reject,
  # report the one that rejects first; if tied, report the larger margin
  # beta_allowed - T2_.  Without rejection, report the direction with larger
  # final margin, only for diagnostics.
  if (isTRUE(vf == 1) || isTRUE(vr == 1)) {
    candidates <- data.frame(
      direction = c("forward", "reverse"),
      violation = c(vf, vr),
      stopped_at_n = c(sf, sr),
      margin = c(mf, mr),
      stringsAsFactors = FALSE
    )
    cand_ok <- !is.na(candidates$violation) & candidates$violation == 1 & is.finite(candidates$stopped_at_n)
    candidates <- candidates[cand_ok, , drop = FALSE]
    if (!nrow(candidates)) {
      # Fallback in the unlikely case of malformed outputs.
      candidates <- data.frame(
        direction = c("forward", "reverse"),
        stopped_at_n = c(sf, sr),
        margin = c(mf, mr),
        stringsAsFactors = FALSE
      )
    }
    candidates <- candidates[order(candidates$stopped_at_n, -candidates$margin), , drop = FALSE]
    selected_direction <- candidates$direction[1L]
    selected <- if (identical(selected_direction, "forward")) forward else reverse
    selected$violation <- 1
    selected$stopped_at_n <- candidates$stopped_at_n[1L]
    selected$reason <- paste0("violation in ", selected_direction, " orientation")
    return(list(direction = selected_direction, result = selected))
  }

  # No rejection in either direction.  Select the larger final margin for
  # diagnostic columns, but keep violation = 0.
  if (is.finite(mf) || is.finite(mr)) {
    selected_direction <- ifelse(!is.finite(mr) || (is.finite(mf) && mf >= mr), "forward", "reverse")
  } else {
    selected_direction <- "forward"
  }
  selected <- if (identical(selected_direction, "forward")) forward else reverse
  selected$violation <- 0
  selected$reason <- "max_iter reached without violation in either orientation"
  list(direction = selected_direction, result = selected)
}

row_from_direction_pair <- function(pair, trial, eps) {
  forward <- pair$forward
  reverse <- pair$reverse
  sel <- select_direction_result(forward, reverse)
  z <- sel$result
  reason <- get_reason1(z)
  violation <- normalize_violation(z)

  base <- data.frame(
    trial              = as.integer(trial),
    stopped_at_n       = safe_num(z$stopped_at_n),
    selected_direction = as.character(sel$direction),
    T1_                = safe_num(z$T1_),
    T2_                = safe_num(z$T2_),
    s1                 = safe_num(z$s1),
    s2                 = safe_num(z$s2),
    beta_allowed       = safe_num(z$beta_allowed),
    margin             = result_margin(z),
    reason             = as.character(reason),
    violation          = as.numeric(violation),
    eps                = as.numeric(eps),
    row.names          = NULL
  )

  cbind(
    base,
    prefixed_result_df(forward, "fwd"),
    prefixed_result_df(reverse, "rev")
  )
}

read_existing_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  # This script writes bidirectional CSVs.  Do not resume from older
  # single-orientation CSVs even if they share a similar name.
  required <- c("trial", "stopped_at_n", "selected_direction", "fwd_stopped_at_n", "rev_stopped_at_n")
  if (!all(required %in% names(df))) return(NULL)
  df
}

write_sorted_csv <- function(df, path) {
  df <- df[order(df$trial), , drop = FALSE]
  utils::write.csv(df, path, row.names = FALSE)
}

get_completed_trials <- function(df, trials) {
  if (is.null(df) || !nrow(df)) return(integer(0L))
  # Require the bidirectional schema.  Older single-orientation CSVs should
  # not be treated as completed for this script.
  required <- c("trial", "stopped_at_n", "selected_direction", "fwd_stopped_at_n", "rev_stopped_at_n")
  if (!all(required %in% names(df))) return(integer(0L))
  ok <- !is.na(df$trial) & df$trial >= 1L & df$trial <= trials & !is.na(df$stopped_at_n)
  unique(as.integer(df$trial[ok]))
}

split_chunks <- function(x, chunk_size) {
  if (!length(x)) return(list())
  split(x, ceiling(seq_along(x) / chunk_size))
}

run_missing_trials_for_eps <- function(
  eps, eps_index, missing_trials, existing_df, outfile,
  q_alpha, classifier, x1, x2,
  M_burn, h, eta_search_max, max_iter, eval_step, refit_threshold,
  boundary_tweaks, mc_cores, show_progress, chunk_size, seed, delta_claim
) {
  Mechanism <- make_nonDP_gaussian(eps)
  claimed_curve <- make_approxdp_tradeoff(eps, delta_claim)

  run_one_direction <- function(trial_id, direction) {
    if (!is.na(seed)) {
      # Trial-level and direction-level seed makes resume independent of
      # chunking/interruption and avoids coupling the two orientations.
      dir_offset <- if (identical(direction, "forward")) 0L else 50000017L
      seed_i <- as.integer((seed + 1000003L * eps_index + 97L * trial_id + dir_offset) %% .Machine$integer.max)
      set.seed(seed_i)
    }

    xx1 <- if (identical(direction, "forward")) x1 else x2
    xx2 <- if (identical(direction, "forward")) x2 else x1

    tryCatch(
      sequential_audit_simple(
        Mechanism = Mechanism,
        x1 = xx1,
        x2 = xx2,
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
        message(sprintf("[eps=%s trial=%d direction=%s] Error: %s", fmt_eps(eps), trial_id, direction, conditionMessage(e)))
        make_empty_result(conditionMessage(e))
      }
    )
  }

  run_one <- function(trial_id) {
    list(
      forward = run_one_direction(trial_id, "forward"),
      reverse = run_one_direction(trial_id, "reverse")
    )
  }

  chunks <- split_chunks(missing_trials, chunk_size)
  df <- existing_df
  if (is.null(df)) df <- data.frame()

  for (j in seq_along(chunks)) {
    trial_chunk <- chunks[[j]]
    cat(sprintf(
      "eps=%s: running chunk %d/%d with %d trials, both orientations\n",
      fmt_eps(eps), j, length(chunks), length(trial_chunk)
    ))

    if (isTRUE(show_progress) && requireNamespace("pbmcapply", quietly = TRUE)) {
      res_list <- pbmcapply::pbmclapply(
        trial_chunk,
        run_one,
        mc.cores = mc_cores,
        ignore.interactive = TRUE
      )
    } else {
      res_list <- parallel::mclapply(
        trial_chunk,
        run_one,
        mc.cores = mc_cores,
        mc.set.seed = TRUE
      )
    }

    new_df <- do.call(
      rbind,
      Map(function(z, id) row_from_direction_pair(z, trial = id, eps = eps), res_list, trial_chunk)
    )

    if (nrow(df)) {
      # Replace any stale duplicate rows for the same trials.
      df <- df[!(df$trial %in% new_df$trial), , drop = FALSE]
      df <- rbind(df, new_df)
    } else {
      df <- new_df
    }
    write_sorted_csv(df, outfile)
    cat("Saved checkpoint:", outfile, "\n")
  }

  invisible(df)
}

summarize_one_csv <- function(path, eps, trials_expected) {
  df <- read_existing_csv(path)
  if (is.null(df)) {
    return(data.frame(
      eps = eps, n_done = 0L, n_expected = trials_expected, n_reject = 0L,
      rejection_rate = NA_real_, rejection_se = NA_real_, rejection_ci_low = NA_real_, rejection_ci_high = NA_real_,
      mean_stop_reject = NA_real_, sd_stop_reject = NA_real_, se_stop_reject = NA_real_,
      stop_reject_ci_low = NA_real_, stop_reject_ci_high = NA_real_,
      mean_stop_all = NA_real_, sd_stop_all = NA_real_, se_stop_all = NA_real_,
      csv = path,
      row.names = NULL
    ))
  }

  if ("violation" %in% names(df)) {
    rejected <- df$violation == 1
  } else {
    rejected <- df$reason %in% c("violation", "violation at burn-in")
  }
  rejected[is.na(rejected)] <- FALSE

  stops <- as.numeric(df$stopped_at_n)
  n_done <- sum(!is.na(stops))
  n_rej <- sum(rejected & !is.na(stops))
  p_hat <- if (n_done > 0L) n_rej / n_done else NA_real_
  p_se <- if (n_done > 0L && is.finite(p_hat)) sqrt(p_hat * (1 - p_hat) / n_done) else NA_real_

  stop_rej <- stops[rejected & !is.na(stops)]
  mean_rej <- if (length(stop_rej)) mean(stop_rej) else NA_real_
  sd_rej <- if (length(stop_rej) > 1L) stats::sd(stop_rej) else NA_real_
  se_rej <- if (length(stop_rej) > 1L) sd_rej / sqrt(length(stop_rej)) else NA_real_

  stop_all <- stops[!is.na(stops)]
  mean_all <- if (length(stop_all)) mean(stop_all) else NA_real_
  sd_all <- if (length(stop_all) > 1L) stats::sd(stop_all) else NA_real_
  se_all <- if (length(stop_all) > 1L) sd_all / sqrt(length(stop_all)) else NA_real_

  data.frame(
    eps = eps,
    n_done = n_done,
    n_expected = trials_expected,
    n_reject = n_rej,
    rejection_rate = p_hat,
    rejection_se = p_se,
    rejection_ci_low = if (is.finite(p_hat) && is.finite(p_se)) max(0, p_hat - 1.96 * p_se) else NA_real_,
    rejection_ci_high = if (is.finite(p_hat) && is.finite(p_se)) min(1, p_hat + 1.96 * p_se) else NA_real_,
    mean_stop_reject = mean_rej,
    sd_stop_reject = sd_rej,
    se_stop_reject = se_rej,
    stop_reject_ci_low = if (is.finite(mean_rej) && is.finite(se_rej)) max(0, mean_rej - 1.96 * se_rej) else mean_rej,
    stop_reject_ci_high = if (is.finite(mean_rej) && is.finite(se_rej)) mean_rej + 1.96 * se_rej else mean_rej,
    mean_stop_all = mean_all,
    sd_stop_all = sd_all,
    se_stop_all = se_all,
    csv = path,
    row.names = NULL
  )
}

save_summary <- function(eps_grid, trials, outdir, M_burn) {
  rows <- lapply(eps_grid, function(eps) {
    path <- file.path(outdir, sprintf("stops_burn%d_%sNONDP_bidir.csv", M_burn, fmt_eps(eps)))
    summarize_one_csv(path, eps = eps, trials_expected = trials)
  })
  summary_df <- do.call(rbind, rows)
  summary_df <- summary_df[order(summary_df$eps), , drop = FALSE]
  out <- file.path(outdir, "nondp_gaussian_auditor_summary_bidir.csv")
  utils::write.csv(summary_df, out, row.names = FALSE)
  cat("Saved summary:", out, "\n")
  summary_df
}

plot_with_errorbars <- function(x, y, lo, hi, xlab, ylab, main, outfile, ylim = NULL) {
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) {
    warning("No finite data for plot: ", outfile)
    return(invisible(FALSE))
  }
  x <- x[ok]; y <- y[ok]; lo <- lo[ok]; hi <- hi[ok]
  if (is.null(ylim)) {
    ymin <- min(c(y, lo), na.rm = TRUE)
    ymax <- max(c(y, hi), na.rm = TRUE)
    pad <- 0.06 * max(1e-9, ymax - ymin)
    ylim <- c(max(0, ymin - pad), ymax + pad)
  }

  grDevices::png(outfile, width = 1200, height = 850, res = 150)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(oldpar); grDevices::dev.off() }, add = TRUE)
  graphics::par(mar = c(5.2, 5.4, 3.2, 1.5), cex.axis = 1.1, cex.lab = 1.25)
  graphics::plot(
    x, y,
    type = "b", pch = 19, lwd = 2.5,
    xlab = xlab, ylab = ylab, main = main,
    ylim = ylim
  )
  graphics::grid(col = "gray85", lty = "dotted")
  graphics::lines(x, y, type = "b", pch = 19, lwd = 2.5, col = "steelblue")
  graphics::arrows(
    x0 = x, y0 = lo,
    x1 = x, y1 = hi,
    angle = 90, code = 3, length = 0.05,
    col = "gray25", lwd = 1.2
  )
  invisible(TRUE)
}

make_figures <- function(summary_df, outdir) {
  stop_png <- file.path(outdir, "NonDPGaussian1_bidir_avg_stop_to_reject_by_eps.png")
  rej_png  <- file.path(outdir, "NonDPGaussian1_bidir_rejection_rate_by_eps.png")

  plot_with_errorbars(
    x = summary_df$eps,
    y = summary_df$mean_stop_reject,
    lo = summary_df$stop_reject_ci_low,
    hi = summary_df$stop_reject_ci_high,
    xlab = expression(Claimed~epsilon),
    ylab = "Average stopping time to rejection",
    main = "NonDPGaussian1: sample efficiency, both orientations",
    outfile = stop_png
  )

  plot_with_errorbars(
    x = summary_df$eps,
    y = summary_df$rejection_rate,
    lo = summary_df$rejection_ci_low,
    hi = summary_df$rejection_ci_high,
    xlab = expression(Claimed~epsilon),
    ylab = "Empirical rejection rate",
    main = "NonDPGaussian1: empirical rejection rate, both orientations",
    outfile = rej_png,
    ylim = c(0, 1)
  )

  cat("Saved figure:", stop_png, "\n")
  cat("Saved figure:", rej_png, "\n")
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

src_dir_arg <- get_flag_value(args, "--src-dir", NULL)
src_dir <- find_source_dir(src_dir_arg)
source_audit_files(src_dir)

outdir <- get_flag_value(args, "--outdir", "results/nondp_gaussian_auditor_sweep")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

eps_grid <- parse_eps_grid(get_flag_value(args, "--eps-grid", NULL))
trials <- as.integer(get_flag_value(args, "--trials", "1000"))
M_burn <- as.integer(get_flag_value(args, "--burn", "50"))
mc_cores <- as.integer(get_flag_value(args, "--cores", "8"))
chunk_size <- as.integer(get_flag_value(args, "--chunk-size", "0"))
if (is.na(chunk_size) || chunk_size <= 0L) chunk_size <- max(25L, 2L * mc_cores)

h <- as.numeric(get_flag_value(args, "--h", "0.1"))
alpha <- as.numeric(get_flag_value(args, "--alpha", "0.05"))
delta_claim <- as.numeric(get_flag_value(args, "--delta", "1e-5"))
eta_search_max <- as.numeric(get_flag_value(args, "--eta-search-max", "15"))
max_iter <- as.integer(get_flag_value(args, "--max-iter", "10000"))
eval_step <- as.integer(get_flag_value(args, "--eval-step", "10"))
refit_threshold <- as.numeric(get_flag_value(args, "--refit-threshold", "0.10"))
qalpha_sims <- as.integer(get_flag_value(args, "--qalpha-sims", "10000"))
qalpha_kmax <- as.integer(get_flag_value(args, "--qalpha-kmax", "10000"))
seed <- as.integer(get_flag_value(args, "--seed", "20260524"))
if (is.na(seed)) seed <- NA_integer_

show_progress <- get_flag_bool(args, "--show-progress", TRUE)
overwrite <- get_flag_bool(args, "--overwrite", FALSE)
plot_only <- get_flag_bool(args, "--plot-only", FALSE)
boundary_tweaks <- get_flag_bool(args, "--boundary-tweaks", FALSE)
split_alpha_orientations <- get_flag_bool(args, "--split-alpha-orientations", TRUE)

if (is.na(trials) || trials < 1L) stop("--trials must be a positive integer")
if (is.na(M_burn) || M_burn < 1L) stop("--burn must be a positive integer")
if (is.na(mc_cores) || mc_cores < 1L) stop("--cores must be a positive integer")

avail <- parallel::detectCores(logical = TRUE)
mc_cores <- min(mc_cores, avail)

cat("Source dir:", src_dir, "\n")
cat("Output dir:", outdir, "\n")
cat("Epsilon grid:", paste(fmt_eps(eps_grid), collapse = ", "), "\n")
cat("Trials:", trials, "burn:", M_burn, "cores:", mc_cores, "chunk-size:", chunk_size, "\n")
cat("Claim delta:", delta_claim, "\n")
cat("Both orientations: TRUE; split alpha over orientations:", split_alpha_orientations, "\n")

x1 <- c(0)
x2 <- c(0, 1)
classifier <- make_kde_classifier()

if (!isTRUE(plot_only)) {
  qalpha_file <- file.path(
    outdir,
    sprintf(
      "qalpha_burn%d_alpha%s_dirs%s_sims%d_kmax%d_eval%d.txt",
      M_burn,
      formatC(alpha, format = "fg", digits = 8),
      ifelse(split_alpha_orientations, "2split", "2nosplit"),
      qalpha_sims,
      qalpha_kmax,
      eval_step
    )
  )

  qalpha_tail <- if (isTRUE(split_alpha_orientations)) alpha / 4 else alpha / 2
  cat("q_alpha Gaussian-tail probability:", qalpha_tail, "\n")

  if (file.exists(qalpha_file)) {
    q_alpha <- as.numeric(readLines(qalpha_file, warn = FALSE)[1L])
    cat("Loaded cached q_alpha:", q_alpha, "from", qalpha_file, "\n")
  } else {
    if (!is.na(seed)) set.seed(seed)
    cat("Computing q_alpha ...\n")
    q_alpha <- simulate_gaussian_sup_quantile(
      M = M_burn,
      alpha = qalpha_tail,
      sims = qalpha_sims,
      k_max = qalpha_kmax,
      eval_step = eval_step
    )
    writeLines(as.character(q_alpha), qalpha_file)
    cat("Saved q_alpha:", q_alpha, "to", qalpha_file, "\n")
  }

  for (ii in seq_along(eps_grid)) {
    eps <- eps_grid[ii]
    outfile <- file.path(outdir, sprintf("stops_burn%d_%sNONDP_bidir.csv", M_burn, fmt_eps(eps)))

    if (isTRUE(overwrite) && file.exists(outfile)) {
      file.remove(outfile)
    }

    existing <- read_existing_csv(outfile)
    done <- get_completed_trials(existing, trials)
    missing <- setdiff(seq_len(trials), done)

    cat("\n--- NonDPGaussian1 eps=", fmt_eps(eps), " ---\n", sep = "")
    cat("Completed:", length(done), "/", trials, " Missing:", length(missing), "\n")

    if (length(missing) > 0L) {
      run_missing_trials_for_eps(
        eps = eps,
        eps_index = ii,
        missing_trials = missing,
        existing_df = existing,
        outfile = outfile,
        q_alpha = q_alpha,
        classifier = classifier,
        x1 = x1,
        x2 = x2,
        M_burn = M_burn,
        h = h,
        eta_search_max = eta_search_max,
        max_iter = max_iter,
        eval_step = eval_step,
        refit_threshold = refit_threshold,
        boundary_tweaks = boundary_tweaks,
        mc_cores = mc_cores,
        show_progress = show_progress,
        chunk_size = chunk_size,
        seed = seed,
        delta_claim = delta_claim
      )
    } else {
      cat("Resume: nothing to run for eps=", fmt_eps(eps), "\n", sep = "")
    }
  }
}

summary_df <- save_summary(eps_grid, trials, outdir, M_burn)
make_figures(summary_df, outdir)

cat("Done.\n")
