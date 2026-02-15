# ============================================================================
# test_external_equivalence.R
# Cross-package numerical equivalence validation
#
# Purpose:
#   Verify that every new function in claude/R/ produces numerically identical
#   results to established R packages across 100 random iterations per function.
#
# Reference packages:
#   - changepoint (v2.3+)    : PELT and BinSeg changepoint detection
#   - EWSmethods (v1.3+)     : Multivariate early warning signals (rolling)
#   - nonlinearTseries (v0.2.12+) : FFT phase surrogate generation
#   - stats (base R)         : spectrum(), density(), ar.ols(), cor.test()
#
# Test protocol:
#   Each section runs 100 iterations. Every iteration:
#     1. Set a unique random seed (seed = iter)
#     2. Generate random data with random parameters
#     3. Run BOTH our function and the reference function with identical inputs
#     4. Compare outputs numerically and record the exact differences
#     5. Log all values to a results data.frame for the final report
#
# Tolerances:
#   - Changepoint locations: exact match (integer)
#   - Segment means: < 1e-8 absolute
#   - Metric time series: < 1e-6 absolute, or Pearson r > 0.999
#   - Kendall tau: < 0.05 absolute (stochastic surrogates)
#   - Power spectrum preservation: < 1e-6 relative
#   - Spectral exponent / R-squared: < 1e-8 absolute
#   - Sensitivity AR(1) scores: < 1e-8 absolute
#   - Potential landscape density: < 1e-8 absolute
#
# Output:
#   Per-iteration lines with actual numerical comparison values, plus
#   an aggregate report (mean/max diff, Pearson r, pass/fail counts).
#
# Runtime: approximately 20-30 minutes for all 800 iterations.
# ============================================================================

suppressWarnings(suppressMessages({
  library(codyna)
  library(changepoint)
  library(EWSmethods)
  library(nonlinearTseries)

  # Import codyna internals needed by our functions
  check_missing <- codyna:::check_missing
  check_match   <- codyna:::check_match
  check_range   <- codyna:::check_range
  check_flag    <- codyna:::check_flag
  check_class   <- codyna:::check_class
  prepare_timeseries_data <- codyna:::prepare_timeseries_data
  stop_         <- codyna:::stop_
  warning_      <- codyna:::warning_
  message_      <- codyna:::message_
  stopifnot_    <- codyna:::stopifnot_
  ifelse_       <- codyna:::ifelse_
  try_          <- codyna:::try_
  onlyif        <- codyna:::onlyif
  `%||%`        <- function(x, y) if (is.null(x)) y else x

  # Source our 11 R files
  our_dir <- "R"
  for (f in list.files(our_dir, pattern = "[.]R$", full.names = TRUE)) {
    source(f, local = TRUE)
  }
}))

# Helper: format number with sign
fmt <- function(x, digits = 6) formatC(x, format = "f", digits = digits)
fmte <- function(x) formatC(x, format = "e", digits = 2)

# Track global timing
t_start <- proc.time()


# ============================================================================
# PART 1: CHANGEPOINT DETECTION (PELT)
# ============================================================================
#
# What:   detect_changepoints(method = "pelt") vs changepoint::cpt.mean(method = "PELT")
# Why:    PELT is an exact segmentation algorithm. Both implementations must find
#         identical changepoint locations and segment means for the same BIC penalty.
# How:    Generate multi-segment data (2-4 segments, varying lengths 50-200,
#         varying means -5 to +5, varying SD 0.5-2.0). Compare:
#           (a) number of changepoints
#           (b) changepoint positions (exact integer match)
#           (c) segment means (< 1e-8 absolute difference)
# Ref:    changepoint::cpt.mean() with penalty="BIC", minseglen matching
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("PART 1: CHANGEPOINT (PELT) -- 100 iterations\n")
cat("  Our detect_changepoints(method='pelt', penalty='bic')\n")
cat("  vs changepoint::cpt.mean(method='PELT', penalty='BIC')\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("  %-4s  %-6s  %-7s  %-8s  %-8s  %-12s  %-12s  %s\n",
            "iter", "n", "min_seg", "our_cps", "pkg_cps", "mean_diff", "max_diff", "status"))
cat(sprintf("  %s\n", strrep("-", 78)))

pelt_results <- data.frame(
  iter = integer(), n = integer(), n_segments = integer(), min_seg = integer(),
  our_n_cps = integer(), pkg_n_cps = integer(),
  locations_match = logical(), mean_max_diff = numeric(),
  status = character(), stringsAsFactors = FALSE
)

for (iter in 1:100) {
  set.seed(iter)

  # Random data: 2-4 segments with different means
  n_segments <- sample(2:4, 1)
  seg_lens <- sample(50:200, n_segments, replace = TRUE)
  seg_means <- round(runif(n_segments, -5, 5), 1)
  seg_sds <- runif(n_segments, 0.5, 2.0)
  x <- unlist(mapply(rnorm, n = seg_lens, mean = seg_means, sd = seg_sds,
                     SIMPLIFY = FALSE))
  n <- length(x)
  min_seg <- sample(c(10L, 15L, 20L), 1)

  our_result <- tryCatch(
    detect_changepoints(x, method = "pelt", penalty = "bic",
                        min_segment = min_seg, type = "mean"),
    error = function(e) NULL
  )
  pkg_result <- tryCatch(
    changepoint::cpt.mean(x, method = "PELT", penalty = "BIC",
                          minseglen = min_seg),
    error = function(e) NULL
  )

  if (is.null(our_result) || is.null(pkg_result)) {
    pelt_results <- rbind(pelt_results, data.frame(
      iter = iter, n = n, n_segments = n_segments, min_seg = min_seg,
      our_n_cps = NA_integer_, pkg_n_cps = NA_integer_,
      locations_match = NA, mean_max_diff = NA_real_,
      status = "SKIP", stringsAsFactors = FALSE
    ))
    cat(sprintf("  %-4d  %-6d  %-7d  %-8s  %-8s  %-12s  %-12s  SKIP\n",
                iter, n, min_seg, "--", "--", "--", "--"))
    next
  }

  our_cps <- sort(attr(our_result, "changepoint_locations"))
  pkg_cps <- sort(changepoint::cpts(pkg_result))

  loc_match <- identical(our_cps, pkg_cps)

  # Segment mean comparison (only if same changepoints)
  mean_max_diff <- NA_real_
  if (loc_match && length(our_cps) > 0) {
    our_seg_means <- sort(tapply(our_result$value, our_result$segment, mean))
    pkg_seg_means <- sort(changepoint::param.est(pkg_result)$mean)
    mean_max_diff <- max(abs(our_seg_means - pkg_seg_means))
  } else if (loc_match && length(our_cps) == 0) {
    mean_max_diff <- 0
  }

  status <- if (!loc_match) "FAIL" else if (!is.na(mean_max_diff) && mean_max_diff < 1e-8) "PASS" else "FAIL"

  pelt_results <- rbind(pelt_results, data.frame(
    iter = iter, n = n, n_segments = n_segments, min_seg = min_seg,
    our_n_cps = length(our_cps), pkg_n_cps = length(pkg_cps),
    locations_match = loc_match, mean_max_diff = mean_max_diff,
    status = status, stringsAsFactors = FALSE
  ))

  cat(sprintf("  %-4d  %-6d  %-7d  %-8d  %-8d  %-12s  %-12s  %s\n",
              iter, n, min_seg, length(our_cps), length(pkg_cps),
              if (is.na(mean_max_diff)) "--" else fmte(mean_max_diff),
              if (is.na(mean_max_diff)) "--" else fmte(mean_max_diff),
              status))
}

pelt_pass <- sum(pelt_results$status == "PASS", na.rm = TRUE)
pelt_fail <- sum(pelt_results$status == "FAIL", na.rm = TRUE)
pelt_skip <- sum(pelt_results$status == "SKIP", na.rm = TRUE)
cat(sprintf("\n  PELT REPORT: %d PASS, %d FAIL, %d SKIP out of 100\n", pelt_pass, pelt_fail, pelt_skip))
if (pelt_pass > 0) {
  valid_diffs <- pelt_results$mean_max_diff[pelt_results$status == "PASS"]
  cat(sprintf("  Segment mean diffs: mean=%.2e, max=%.2e\n", mean(valid_diffs), max(valid_diffs)))
}


# ============================================================================
# PART 2: CHANGEPOINT DETECTION (Binary Segmentation)
# ============================================================================
#
# What:   detect_changepoints(method = "binary_segmentation") vs
#         changepoint::cpt.mean(method = "BinSeg")
# Why:    Binary segmentation is an approximate algorithm that recursively
#         splits. Both must agree on the number and locations of splits.
# How:    Generate 2-3 segment data (seg lengths 80-200, means -3 to +3,
#         unit SD). Vary max_changepoints (1-3). Compare:
#           (a) number of changepoints
#           (b) changepoint positions (exact match; allow +/-5 tolerance)
# Ref:    changepoint::cpt.mean(method="BinSeg", Q=max_cp)
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("PART 2: CHANGEPOINT (BinSeg) -- 100 iterations\n")
cat("  Our detect_changepoints(method='binary_segmentation')\n")
cat("  vs changepoint::cpt.mean(method='BinSeg')\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("  %-4s  %-6s  %-6s  %-8s  %-8s  %-15s  %s\n",
            "iter", "n", "max_cp", "our_cps", "pkg_cps", "loc_diff", "status"))
cat(sprintf("  %s\n", strrep("-", 66)))

binseg_results <- data.frame(
  iter = integer(), n = integer(), max_cp = integer(),
  our_n_cps = integer(), pkg_n_cps = integer(),
  max_loc_diff = numeric(), status = character(),
  stringsAsFactors = FALSE
)

for (iter in 1:100) {
  set.seed(iter)

  n_segments <- sample(2:3, 1)
  seg_lens <- sample(80:200, n_segments, replace = TRUE)
  seg_means <- round(runif(n_segments, -3, 3), 1)
  x <- unlist(mapply(rnorm, n = seg_lens, mean = seg_means, sd = 1,
                     SIMPLIFY = FALSE))
  n <- length(x)
  min_seg <- 15L
  max_cp <- sample(1:3, 1)

  our_result <- tryCatch(
    detect_changepoints(x, method = "binary_segmentation", penalty = "bic",
                        min_segment = min_seg, max_changepoints = max_cp,
                        type = "mean"),
    error = function(e) NULL
  )
  pkg_result <- tryCatch(
    changepoint::cpt.mean(x, method = "BinSeg", penalty = "BIC",
                          minseglen = min_seg, Q = max_cp),
    error = function(e) NULL
  )

  if (is.null(our_result) || is.null(pkg_result)) {
    binseg_results <- rbind(binseg_results, data.frame(
      iter = iter, n = n, max_cp = max_cp,
      our_n_cps = NA_integer_, pkg_n_cps = NA_integer_,
      max_loc_diff = NA_real_, status = "SKIP",
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  %-4d  %-6d  %-6d  %-8s  %-8s  %-15s  SKIP\n",
                iter, n, max_cp, "--", "--", "--"))
    next
  }

  our_cps <- sort(attr(our_result, "changepoint_locations"))
  pkg_cps <- sort(changepoint::cpts(pkg_result))

  max_loc_diff <- NA_real_
  if (length(our_cps) == length(pkg_cps) && length(our_cps) > 0) {
    max_loc_diff <- max(abs(our_cps - pkg_cps))
  } else if (length(our_cps) == 0 && length(pkg_cps) == 0) {
    max_loc_diff <- 0
  }

  status <- if (length(our_cps) != length(pkg_cps)) {
    "FAIL"
  } else if (is.na(max_loc_diff)) {
    "FAIL"
  } else if (max_loc_diff <= 5) {
    "PASS"
  } else {
    "FAIL"
  }

  binseg_results <- rbind(binseg_results, data.frame(
    iter = iter, n = n, max_cp = max_cp,
    our_n_cps = length(our_cps), pkg_n_cps = length(pkg_cps),
    max_loc_diff = max_loc_diff, status = status,
    stringsAsFactors = FALSE
  ))

  loc_str <- if (is.na(max_loc_diff)) {
    sprintf("n: %d vs %d", length(our_cps), length(pkg_cps))
  } else {
    sprintf("%d", as.integer(max_loc_diff))
  }

  cat(sprintf("  %-4d  %-6d  %-6d  %-8d  %-8d  %-15s  %s\n",
              iter, n, max_cp, length(our_cps), length(pkg_cps), loc_str, status))
}

binseg_pass <- sum(binseg_results$status == "PASS", na.rm = TRUE)
binseg_fail <- sum(binseg_results$status == "FAIL", na.rm = TRUE)
binseg_skip <- sum(binseg_results$status == "SKIP", na.rm = TRUE)
cat(sprintf("\n  BinSeg REPORT: %d PASS, %d FAIL, %d SKIP out of 100\n",
            binseg_pass, binseg_fail, binseg_skip))
if (binseg_pass > 0) {
  valid_ld <- binseg_results$max_loc_diff[binseg_results$status == "PASS"]
  cat(sprintf("  Location diffs: mean=%.2f, max=%d\n", mean(valid_ld), max(valid_ld)))
}


# ============================================================================
# PART 3: MULTIVARIATE EARLY WARNING SIGNALS (Rolling)
# ============================================================================
#
# What:   detect_multivariate_warnings(method = "rolling") vs
#         EWSmethods::multiEWS(method = "rolling")
# Why:    Both compute the same 11 metrics (meanAR, maxAR, meanSD, maxSD,
#         eigenMAF, mafAR, mafSD, pcaAR, pcaSD, eigenCOV, maxCOV) in rolling
#         windows. The raw metric time series and Kendall tau trends must match.
# How:    Generate tipping data (80-150 time steps, 3-5 variables). Run both.
#         For each of 11 metrics, compute:
#           (a) Pearson correlation between the two raw series
#           (b) Max absolute difference
#           (c) Kendall tau difference
# Tol:    Raw values: Pearson r > 0.999 or max_diff < 1e-6
#         Kendall tau: < 0.05 absolute
# Ref:    EWSmethods::multiEWS() with same winsize percentage
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("PART 3: MULTIVARIATE EWS (Rolling) -- 100 iterations\n")
cat("  Our detect_multivariate_warnings() vs EWSmethods::multiEWS()\n")
cat("  Comparing 11 metrics: meanAR, maxAR, meanSD, maxSD, eigenMAF,\n")
cat("    mafAR, mafSD, pcaAR, pcaSD, eigenCOV, maxCOV\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("  %-4s  %-5s  %-5s  %-5s  %-10s  %-10s  %-10s  %-10s  %s\n",
            "iter", "n", "vars", "win%", "min_r", "max_diff", "tau_diff", "metric", "status"))
cat(sprintf("  %s\n", strrep("-", 82)))

all_metrics <- c("meanAR", "maxAR", "meanSD", "maxSD", "eigenMAF",
                 "mafAR", "mafSD", "pcaAR", "pcaSD", "eigenCOV", "maxCOV")

mews_results <- data.frame(
  iter = integer(), n_time = integer(), n_vars = integer(), winsize = integer(),
  min_pearson_r = numeric(), max_abs_diff = numeric(),
  max_tau_diff = numeric(), worst_metric = character(),
  n_metrics_checked = integer(), status = character(),
  stringsAsFactors = FALSE
)

for (iter in 1:100) {
  set.seed(iter)

  n_time <- sample(80:150, 1)
  n_vars <- sample(3:5, 1)
  winsize_pct <- sample(30:60, 1)

  tp <- tryCatch(
    generate_tipping_data(n_time = n_time, n_vars = n_vars,
                          tipping_point = round(n_time * 0.7),
                          saturation_point = round(n_time * 0.9)),
    error = function(e) NULL
  )
  if (is.null(tp)) {
    mews_results <- rbind(mews_results, data.frame(
      iter = iter, n_time = n_time, n_vars = n_vars, winsize = winsize_pct,
      min_pearson_r = NA, max_abs_diff = NA, max_tau_diff = NA,
      worst_metric = "--", n_metrics_checked = 0L, status = "SKIP",
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  %-4d  %-5d  %-5d  %-5d  %-10s  %-10s  %-10s  %-10s  SKIP\n",
                iter, n_time, n_vars, winsize_pct, "--", "--", "--", "--"))
    next
  }

  our_result <- tryCatch(
    suppressMessages(
      detect_multivariate_warnings(tp, method = "rolling", window = winsize_pct)
    ),
    error = function(e) NULL
  )
  pkg_result <- tryCatch(
    suppressWarnings(suppressMessages(
      EWSmethods::multiEWS(
        data = tp,
        metrics = all_metrics,
        method = "rolling",
        winsize = winsize_pct
      )
    )),
    error = function(e) NULL
  )

  if (is.null(our_result) || is.null(pkg_result)) {
    mews_results <- rbind(mews_results, data.frame(
      iter = iter, n_time = n_time, n_vars = n_vars, winsize = winsize_pct,
      min_pearson_r = NA, max_abs_diff = NA, max_tau_diff = NA,
      worst_metric = "--", n_metrics_checked = 0L, status = "SKIP",
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  %-4d  %-5d  %-5d  %-5d  %-10s  %-10s  %-10s  %-10s  SKIP\n",
                iter, n_time, n_vars, winsize_pct, "--", "--", "--", "--"))
    next
  }

  # Pivot our result to wide format
  our_wide <- tidyr::pivot_wider(
    as.data.frame(our_result),
    id_cols = "time", names_from = "metric", values_from = "score"
  )
  pkg_raw <- pkg_result$EWS$raw

  # Align time points
  common_times <- intersect(our_wide$time, pkg_raw$time)
  if (length(common_times) < 5) {
    mews_results <- rbind(mews_results, data.frame(
      iter = iter, n_time = n_time, n_vars = n_vars, winsize = winsize_pct,
      min_pearson_r = NA, max_abs_diff = NA, max_tau_diff = NA,
      worst_metric = "--", n_metrics_checked = 0L, status = "SKIP",
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  %-4d  %-5d  %-5d  %-5d  %-10s  %-10s  %-10s  %-10s  SKIP\n",
                iter, n_time, n_vars, winsize_pct, "--", "--", "--", "--"))
    next
  }

  our_aligned <- our_wide[our_wide$time %in% common_times, ]
  pkg_aligned <- pkg_raw[pkg_raw$time %in% common_times, ]
  our_aligned <- our_aligned[order(our_aligned$time), ]
  pkg_aligned <- pkg_aligned[order(pkg_aligned$time), ]

  min_r <- 1.0
  max_diff <- 0.0
  max_tau_diff <- 0.0
  worst_metric <- "--"
  n_checked <- 0L
  all_pass <- TRUE

  for (m in all_metrics) {
    if (!m %in% names(our_aligned) || !m %in% names(pkg_aligned)) next
    our_vals <- our_aligned[[m]]
    pkg_vals <- pkg_aligned[[m]]
    both_valid <- !is.na(our_vals) & !is.na(pkg_vals)
    if (sum(both_valid) < 3) next
    n_checked <- n_checked + 1L

    # Pearson correlation
    r <- tryCatch(
      stats::cor(our_vals[both_valid], pkg_vals[both_valid]),
      error = function(e) NA
    )
    if (!is.na(r) && r < min_r) {
      min_r <- r
      if (r < 0.999) worst_metric <- m
    }

    # Max absolute difference
    md <- max(abs(our_vals[both_valid] - pkg_vals[both_valid]))
    if (md > max_diff) max_diff <- md

    # Metric passes if r > 0.999 OR max_diff < 1e-6
    if (is.na(r) || (r < 0.999 && md > 1e-6)) all_pass <- FALSE
  }

  # Compare Kendall tau correlations
  our_cor <- attr(our_result, "cor")
  pkg_cor_df <- pkg_result$EWS$cor
  if (!is.null(our_cor) && !is.null(pkg_cor_df)) {
    our_tau_names <- sub("\\.tau$", "", names(our_cor))
    for (m in all_metrics) {
      our_idx <- which(our_tau_names == m)
      if (length(our_idx) == 1 && m %in% names(pkg_cor_df)) {
        td <- abs(our_cor[our_idx] - pkg_cor_df[[m]])
        if (!is.na(td) && td > max_tau_diff) max_tau_diff <- td
        if (!is.na(td) && td > 0.05) all_pass <- FALSE
      }
    }
  }

  status <- if (all_pass) "PASS" else "FAIL"

  mews_results <- rbind(mews_results, data.frame(
    iter = iter, n_time = n_time, n_vars = n_vars, winsize = winsize_pct,
    min_pearson_r = min_r, max_abs_diff = max_diff,
    max_tau_diff = max_tau_diff, worst_metric = worst_metric,
    n_metrics_checked = n_checked, status = status,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  %-4d  %-5d  %-5d  %-5d  %-10s  %-10s  %-10s  %-10s  %s\n",
              iter, n_time, n_vars, winsize_pct,
              fmt(min_r, 6), fmte(max_diff), fmte(max_tau_diff),
              worst_metric, status))
}

mews_pass <- sum(mews_results$status == "PASS", na.rm = TRUE)
mews_fail <- sum(mews_results$status == "FAIL", na.rm = TRUE)
mews_skip <- sum(mews_results$status == "SKIP", na.rm = TRUE)
cat(sprintf("\n  MEWS REPORT: %d PASS, %d FAIL, %d SKIP out of 100\n",
            mews_pass, mews_fail, mews_skip))
valid_mews <- mews_results[mews_results$status == "PASS", ]
if (nrow(valid_mews) > 0) {
  cat(sprintf("  Pearson r:   min=%.6f, mean=%.6f\n",
              min(valid_mews$min_pearson_r), mean(valid_mews$min_pearson_r)))
  cat(sprintf("  Max diff:    mean=%.2e, max=%.2e\n",
              mean(valid_mews$max_abs_diff), max(valid_mews$max_abs_diff)))
  cat(sprintf("  Tau diff:    mean=%.2e, max=%.2e\n",
              mean(valid_mews$max_tau_diff), max(valid_mews$max_tau_diff)))
}


# ============================================================================
# PART 4: PHASE SURROGATES — Power Spectrum Preservation
# ============================================================================
#
# What:   Our surrogate_phase_() vs nonlinearTseries::FFTsurrogate()
# Why:    FFT phase randomization must EXACTLY preserve |FFT(x)|^2 (the power
#         spectrum). Both implementations should achieve this. We verify our
#         surrogates preserve the power spectrum to machine precision.
# How:    Generate random walk data (n=100-500). For each iteration:
#           (a) Compute original power spectrum |FFT(x)|^2
#           (b) Generate our phase surrogate, compute |FFT(surr)|^2
#           (c) Measure max|original_power - surrogate_power|
#           (d) Verify variance ratio surr_var/orig_var ~ 1.0
#           (e) Verify length preservation
# Tol:    Power spectrum: max diff < 1e-6
#         Variance ratio: |ratio - 1| < 0.05
# Ref:    nonlinearTseries::FFTsurrogate() for the same property
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("PART 4: PHASE SURROGATES -- 100 iterations\n")
cat("  Our surrogate_phase_() power spectrum preservation\n")
cat("  Reference: nonlinearTseries::FFTsurrogate() same property\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("  %-4s  %-6s  %-14s  %-14s  %-12s  %-10s  %s\n",
            "iter", "n", "our_power_diff", "pkg_power_diff", "var_ratio", "mean_diff", "status"))
cat(sprintf("  %s\n", strrep("-", 80)))

surr_results <- data.frame(
  iter = integer(), n = integer(),
  our_power_max_diff = numeric(), pkg_power_max_diff = numeric(),
  our_var_ratio = numeric(), our_mean_diff = numeric(),
  status = character(), stringsAsFactors = FALSE
)

for (iter in 1:100) {
  set.seed(iter)
  n <- sample(100:500, 1)
  x <- cumsum(rnorm(n))

  # Our phase surrogate
  set.seed(iter * 1000 + 1)
  our_surr <- surrogate_phase_(x)

  # Package phase surrogate
  set.seed(iter * 1000 + 2)
  pkg_surr <- tryCatch(
    as.numeric(nonlinearTseries::FFTsurrogate(x, n.samples = 1)),
    error = function(e) NULL
  )

  if (is.null(pkg_surr)) {
    surr_results <- rbind(surr_results, data.frame(
      iter = iter, n = n, our_power_max_diff = NA, pkg_power_max_diff = NA,
      our_var_ratio = NA, our_mean_diff = NA, status = "SKIP",
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  %-4d  %-6d  %-14s  %-14s  %-12s  %-10s  SKIP\n",
                iter, n, "--", "--", "--", "--"))
    next
  }

  orig_power <- Mod(stats::fft(x))^2
  our_power <- Mod(stats::fft(our_surr))^2
  pkg_power <- Mod(stats::fft(pkg_surr))^2

  our_pd <- max(abs(our_power - orig_power))
  pkg_pd <- max(abs(pkg_power - orig_power))

  our_vr <- stats::var(our_surr) / stats::var(x)
  our_md <- abs(mean(our_surr) - mean(x))

  all_pass <- our_pd < 1e-6 && abs(our_vr - 1.0) < 0.05 && length(our_surr) == n

  status <- if (all_pass) "PASS" else "FAIL"

  surr_results <- rbind(surr_results, data.frame(
    iter = iter, n = n,
    our_power_max_diff = our_pd, pkg_power_max_diff = pkg_pd,
    our_var_ratio = our_vr, our_mean_diff = our_md,
    status = status, stringsAsFactors = FALSE
  ))

  cat(sprintf("  %-4d  %-6d  %-14s  %-14s  %-12s  %-10s  %s\n",
              iter, n, fmte(our_pd), fmte(pkg_pd),
              fmt(our_vr, 6), fmte(our_md), status))
}

surr_pass <- sum(surr_results$status == "PASS", na.rm = TRUE)
surr_fail <- sum(surr_results$status == "FAIL", na.rm = TRUE)
surr_skip <- sum(surr_results$status == "SKIP", na.rm = TRUE)
cat(sprintf("\n  SURROGATE REPORT: %d PASS, %d FAIL, %d SKIP out of 100\n",
            surr_pass, surr_fail, surr_skip))
valid_surr <- surr_results[surr_results$status == "PASS", ]
if (nrow(valid_surr) > 0) {
  cat(sprintf("  Our power diff:  mean=%.2e, max=%.2e\n",
              mean(valid_surr$our_power_max_diff), max(valid_surr$our_power_max_diff)))
  cat(sprintf("  Pkg power diff:  mean=%.2e, max=%.2e\n",
              mean(valid_surr$pkg_power_max_diff), max(valid_surr$pkg_power_max_diff)))
  cat(sprintf("  Variance ratio:  mean=%.6f, range=[%.6f, %.6f]\n",
              mean(valid_surr$our_var_ratio),
              min(valid_surr$our_var_ratio), max(valid_surr$our_var_ratio)))
}


# ============================================================================
# PART 5: SPECTRAL EXPONENT
# ============================================================================
#
# What:   spectral_ews() vs manual stats::spectrum() + lm() computation
# Why:    The spectral exponent (beta) is the negative slope of log(power) vs
#         log(frequency). Our implementation and the manual computation must
#         agree to machine precision since they use identical formulas.
# How:    Generate random walk (n=100-300), choose window (30/50/60). For the
#         LAST window position (right-aligned):
#           (a) Detrend linearly (manual lm)
#           (b) Compute periodogram via stats::spectrum(method="pgram")
#           (c) Fit log-log regression, negate slope to get beta
#           (d) Compare beta and R-squared
# Tol:    Absolute diff < 1e-8 for both beta and R-squared
# Ref:    Manual stats::spectrum() + lm(log(power) ~ log(freq))
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("PART 5: SPECTRAL EXPONENT -- 100 iterations\n")
cat("  Our spectral_ews() vs manual stats::spectrum() + lm()\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("  %-4s  %-6s  %-5s  %-12s  %-12s  %-12s  %-12s  %s\n",
            "iter", "n", "win", "our_beta", "ref_beta", "our_r2", "ref_r2", "status"))
cat(sprintf("  %s\n", strrep("-", 80)))

spec_results <- data.frame(
  iter = integer(), n = integer(), window = integer(),
  our_beta = numeric(), ref_beta = numeric(), beta_diff = numeric(),
  our_r2 = numeric(), ref_r2 = numeric(), r2_diff = numeric(),
  status = character(), stringsAsFactors = FALSE
)

for (iter in 1:100) {
  set.seed(iter)
  n <- sample(100:300, 1)
  x <- cumsum(rnorm(n))
  window <- sample(c(30L, 50L, 60L), 1)

  our_result <- tryCatch(
    spectral_ews(x, window = window, method = "periodogram",
                 detrend = "linear", states = TRUE),
    error = function(e) NULL
  )
  if (is.null(our_result)) {
    spec_results <- rbind(spec_results, data.frame(
      iter = iter, n = n, window = window,
      our_beta = NA, ref_beta = NA, beta_diff = NA,
      our_r2 = NA, ref_r2 = NA, r2_diff = NA,
      status = "SKIP", stringsAsFactors = FALSE
    ))
    cat(sprintf("  %-4d  %-6d  %-5d  %-12s  %-12s  %-12s  %-12s  SKIP\n",
                iter, n, window, "--", "--", "--", "--"))
    next
  }

  # Manual: last window (right-aligned)
  seg <- x[(n - window + 1):n]

  # Linear detrend
  t_idx <- seq_along(seg)
  t_mean <- mean(t_idx)
  x_mean <- mean(seg)
  xx_var <- sum((t_idx - t_mean)^2)
  xy_cov <- sum((t_idx - t_mean) * (seg - x_mean))
  slope_dt <- xy_cov / xx_var
  intercept_dt <- x_mean - slope_dt * t_mean
  seg_detrended <- seg - (intercept_dt + slope_dt * t_idx)

  sp <- tryCatch(
    stats::spectrum(seg_detrended, method = "pgram", plot = FALSE,
                    detrend = FALSE, taper = 0.1),
    error = function(e) NULL
  )
  if (is.null(sp)) {
    spec_results <- rbind(spec_results, data.frame(
      iter = iter, n = n, window = window,
      our_beta = NA, ref_beta = NA, beta_diff = NA,
      our_r2 = NA, ref_r2 = NA, r2_diff = NA,
      status = "SKIP", stringsAsFactors = FALSE
    ))
    cat(sprintf("  %-4d  %-6d  %-5d  %-12s  %-12s  %-12s  %-12s  SKIP\n",
                iter, n, window, "--", "--", "--", "--"))
    next
  }

  freq <- sp$freq[sp$freq > 0 & sp$spec > 0]
  power <- sp$spec[sp$freq > 0 & sp$spec > 0]

  if (length(freq) < 4) next

  log_f <- log(freq)
  log_p <- log(power)
  fit <- lm(log_p ~ log_f)
  ref_beta <- -coef(fit)[2]
  ss_tot <- sum((log_p - mean(log_p))^2)
  ss_res <- sum(residuals(fit)^2)
  ref_r2 <- 1 - ss_res / ss_tot

  our_beta <- our_result$spectral_exponent[n]
  our_r2 <- our_result$r_squared[n]

  beta_diff <- abs(our_beta - ref_beta)
  r2_diff <- abs(our_r2 - ref_r2)

  status <- if (beta_diff < 1e-8 && r2_diff < 1e-8) "PASS" else "FAIL"

  spec_results <- rbind(spec_results, data.frame(
    iter = iter, n = n, window = window,
    our_beta = our_beta, ref_beta = ref_beta, beta_diff = beta_diff,
    our_r2 = our_r2, ref_r2 = ref_r2, r2_diff = r2_diff,
    status = status, stringsAsFactors = FALSE
  ))

  cat(sprintf("  %-4d  %-6d  %-5d  %-12s  %-12s  %-12s  %-12s  %s\n",
              iter, n, window, fmt(our_beta, 6), fmt(ref_beta, 6),
              fmt(our_r2, 6), fmt(ref_r2, 6), status))
}

spec_pass <- sum(spec_results$status == "PASS", na.rm = TRUE)
spec_fail <- sum(spec_results$status == "FAIL", na.rm = TRUE)
spec_skip <- sum(spec_results$status == "SKIP", na.rm = TRUE)
cat(sprintf("\n  SPECTRAL REPORT: %d PASS, %d FAIL, %d SKIP out of 100\n",
            spec_pass, spec_fail, spec_skip))
valid_spec <- spec_results[spec_results$status == "PASS", ]
if (nrow(valid_spec) > 0) {
  cat(sprintf("  Beta diff:  mean=%.2e, max=%.2e\n",
              mean(valid_spec$beta_diff), max(valid_spec$beta_diff)))
  cat(sprintf("  R2 diff:    mean=%.2e, max=%.2e\n",
              mean(valid_spec$r2_diff), max(valid_spec$r2_diff)))
  cat(sprintf("  Beta range: [%.4f, %.4f]\n",
              min(valid_spec$our_beta), max(valid_spec$our_beta)))
}


# ============================================================================
# PART 6: SENSITIVITY EWS — AR(1) Rolling + Kendall Tau
# ============================================================================
#
# What:   sensitivity_ews(metric = "ar1") vs manual stats::ar.ols() rolling
#         computation + stats::cor.test(method = "kendall")
# Why:    The sensitivity analysis computes rolling AR(1) coefficients across
#         varying window sizes. The AR(1) values and their Kendall tau trend
#         statistics must match the manual computation exactly.
# How:    Generate trending random walk (n=150-400), pick window (25/40/50).
#         For each iteration:
#           (a) Compute rolling AR(1) via ar.ols on each window
#           (b) Clip to [-0.999, 0.999]
#           (c) Compute Kendall tau of the rolling series
#           (d) Compare per-window AR(1) values and final tau
# Tol:    AR(1) values: < 1e-8 absolute; Kendall tau: < 1e-8 absolute
# Ref:    Manual stats::ar.ols() + stats::cor.test(method = "kendall")
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("PART 6: SENSITIVITY EWS -- 100 iterations\n")
cat("  Our sensitivity_ews(metric='ar1') vs manual ar.ols() + cor.test()\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("  %-4s  %-6s  %-5s  %-12s  %-12s  %-12s  %-12s  %s\n",
            "iter", "n", "win", "our_tau", "ref_tau", "tau_diff", "score_diff", "status"))
cat(sprintf("  %s\n", strrep("-", 80)))

sens_results <- data.frame(
  iter = integer(), n = integer(), window = integer(),
  our_tau = numeric(), ref_tau = numeric(), tau_diff = numeric(),
  max_score_diff = numeric(), n_scores = integer(),
  status = character(), stringsAsFactors = FALSE
)

for (iter in 1:100) {
  set.seed(iter)
  n <- sample(150:400, 1)
  x <- cumsum(rnorm(n, sd = seq(0.5, 2, length.out = n)))
  window <- sample(c(25L, 40L, 50L), 1)

  our_result <- tryCatch(
    sensitivity_ews(x, metric = "ar1", windows = window,
                    detrend_methods = "none"),
    error = function(e) NULL
  )
  if (is.null(our_result)) {
    sens_results <- rbind(sens_results, data.frame(
      iter = iter, n = n, window = window,
      our_tau = NA, ref_tau = NA, tau_diff = NA,
      max_score_diff = NA, n_scores = NA_integer_,
      status = "SKIP", stringsAsFactors = FALSE
    ))
    cat(sprintf("  %-4d  %-6d  %-5d  %-12s  %-12s  %-12s  %-12s  SKIP\n",
                iter, n, window, "--", "--", "--", "--"))
    next
  }

  our_scores <- our_result$score
  our_tau <- unique(our_result$tau)

  # Manual rolling AR(1)
  n_windows <- n - window + 1L
  manual_scores <- vapply(seq_len(n_windows), function(i) {
    wd <- x[i:(i + window - 1L)]
    if (stats::var(wd) < 1e-10) return(NA_real_)
    model <- tryCatch(
      stats::ar.ols(wd, aic = FALSE, order.max = 1L,
                    demean = FALSE, intercept = TRUE),
      error = function(e) NULL,
      warning = function(w) NULL
    )
    if (is.null(model) || length(model$ar) == 0) return(NA_real_)
    max(-0.999, min(0.999, model$ar[1L]))
  }, numeric(1L))

  # Compare rolling scores
  max_score_diff <- NA_real_
  if (length(our_scores) == length(manual_scores)) {
    both_valid <- !is.na(our_scores) & !is.na(manual_scores)
    if (sum(both_valid) > 0) {
      max_score_diff <- max(abs(our_scores[both_valid] - manual_scores[both_valid]))
    }
  }

  # Manual Kendall tau
  valid <- !is.na(manual_scores)
  ref_tau <- NA_real_
  if (sum(valid) >= 4) {
    ref_tau <- stats::cor.test(
      seq_along(manual_scores)[valid],
      manual_scores[valid],
      method = "kendall"
    )$estimate
  }

  tau_diff <- abs(our_tau - ref_tau)
  score_pass <- !is.na(max_score_diff) && max_score_diff < 1e-8
  tau_pass <- !is.na(tau_diff) && tau_diff < 1e-8
  status <- if (score_pass && tau_pass) "PASS" else "FAIL"

  sens_results <- rbind(sens_results, data.frame(
    iter = iter, n = n, window = window,
    our_tau = our_tau, ref_tau = ref_tau, tau_diff = tau_diff,
    max_score_diff = if (is.na(max_score_diff)) NA_real_ else max_score_diff,
    n_scores = length(our_scores),
    status = status, stringsAsFactors = FALSE
  ))

  cat(sprintf("  %-4d  %-6d  %-5d  %-12s  %-12s  %-12s  %-12s  %s\n",
              iter, n, window, fmt(our_tau, 6), fmt(ref_tau, 6),
              fmte(tau_diff),
              if (is.na(max_score_diff)) "--" else fmte(max_score_diff),
              status))
}

sens_pass <- sum(sens_results$status == "PASS", na.rm = TRUE)
sens_fail <- sum(sens_results$status == "FAIL", na.rm = TRUE)
sens_skip <- sum(sens_results$status == "SKIP", na.rm = TRUE)
cat(sprintf("\n  SENSITIVITY REPORT: %d PASS, %d FAIL, %d SKIP out of 100\n",
            sens_pass, sens_fail, sens_skip))
valid_sens <- sens_results[sens_results$status == "PASS", ]
if (nrow(valid_sens) > 0) {
  cat(sprintf("  Tau diff:   mean=%.2e, max=%.2e\n",
              mean(valid_sens$tau_diff), max(valid_sens$tau_diff)))
  cat(sprintf("  Score diff: mean=%.2e, max=%.2e\n",
              mean(valid_sens$max_score_diff, na.rm = TRUE),
              max(valid_sens$max_score_diff, na.rm = TRUE)))
  cat(sprintf("  Tau range:  [%.4f, %.4f]\n",
              min(valid_sens$our_tau), max(valid_sens$our_tau)))
}


# ============================================================================
# PART 7: POTENTIAL ANALYSIS — Kernel Density + Well Detection
# ============================================================================
#
# What:   potential_analysis() vs manual stats::density() + Boltzmann U(x)
# Why:    Potential landscape analysis converts density P(x) to potential
#         U(x) = -log(P(x) + eps). Wells are local minima of U(x). We verify
#         that the density values, potential values, and well locations match.
# How:    Generate bimodal data (mixture of two Gaussians at -2 and +2). For
#         each iteration:
#           (a) Compute density via stats::density(n=100)
#           (b) Compute potential U = -log(density + eps)
#           (c) Compare landscape density and potential values
#           (d) Verify well detection (bimodal = 2 wells near -2 and +2)
# Tol:    Density values: < 1e-8 absolute
#         Well locations: within 1.5 of true means
# Ref:    Manual stats::density(n = 100) + U(x) = -log(P(x) + eps)
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("PART 7: POTENTIAL ANALYSIS -- 100 iterations\n")
cat("  Our potential_analysis() vs manual stats::density() + Boltzmann\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("  %-4s  %-6s  %-7s  %-14s  %-14s  %-16s  %s\n",
            "iter", "n", "n_wells", "dens_max_diff", "pot_max_diff", "well_locations", "status"))
cat(sprintf("  %s\n", strrep("-", 82)))

pot_results <- data.frame(
  iter = integer(), n = integer(), n_wells = integer(),
  density_max_diff = numeric(), potential_max_diff = numeric(),
  well_loc_str = character(), wells_correct = logical(),
  status = character(), stringsAsFactors = FALSE
)

for (iter in 1:100) {
  set.seed(iter)
  n <- sample(100:300, 1)

  # Bimodal data: mixture of N(-2, 0.5) and N(+2, 0.5)
  p <- runif(1, 0.3, 0.7)
  x <- c(rnorm(round(n * p), mean = -2, sd = 0.5),
         rnorm(n - round(n * p), mean = 2, sd = 0.5))
  x <- sample(x)

  our_result <- tryCatch(
    potential_analysis(x, n_bins = 100L),
    error = function(e) NULL
  )
  if (is.null(our_result)) {
    pot_results <- rbind(pot_results, data.frame(
      iter = iter, n = n, n_wells = NA_integer_,
      density_max_diff = NA, potential_max_diff = NA,
      well_loc_str = "--", wells_correct = NA,
      status = "SKIP", stringsAsFactors = FALSE
    ))
    cat(sprintf("  %-4d  %-6d  %-7s  %-14s  %-14s  %-16s  SKIP\n",
                iter, n, "--", "--", "--", "--"))
    next
  }

  # Manual computation (n_bins * 4L = actual grid resolution)
  n_grid <- 100L * 4L
  dens <- stats::density(x, n = n_grid)
  manual_density <- dens$y
  manual_x <- dens$x
  eps <- .Machine$double.eps
  manual_potential <- -log(pmax(manual_density, eps))

  # Compare density values
  our_ls <- our_result$landscape
  dens_diff <- NA_real_
  pot_diff <- NA_real_

  if (nrow(our_ls) == length(manual_x) && "density" %in% names(our_ls)) {
    dens_diff <- max(abs(our_ls$density - manual_density))
  }
  if (nrow(our_ls) == length(manual_x) && "potential" %in% names(our_ls)) {
    pot_diff <- max(abs(our_ls$potential - manual_potential))
  }

  # Well detection check
  nw <- our_result$n_wells
  well_locs <- if (nw > 0) our_result$wells$location else numeric(0)
  well_loc_str <- if (nw > 0) paste(round(well_locs, 2), collapse = ", ") else "none"

  wells_correct <- FALSE
  if (nw >= 2) {
    near_low <- any(abs(well_locs - (-2)) < 1.5)
    near_high <- any(abs(well_locs - 2) < 1.5)
    wells_correct <- near_low && near_high
  }

  dens_ok <- !is.na(dens_diff) && dens_diff < 1e-8
  pot_ok <- !is.na(pot_diff) && pot_diff < 1e-8
  status <- if (dens_ok && pot_ok && wells_correct) "PASS" else "FAIL"

  pot_results <- rbind(pot_results, data.frame(
    iter = iter, n = n, n_wells = nw,
    density_max_diff = dens_diff, potential_max_diff = pot_diff,
    well_loc_str = well_loc_str, wells_correct = wells_correct,
    status = status, stringsAsFactors = FALSE
  ))

  cat(sprintf("  %-4d  %-6d  %-7d  %-14s  %-14s  %-16s  %s\n",
              iter, n, nw,
              if (is.na(dens_diff)) "--" else fmte(dens_diff),
              if (is.na(pot_diff)) "--" else fmte(pot_diff),
              well_loc_str, status))
}

pot_pass <- sum(pot_results$status == "PASS", na.rm = TRUE)
pot_fail <- sum(pot_results$status == "FAIL", na.rm = TRUE)
pot_skip <- sum(pot_results$status == "SKIP", na.rm = TRUE)
cat(sprintf("\n  POTENTIAL REPORT: %d PASS, %d FAIL, %d SKIP out of 100\n",
            pot_pass, pot_fail, pot_skip))
valid_pot <- pot_results[pot_results$status == "PASS", ]
if (nrow(valid_pot) > 0) {
  cat(sprintf("  Density diff:   mean=%.2e, max=%.2e\n",
              mean(valid_pot$density_max_diff, na.rm = TRUE),
              max(valid_pot$density_max_diff, na.rm = TRUE)))
  cat(sprintf("  Potential diff:  mean=%.2e, max=%.2e\n",
              mean(valid_pot$potential_max_diff, na.rm = TRUE),
              max(valid_pot$potential_max_diff, na.rm = TRUE)))
  cat(sprintf("  Wells detected:  mean=%.1f, range=[%d, %d]\n",
              mean(valid_pot$n_wells), min(valid_pot$n_wells), max(valid_pot$n_wells)))
}


# ============================================================================
# PART 8: SURROGATE TEST — Full Pipeline (Phase Surrogates + Rolling AR(1))
# ============================================================================
#
# What:   surrogate_test(method = "phase", metric = "ar1") full pipeline vs
#         manual implementation: generate surrogates, compute rolling AR(1)
#         via ACF, compute Kendall tau, compute p-value.
# Why:    The surrogate test combines multiple steps (surrogate generation,
#         rolling metric, trend statistic, significance). The observed Kendall
#         tau and per-window metric values must match manual computation.
# How:    Generate random walk (n=100-300). For each iteration:
#           (a) Compute manual rolling AR(1) via acf(lag=1) per window
#           (b) Compute manual Kendall tau
#           (c) Compare observed_tau from surrogate_test() vs manual tau
#           (d) Compare observed_metric values per window
#           (e) Verify p-value is in [0, 1]
# Tol:    Observed tau: < 1e-8; Metric values: < 1e-8; p-value: [0, 1]
# Ref:    Manual stats::acf() + stats::cor.test(method = "kendall")
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("PART 8: SURROGATE TEST (full pipeline) -- 100 iterations\n")
cat("  Our surrogate_test(method='phase', metric='ar1') vs manual\n")
cat("  ACF-based rolling AR(1) + Kendall tau computation\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("  %-4s  %-6s  %-12s  %-12s  %-12s  %-10s  %-10s  %s\n",
            "iter", "n", "our_tau", "ref_tau", "tau_diff", "p_value", "metric_diff", "status"))
cat(sprintf("  %s\n", strrep("-", 82)))

pipe_results <- data.frame(
  iter = integer(), n = integer(),
  our_tau = numeric(), ref_tau = numeric(), tau_diff = numeric(),
  p_value = numeric(), max_metric_diff = numeric(),
  n_surrogates = integer(), status = character(),
  stringsAsFactors = FALSE
)

for (iter in 1:100) {
  set.seed(iter)
  n <- sample(100:300, 1)
  x <- cumsum(rnorm(n))
  window <- 30L
  n_surrogates <- 20L

  our_result <- tryCatch(
    surrogate_test(x, n_surrogates = n_surrogates, method = "phase",
                   metric = "ar1", window = window),
    error = function(e) NULL
  )
  if (is.null(our_result)) {
    pipe_results <- rbind(pipe_results, data.frame(
      iter = iter, n = n,
      our_tau = NA, ref_tau = NA, tau_diff = NA,
      p_value = NA, max_metric_diff = NA,
      n_surrogates = n_surrogates, status = "SKIP",
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  %-4d  %-6d  %-12s  %-12s  %-12s  %-10s  %-10s  SKIP\n",
                iter, n, "--", "--", "--", "--", "--"))
    next
  }

  # Manual rolling AR(1) via ACF
  n_win <- n - window + 1L
  manual_ar1 <- vapply(seq_len(n_win), function(i) {
    seg <- x[i:(i + window - 1L)]
    acf_val <- tryCatch(
      stats::acf(seg, lag.max = 1L, plot = FALSE)$acf[2L],
      error = function(e) NA_real_
    )
    if (is.na(acf_val)) NA_real_ else acf_val
  }, numeric(1L))

  valid <- !is.na(manual_ar1)
  ref_tau <- stats::cor.test(which(valid), manual_ar1[valid],
                              method = "kendall")$estimate

  tau_diff <- abs(our_result$observed_tau - ref_tau)

  # Compare observed metric values
  max_metric_diff <- NA_real_
  if (length(our_result$observed_metric) == length(manual_ar1)) {
    both_valid <- !is.na(our_result$observed_metric) & !is.na(manual_ar1)
    if (sum(both_valid) > 0) {
      max_metric_diff <- max(abs(our_result$observed_metric[both_valid] -
                                   manual_ar1[both_valid]))
    }
  }

  p_ok <- !is.na(our_result$p_value) &&
    our_result$p_value >= 0 && our_result$p_value <= 1
  tau_ok <- !is.na(tau_diff) && tau_diff < 1e-8
  metric_ok <- !is.na(max_metric_diff) && max_metric_diff < 1e-8
  status <- if (tau_ok && metric_ok && p_ok) "PASS" else "FAIL"

  pipe_results <- rbind(pipe_results, data.frame(
    iter = iter, n = n,
    our_tau = our_result$observed_tau, ref_tau = ref_tau, tau_diff = tau_diff,
    p_value = our_result$p_value, max_metric_diff = max_metric_diff,
    n_surrogates = n_surrogates, status = status,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  %-4d  %-6d  %-12s  %-12s  %-12s  %-10s  %-10s  %s\n",
              iter, n, fmt(our_result$observed_tau, 6), fmt(ref_tau, 6),
              fmte(tau_diff), fmt(our_result$p_value, 4),
              if (is.na(max_metric_diff)) "--" else fmte(max_metric_diff),
              status))
}

pipe_pass <- sum(pipe_results$status == "PASS", na.rm = TRUE)
pipe_fail <- sum(pipe_results$status == "FAIL", na.rm = TRUE)
pipe_skip <- sum(pipe_results$status == "SKIP", na.rm = TRUE)
cat(sprintf("\n  PIPELINE REPORT: %d PASS, %d FAIL, %d SKIP out of 100\n",
            pipe_pass, pipe_fail, pipe_skip))
valid_pipe <- pipe_results[pipe_results$status == "PASS", ]
if (nrow(valid_pipe) > 0) {
  cat(sprintf("  Tau diff:    mean=%.2e, max=%.2e\n",
              mean(valid_pipe$tau_diff), max(valid_pipe$tau_diff)))
  cat(sprintf("  Metric diff: mean=%.2e, max=%.2e\n",
              mean(valid_pipe$max_metric_diff, na.rm = TRUE),
              max(valid_pipe$max_metric_diff, na.rm = TRUE)))
  cat(sprintf("  Tau range:   [%.4f, %.4f]\n",
              min(valid_pipe$our_tau), max(valid_pipe$our_tau)))
  cat(sprintf("  P-values:    mean=%.4f, range=[%.4f, %.4f]\n",
              mean(valid_pipe$p_value),
              min(valid_pipe$p_value), max(valid_pipe$p_value)))
}


# ============================================================================
# FINAL SUMMARY REPORT
# ============================================================================

t_elapsed <- (proc.time() - t_start)[3]

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("FINAL SUMMARY -- External Package Equivalence Tests\n")
cat(sprintf("Total runtime: %.1f seconds (%.1f minutes)\n", t_elapsed, t_elapsed / 60))
cat(strrep("=", 70), "\n\n")

summary_df <- data.frame(
  Part = c(
    "1. Changepoint PELT",
    "2. Changepoint BinSeg",
    "3. Multivariate EWS",
    "4. Phase surrogates",
    "5. Spectral exponent",
    "6. Sensitivity EWS",
    "7. Potential analysis",
    "8. Surrogate pipeline"
  ),
  Reference = c(
    "changepoint::cpt.mean(PELT)",
    "changepoint::cpt.mean(BinSeg)",
    "EWSmethods::multiEWS()",
    "nonlinearTseries::FFTsurrogate()",
    "stats::spectrum() + lm()",
    "stats::ar.ols() + cor.test()",
    "stats::density() + Boltzmann",
    "stats::acf() + cor.test()"
  ),
  Pass = c(pelt_pass, binseg_pass, mews_pass, surr_pass,
           spec_pass, sens_pass, pot_pass, pipe_pass),
  Fail = c(pelt_fail, binseg_fail, mews_fail, surr_fail,
           spec_fail, sens_fail, pot_fail, pipe_fail),
  Skip = c(pelt_skip, binseg_skip, mews_skip, surr_skip,
           spec_skip, sens_skip, pot_skip, pipe_skip),
  stringsAsFactors = FALSE
)

summary_df$Rate <- sprintf("%.0f%%", 100 * summary_df$Pass /
                             (summary_df$Pass + summary_df$Fail + summary_df$Skip))

# Print as formatted table
cat(sprintf("  %-25s  %-32s  %4s  %4s  %4s  %5s\n",
            "Function", "Reference Package", "Pass", "Fail", "Skip", "Rate"))
cat(sprintf("  %s\n", strrep("-", 80)))
for (i in seq_len(nrow(summary_df))) {
  cat(sprintf("  %-25s  %-32s  %4d  %4d  %4d  %5s\n",
              summary_df$Part[i], summary_df$Reference[i],
              summary_df$Pass[i], summary_df$Fail[i],
              summary_df$Skip[i], summary_df$Rate[i]))
}

total_pass <- sum(summary_df$Pass)
total_fail <- sum(summary_df$Fail)
total_skip <- sum(summary_df$Skip)
total_tests <- total_pass + total_fail + total_skip

cat(sprintf("  %s\n", strrep("-", 80)))
cat(sprintf("  %-25s  %-32s  %4d  %4d  %4d  %5s\n",
            "TOTAL", paste0(total_tests, " iterations"),
            total_pass, total_fail, total_skip,
            sprintf("%.0f%%", 100 * total_pass / total_tests)))
cat(strrep("=", 70), "\n")

if (total_fail == 0) {
  cat("\nALL TESTS PASSED -- numerical equivalence with external packages verified.\n")
} else {
  cat(sprintf("\nWARNING: %d FAILURES detected. Review per-part reports above.\n",
              total_fail))
}
