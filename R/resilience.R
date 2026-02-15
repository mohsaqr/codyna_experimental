# ============================================================================
# resilience.R --Resilience metrics, classification, and visualization
# Metric computations are numerically identical to Resilience4_Functions/
# 05_resilience_metrics_functions.R (calculate_rolling_resilience_metrics).
# ============================================================================

# --------------------------------------------------------------------------
# Exported: resilience()
# --------------------------------------------------------------------------

#' Resilience Metrics for Time Series Data
#'
#' Computes rolling-window resilience metrics for univariate time series data,
#' quantifying stability, volatility, variability, memory, complexity,
#' and recovery characteristics. Metric computations are identical to
#' `calculate_rolling_resilience_metrics()` from the Resilience4 toolkit.
#'
#' @details
#' The function slides a fixed-size window across the time series. Within each
#' window position, eight resilience-relevant metrics are computed on the local
#' segment. The resulting metric trajectories reveal *when* and *how* resilience
#' changes over time, not just whether the system is resilient on average.
#'
#' The eight metrics and their mathematical definitions are:
#'
#' \describe{
#'   \item{vsi (Variance of Sub-window Variances)}{The variance of rolling
#'     variances computed on sub-windows within each main window
#'     (\eqn{\mathrm{Var}(\sigma^2_{w_1}, \ldots, \sigma^2_{w_k})}). Captures
#'     instability in volatility. Lower values indicate more stable variance
#'     and thus better resilience.}
#'   \item{arch_lm (ARCH-LM Statistic)}{The Engle (1982) test statistic
#'     \eqn{n \cdot R^2} from regressing squared residuals on their own lags.
#'     Detects autoregressive conditional heteroscedasticity (volatility
#'     clustering). Lower values indicate more homogeneous variance and better
#'     resilience.}
#'   \item{cv (Coefficient of Variation)}{\eqn{\mathrm{CV} = \sigma / |\mu|},
#'     the ratio of standard deviation to the absolute mean. Measures relative
#'     dispersion. Lower values indicate tighter fluctuations around the mean
#'     and better resilience.}
#'   \item{recovery_time}{Median number of time steps required to return within
#'     a threshold of the baseline after a shock event. Shocks are detected as
#'     deviations exceeding `shock_threshold` baseline standard deviations.
#'     Lower values indicate faster recovery and better resilience.}
#'   \item{recovery_slope}{Median OLS slope of the series over a short window
#'     immediately following each detected shock. Positive slopes indicate
#'     upward recovery toward baseline; the sign and magnitude together convey
#'     recovery direction and speed. Higher (more positive) values indicate
#'     better resilience.}
#'   \item{sample_entropy}{\eqn{\mathrm{SampEn} = -\ln(A / B)}, where \eqn{A}
#'     is the number of template matches at embedding dimension \eqn{m + 1} and
#'     \eqn{B} is the number of matches at dimension \eqn{m} (Richman &
#'     Moorman, 2000). Quantifies signal complexity and regularity. Higher
#'     values indicate greater unpredictability and better resilience (overly
#'     regular dynamics may signal rigidity).}
#'   \item{dfa_alpha (DFA Scaling Exponent)}{The slope of the log-log
#'     fluctuation plot from Detrended Fluctuation Analysis (Peng et al.,
#'     1994). Values near 0.5 indicate uncorrelated noise; values near 1.0
#'     indicate \eqn{1/f} noise; values near 1.5 indicate Brownian motion.
#'     Lower values (closer to white noise) indicate better resilience.}
#'   \item{ac_ratio (Autocorrelation Ratio)}{\eqn{|\mathrm{acf}[1]| /
#'     \mathrm{mean}(|\mathrm{acf}[2{:}k]|)}, the ratio of lag-1
#'     autocorrelation magnitude to the mean magnitude of higher-order lags.
#'     High values suggest strong short-term memory with faster decorrelation.
#'     Higher values indicate better resilience.}
#' }
#'
#' Note that metrics are directional: for some (vsi, arch_lm, cv,
#' recovery_time, dfa_alpha) lower values signal better resilience, while for
#' others (recovery_slope, sample_entropy, ac_ratio) higher values are
#' preferable. The [classify_resilience()] function handles this directionality
#' automatically.
#'
#' @export
#' @param data \[`ts`, `numeric()`, `data.frame`\]\cr
#'   Univariate time series data.
#' @param window \[`integer(1)`: `50L`\]\cr
#'   A positive `integer` specifying the rolling window size.
#'   Must be at least `20`.
#' @param align \[`character(1)`: `"right"`\]\cr
#'   Alignment of the window. The available options are:
#'   `"right"`, `"center"`, and `"left"`.
#' @param metrics \[`character()`\]\cr
#'   A vector of metrics to calculate, or `"all"` for all eight metrics.
#' @param shock_threshold \[`numeric(1)`: `2`\]\cr
#'   Number of baseline standard deviations for shock detection.
#' @param recovery_threshold \[`numeric(1)`: `0.1`\]\cr
#'   Fraction of baseline SD within which recovery is declared.
#' @param recovery_window \[`integer(1)`: `5L`\]\cr
#'   Number of points after a shock to fit recovery slope.
#' @param baseline_window \[`integer(1)`: `10L`\]\cr
#'   Number of initial points used to establish baseline.
#' @param edim \[`integer(1)`: `2L`\]\cr
#'   Embedding dimension for sample entropy.
#' @param r \[`numeric(1)`: `NULL`\]\cr
#'   Tolerance for sample entropy template matching.
#'   If `NULL`, uses `0.2 * sd(x)`.
#' @param tau \[`integer(1)`: `1L`\]\cr
#'   Time delay for sample entropy embedding.
#' @param max_lag \[`integer(1)`: `5L`\]\cr
#'   Maximum lag for autocorrelation ratio.
#' @param lags \[`integer(1)`: `1L`\]\cr
#'   Number of lags for ARCH-LM test.
#' @param demean \[`logical(1)`: `FALSE`\]\cr
#'   Whether to demean the series before ARCH-LM computation.
#' @return An object of class `"resilience"` (inherits from
#'   [tibble::tbl_df]) with the following columns:
#'   \itemize{
#'     \item `time` -- integer or Date index of each observation.
#'     \item `value` -- the original time series values.
#'     \item `vsi` -- Variance of Sub-window Variances (if requested).
#'     \item `arch_lm` -- ARCH-LM statistic (if requested).
#'     \item `cv` -- Coefficient of Variation (if requested).
#'     \item `recovery_time` -- median recovery time in steps (if requested).
#'     \item `recovery_slope` -- median post-shock OLS slope (if requested).
#'     \item `sample_entropy` -- Sample Entropy (if requested).
#'     \item `dfa_alpha` -- DFA scaling exponent (if requested).
#'     \item `ac_ratio` -- Autocorrelation Ratio (if requested).
#'   }
#'
#'   The returned object carries three attributes:
#'   \itemize{
#'     \item `window` -- the window size used.
#'     \item `align` -- the alignment method (`"right"`, `"center"`, or
#'       `"left"`).
#'     \item `metrics` -- character vector of computed metric names.
#'   }
#'
#'   Leading or trailing rows contain `NA` values for the metrics depending on
#'   the alignment setting, because the rolling window cannot be computed at
#'   those positions.
#'
#' @references
#' Richman, J. S. & Moorman, J. R. (2000). Physiological time-series analysis
#' using approximate entropy and sample entropy. *American Journal of
#' Physiology-Heart and Circulatory Physiology*, 278(6), H2039--H2049.
#' \doi{10.1152/ajpheart.2000.278.6.H2039}
#'
#' Peng, C.-K., Buldyrev, S. V., Havlin, S., Simons, M., Stanley, H. E., &
#' Goldberger, A. L. (1994). Mosaic organization of DNA nucleotides.
#' *Physical Review E*, 49(2), 1685--1689.
#' \doi{10.1103/PhysRevE.49.1685}
#'
#' Engle, R. F. (1982). Autoregressive Conditional Heteroscedasticity with
#' Estimates of the Variance of United Kingdom Inflation. *Econometrica*,
#' 50(4), 987--1007. \doi{10.2307/1912773}
#'
#' @seealso [classify_resilience()] for directional normalization and
#'   categorical scoring; [plot.resilience()] for visualization.
#' @family resilience
#' @concept time series
#' @concept resilience
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' res <- resilience(x, window = 50L)
#' res
#'
#' # Visualize raw metric trajectories
#' plot(res, type = "lines")
#' }
resilience <- function(data,
                       window = 50L,
                       align = "right",
                       metrics = "all",
                       shock_threshold = 2,
                       recovery_threshold = 0.1,
                       recovery_window = 5L,
                       baseline_window = 10L,
                       edim = 2L,
                       r = NULL,
                       tau = 1L,
                       max_lag = 5L,
                       lags = 1L,
                       demean = FALSE) {
  check_missing(data)
  data <- prepare_timeseries_data(data)
  values <- data$values
  time <- data$time
  n <- length(values)
  align <- check_match(align, c("right", "center", "left"))
  available_metrics <- c(
    "vsi", "arch_lm", "cv", "recovery_time", "recovery_slope",
    "sample_entropy", "dfa_alpha", "ac_ratio"
  )
  metrics <- check_match(
    metrics,
    c(available_metrics, "all"),
    several.ok = TRUE
  )
  metrics <- ifelse_("all" %in% metrics, available_metrics, metrics)
  check_range(window, type = "integer", min = 20L, max = n)

  # --- Rolling window computation (matches Resilience4 exactly) ---
  total_windows <- n - window + 1L
  results_list <- vector("list", total_windows)

  for (i in seq_len(total_windows)) {
    window_data <- values[i:(i + window - 1L)]
    window_data <- window_data[!is.na(window_data)]
    wd_len <- length(window_data)

    if (wd_len < (window / 2)) {
      row <- stats::setNames(
        as.list(rep(NA_real_, length(metrics))),
        metrics
      )
    } else {
      # Adaptive DFA parameters (matching Resilience4)
      adaptive_min_scale <- max(3L, wd_len %/% 20L)
      adaptive_max_scale <- max(adaptive_min_scale + 2L, wd_len %/% 4L)

      row <- list()
      for (m in metrics) {
        row[[m]] <- switch(
          m,
          vsi            = resilience_vsi(window_data),
          arch_lm        = resilience_arch_lm(window_data, lags, demean),
          cv             = resilience_cv(window_data),
          recovery_time  = resilience_recovery_time(
            window_data, shock_threshold, recovery_threshold, baseline_window
          ),
          recovery_slope = resilience_recovery_slope(
            window_data, shock_threshold, recovery_window, baseline_window
          ),
          sample_entropy = resilience_sample_entropy(
            window_data, edim, r, tau
          ),
          dfa_alpha      = resilience_dfa_alpha(
            window_data,
            min_scale     = adaptive_min_scale,
            max_scale     = adaptive_max_scale,
            min_obs       = max(15L, wd_len %/% 3L),
            min_r_squared = 0.3
          ),
          ac_ratio       = resilience_ac_ratio(window_data, max_lag)
        )
      }
    }
    results_list[[i]] <- row
  }

  # Convert to data frame
  metrics_df <- do.call(
    rbind,
    lapply(results_list, function(r) as.data.frame(r, stringsAsFactors = FALSE))
  )

  # --- Alignment padding ---
  pad_start <- 0L
  pad_end <- 0L
  if (align == "right") {
    pad_start <- window - 1L
  } else if (align == "center") {
    pad_start <- floor((window - 1L) / 2)
    pad_end <- ceiling((window - 1L) / 2)
  } else if (align == "left") {
    pad_end <- window - 1L
  }

  if (pad_start > 0L) {
    na_pad <- data.frame(
      matrix(NA_real_, nrow = pad_start, ncol = ncol(metrics_df))
    )
    names(na_pad) <- names(metrics_df)
    metrics_df <- rbind(na_pad, metrics_df)
  }
  if (pad_end > 0L) {
    na_pad <- data.frame(
      matrix(NA_real_, nrow = pad_end, ncol = ncol(metrics_df))
    )
    names(na_pad) <- names(metrics_df)
    metrics_df <- rbind(metrics_df, na_pad)
  }

  out <- tibble::as_tibble(
    cbind(data.frame(time = time, value = values), metrics_df)
  )
  structure(
    out,
    window = window,
    align = align,
    metrics = metrics,
    class = c("resilience", "tbl_df", "tbl", "data.frame")
  )
}

# --------------------------------------------------------------------------
# Internal metric helpers (@noRd)
# --------------------------------------------------------------------------

#' @noRd
resilience_vsi <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 10L) return(NA_real_)
  window_size <- max(5L, floor(n / 4))
  if (n < window_size * 2L) return(stats::var(x, na.rm = TRUE))
  n_windows <- n - window_size + 1L
  rolling_vars <- numeric(n_windows)
  for (i in seq_len(n_windows)) {
    rolling_vars[i] <- stats::var(x[i:(i + window_size - 1L)], na.rm = TRUE)
  }
  stats::var(rolling_vars, na.rm = TRUE)
}

#' @noRd
resilience_arch_lm <- function(x, lags = 1L, demean = FALSE) {
  x <- as.numeric(x[!is.na(x)])
  n <- length(x)
  if (n < max(10L, lags + 3L)) return(NA_real_)
  if (demean) x <- x - mean(x)
  x2 <- x^2
  X <- matrix(NA_real_, nrow = n, ncol = lags)
  for (i in seq_len(lags)) {
    X[(i + 1L):n, i] <- x2[1L:(n - i)]
  }
  y <- x2[(lags + 1L):n]
  X_clean <- X[(lags + 1L):n, , drop = FALSE]
  if (length(y) < 3L) return(NA_real_)
  fit <- try_(stats::lm(y ~ X_clean))
  if (inherits(fit, "try-error")) return(NA_real_)
  r_sq <- summary(fit)$r.squared
  length(y) * r_sq
}

#' @noRd
resilience_cv <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3L) return(NA_real_)
  mu <- mean(x)
  if (abs(mu) < .Machine$double.eps) return(NA_real_)
  stats::sd(x) / abs(mu)
}

#' @noRd
resilience_recovery_time <- function(x,
                                     shock_threshold = 2,
                                     recovery_threshold = 0.1,
                                     baseline_window = 10L) {
  x <- x[!is.na(x) & is.finite(x)]
  n <- length(x)
  if (n < 20L) return(NA_real_)
  baseline_data <- x[seq_len(min(baseline_window, n))]
  if (length(baseline_data) < 3L) return(NA_real_)
  baseline <- mean(baseline_data)
  baseline_sd <- stats::sd(baseline_data)
  if (is.na(baseline_sd) || baseline_sd == 0) return(NA_real_)
  deviations <- abs(x - baseline)
  shock_points <- which(deviations > shock_threshold * baseline_sd)
  if (length(shock_points) == 0L) return(0)
  recovery_times <- numeric(length(shock_points))
  for (i in seq_along(shock_points)) {
    shock_idx <- shock_points[i]
    recovery_threshold_value <- recovery_threshold * baseline_sd
    if (shock_idx < n) {
      remaining <- abs(x[(shock_idx + 1L):n] - baseline)
      recovered_indices <- which(remaining <= recovery_threshold_value)
      if (length(recovered_indices) > 0L) {
        recovery_times[i] <- recovered_indices[1L]
      } else {
        recovery_times[i] <- n - shock_idx
      }
    } else {
      recovery_times[i] <- 0
    }
  }
  stats::median(recovery_times, na.rm = TRUE)
}

#' @noRd
resilience_recovery_slope <- function(x,
                                      shock_threshold = 2,
                                      recovery_window = 5L,
                                      baseline_window = 10L) {
  x <- x[!is.na(x) & is.finite(x)]
  n <- length(x)
  if (n < 20L) return(NA_real_)
  baseline_data <- x[seq_len(min(baseline_window, n))]
  if (length(baseline_data) < 3L) return(NA_real_)
  baseline <- mean(baseline_data)
  baseline_sd <- stats::sd(baseline_data)
  if (is.na(baseline_sd) || baseline_sd == 0) return(NA_real_)
  deviations <- abs(x - baseline)
  shock_points <- which(deviations > shock_threshold * baseline_sd)
  if (length(shock_points) == 0L) return(0)
  recovery_slopes <- numeric(length(shock_points))
  for (i in seq_along(shock_points)) {
    shock_idx <- shock_points[i]
    end_idx <- min(shock_idx + recovery_window, n)
    if (end_idx - shock_idx < 2L) {
      recovery_slopes[i] <- NA_real_
      next
    }
    time_points <- shock_idx:end_idx
    values <- x[time_points]
    valid_idx <- !is.na(values) & is.finite(values)
    if (sum(valid_idx) < 2L) {
      recovery_slopes[i] <- NA_real_
      next
    }
    fit <- try_(stats::lm(values[valid_idx] ~ time_points[valid_idx]))
    if (inherits(fit, "try-error")) {
      recovery_slopes[i] <- NA_real_
    } else {
      recovery_slopes[i] <- stats::coef(fit)[2L]
    }
  }
  stats::median(recovery_slopes, na.rm = TRUE)
}

#' @noRd
resilience_sample_entropy <- function(x, edim = 2L, r = NULL, tau = 1L) {
  x <- as.numeric(x[!is.na(x)])
  n <- length(x)
  if (n < 10L) return(NA_real_)
  if (is.null(r)) r <- 0.2 * stats::sd(x)
  if (is.na(r) || r < 0 || !is.finite(r)) return(NA_real_)
  if (n <= edim * tau) return(NA_real_)
  # Create embedded vectors
  create_embedded <- function(data, m, delay) {
    nn <- length(data)
    num_vec <- nn - (m - 1L) * delay
    if (num_vec <= 0L) return(matrix(nrow = 0L, ncol = m))
    emb <- matrix(NA_real_, nrow = num_vec, ncol = m)
    for (i in seq_len(num_vec)) {
      indices <- i + (0L:(m - 1L)) * delay
      emb[i, ] <- data[indices]
    }
    emb
  }
  # Count template matches (Chebyshev distance <= r)
  count_matches <- function(vectors, tol) {
    nv <- nrow(vectors)
    if (nv < 2L) return(0L)
    matches <- 0L
    for (i in seq_len(nv - 1L)) {
      for (j in (i + 1L):nv) {
        if (max(abs(vectors[i, ] - vectors[j, ])) <= tol) {
          matches <- matches + 1L
        }
      }
    }
    matches
  }
  vectors_m <- create_embedded(x, edim, tau)
  vectors_m1 <- create_embedded(x, edim + 1L, tau)
  n_m1 <- nrow(vectors_m1)
  if (n_m1 < 2L) return(NA_real_)
  vectors_m_trunc <- vectors_m[seq_len(n_m1), , drop = FALSE]
  A <- count_matches(vectors_m1, r)
  B <- count_matches(vectors_m_trunc, r)
  -log(A / B)
}

#' @noRd
resilience_dfa_alpha <- function(x,
                                 min_scale = 4L,
                                 max_scale = NULL,
                                 scales_per_octave = 8L,
                                 poly_order = 1L,
                                 min_obs = 30L,
                                 min_r_squared = 0.6) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < min_obs) return(NA_real_)
  if (is.null(max_scale)) max_scale <- floor(n / 6)
  if (min_scale >= max_scale) return(NA_real_)
  # Integrated profile
  profile <- cumsum(x - mean(x, na.rm = TRUE))
  # Log-spaced scales
  scales <- unique(round(
    exp(seq(log(min_scale), log(max_scale), length.out = 20L))
  ))
  scales <- scales[scales >= min_scale & scales <= max_scale]
  if (length(scales) < 4L) return(NA_real_)
  fluct <- numeric(length(scales))
  for (si in seq_along(scales)) {
    s <- scales[si]
    n_seg <- floor(n / s)
    if (n_seg < 3L) {
      fluct[si] <- NA_real_
      next
    }
    all_residuals <- c()
    for (seg in seq_len(n_seg)) {
      start_idx <- (seg - 1L) * s + 1L
      end_idx <- seg * s
      segment <- profile[start_idx:end_idx]
      tp <- seq_along(segment)
      fit <- try_(stats::lm(
        segment ~ stats::poly(tp, degree = poly_order, raw = TRUE)
      ))
      if (!inherits(fit, "try-error")) {
        all_residuals <- c(all_residuals, stats::residuals(fit))
      }
    }
    if (length(all_residuals) > 0L) {
      fluct[si] <- sqrt(mean(all_residuals^2))
    } else {
      fluct[si] <- NA_real_
    }
  }
  valid <- !is.na(fluct) & fluct > 0
  if (sum(valid) < 4L) return(NA_real_)
  log_s <- log10(scales[valid])
  log_f <- log10(fluct[valid])
  fit <- try_(stats::lm(log_f ~ log_s))
  if (inherits(fit, "try-error")) return(NA_real_)
  alpha <- stats::coef(fit)[2L]
  r_sq <- summary(fit)$r.squared
  if (r_sq < min_r_squared) {
    warning_(
      "Poor linear fit in DFA (R^2 = {round(r_sq, 3)})."
    )
  }
  alpha
}

#' @noRd
resilience_ac_ratio <- function(x, max_lag = 5L) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 10L) return(NA_real_)
  if (n <= max_lag) return(NA_real_)
  ac_result <- try_(stats::acf(x, lag.max = max_lag, plot = FALSE,
                                na.action = stats::na.pass))
  if (inherits(ac_result, "try-error")) return(NA_real_)
  ac_vals <- as.numeric(ac_result$acf)[-1L]
  if (length(ac_vals) < 2L) return(NA_real_)
  ac_lag1 <- ac_vals[1L]
  ac_higher <- mean(abs(ac_vals[-1L]), na.rm = TRUE)
  if (is.na(ac_higher) || ac_higher == 0) return(Inf)
  abs(ac_lag1) / ac_higher
}

# --------------------------------------------------------------------------
# Exported: classify_resilience()
# --------------------------------------------------------------------------

#' Classify Resilience Scores
#'
#' Applies directional normalization to resilience metrics and computes
#' a weighted composite score with category labels.
#'
#' @details
#' Raw metric values from [resilience()] are not directly comparable because
#' they live on different scales and some are "lower is better" while others
#' are "higher is better." This function resolves both problems in three
#' steps:
#'
#' **1. Directional normalization.** Each metric is rescaled to a \[0, 1\]
#' score where 0 means "Excellent resilience" and 1 means "Troubled
#' resilience," regardless of the original direction. Metrics where lower raw
#' values indicate better resilience (vsi, arch_lm, cv, recovery_time,
#' dfa_alpha) are scaled so that low raw values map to scores near 0.
#' Metrics where higher raw values indicate better resilience
#' (recovery_slope, sample_entropy, ac_ratio) are inverted accordingly.
#'
#' The `"empirical_state_aware"` method (the default) derives scaling
#' thresholds from empirically observed means for known resilience phases
#' (Stable, Recovery, Shock, Turbulent, Volatile). "Excellent" and
#' "Troubled" anchor points are set from quantiles of the good-state and
#' bad-state phase means, respectively. Alternative methods
#' (`"empirical_percentile"`, `"percentile"`, `"minmax"`, `"zscore"`,
#' `"divide_max"`) use purely data-driven scaling.
#'
#' **2. Smoothing.** After scaling, each metric's score trajectory is
#' smoothed with a rolling median (or mean) of width `smooth_window` to
#' reduce transient noise in the classification. Set `smooth_window = 1L`
#' to disable smoothing.
#'
#' **3. Weighted composite.** The per-metric scores are combined into a
#' single `composite_score` using empirically derived weights that reflect
#' each metric's importance for distinguishing resilience phases.
#' The default weights are: vsi = 0.256, cv = 0.226, arch_lm = 0.188,
#' ac_ratio = 0.070, recovery_slope = 0.058, dfa_alpha = 0.036,
#' sample_entropy = 0.028, recovery_time = 0.016.
#'
#' Both per-metric and composite scores are binned into six ordered
#' categories:
#' \itemize{
#'   \item \strong{Excellent} (0.00--0.15)
#'   \item \strong{Solid} (0.15--0.30)
#'   \item \strong{Fair} (0.30--0.50)
#'   \item \strong{Vulnerable} (0.50--0.70)
#'   \item \strong{Failing} (0.70--0.85)
#'   \item \strong{Troubled} (0.85--1.00)
#' }
#'
#' @export
#' @param data \[`resilience`\]\cr
#'   An object of class `resilience` as returned by [resilience()].
#' @param method \[`character(1)`: `"empirical_state_aware"`\]\cr
#'   Normalization method for scaling raw metric values. The available
#'   options are:
#'   * `"empirical_state_aware"`: State-aware scaling using empirical
#'     thresholds derived from known resilience phases.
#'   * `"empirical_percentile"`: Hardcoded percentile-bin scaling.
#'   * `"percentile"`: Empirical CDF (rank-based) scaling.
#'   * `"minmax"`: Min-max scaling to \[0, 1\].
#'   * `"zscore"`: Z-score standardization.
#'   * `"divide_max"`: Division by absolute maximum.
#' @param weights \[`numeric()`: `NULL`\]\cr
#'   Named numeric vector of metric weights for the composite score.
#'   If `NULL`, empirically derived default weights are used.
#' @param smooth_window \[`integer(1)`: `10L`\]\cr
#'   Window size for post-scaling median smoothing (applied per-metric).
#'   Set to `1L` to disable smoothing.
#' @param smooth_method \[`character(1)`: `"median"`\]\cr
#'   Smoothing method: `"median"` or `"mean"`.
#' @return The input `"resilience"` object augmented with additional columns:
#'   \itemize{
#'     \item `{metric}_score` -- directionally normalized score in \[0, 1\]
#'       for each computed metric (e.g., `vsi_score`, `cv_score`).
#'     \item `{metric}_category` -- ordered factor label
#'       (`"Excellent"` through `"Troubled"`) for each metric.
#'     \item `composite_score` -- weighted average of all per-metric scores.
#'     \item `composite_category` -- ordered factor label for the composite
#'       score.
#'   }
#'
#'   All original columns (`time`, `value`, and the raw metric columns) are
#'   preserved. The returned object retains the `"resilience"` class and
#'   its attributes.
#'
#' @seealso [resilience()] for computing the raw metrics;
#'   [plot.resilience()] for visualization of classified results.
#' @family resilience
#' @concept time series
#' @concept resilience
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' res <- resilience(x, window = 50L)
#' cls <- classify_resilience(res)
#' cls
#'
#' # Ribbon plot with classified scores
#' plot(cls)
#' }
classify_resilience <- function(data,
                                method = "empirical_state_aware",
                                weights = NULL,
                                smooth_window = 10L,
                                smooth_method = "median") {
  check_missing(data)
  check_class(data, "resilience")
  method <- check_match(
    method,
    c("empirical_state_aware", "empirical_percentile",
      "percentile", "minmax", "zscore", "divide_max")
  )
  smooth_method <- check_match(smooth_method, c("median", "mean"))
  metrics <- attr(data, "metrics")
  # Set up weights (NOT pre-normalized --matches original behavior)
  default_w <- default_resilience_weights()
  if (is.null(weights)) {
    weights <- default_w[metrics]
  } else {
    stopifnot_(
      is.numeric(weights) && length(weights) == length(metrics),
      "Argument {.arg weights} must be a named numeric vector of
      length {length(metrics)}."
    )
    if (is.null(names(weights))) names(weights) <- metrics
  }
  # Category breaks and labels
  breaks <- c(0, 0.15, 0.30, 0.50, 0.70, 0.85, 1.0)
  labels <- c(
    "Excellent", "Solid", "Fair", "Vulnerable", "Failing", "Troubled"
  )
  # Scale each metric directionally
  n <- nrow(data)
  weighted_mat <- matrix(0, nrow = n, ncol = length(metrics))
  total_weight <- 0
  for (mi in seq_along(metrics)) {
    m <- metrics[mi]
    raw <- data[[m]]
    scored <- scale_directional(raw, m, method)
    # Clamp to [0, 1]
    scored <- pmax(0, pmin(1, scored))
    # Apply smoothing
    if (smooth_window > 1L && length(scored) >= smooth_window) {
      scored <- resilience_smooth_(scored, smooth_window, smooth_method)
    }
    data[[paste0(m, "_score")]] <- scored
    data[[paste0(m, "_category")]] <- cut(
      scored,
      breaks = breaks,
      labels = labels,
      include.lowest = TRUE
    )
    w <- if (m %in% names(weights)) weights[[m]] else 0.01
    weighted_mat[, mi] <- ifelse(is.na(scored), 0, scored * w)
    total_weight <- total_weight + w
  }
  data$composite_score <- rowSums(weighted_mat, na.rm = TRUE) / total_weight
  data$composite_category <- cut(
    data$composite_score,
    breaks = breaks,
    labels = labels,
    include.lowest = TRUE
  )
  structure(
    data,
    class = c("resilience", "tbl_df", "tbl", "data.frame")
  )
}

# --------------------------------------------------------------------------
# Internal classification helpers (@noRd)
# --------------------------------------------------------------------------

#' @noRd
is_lower_better <- function(metric) {
  metric %in% c("vsi", "arch_lm", "cv", "recovery_time", "dfa_alpha")
}

#' @noRd
default_resilience_weights <- function() {
  c(
    vsi             = 0.256,
    cv              = 0.226,
    arch_lm         = 0.188,
    ac_ratio        = 0.070,
    recovery_slope  = 0.058,
    dfa_alpha       = 0.036,
    sample_entropy  = 0.028,
    recovery_time   = 0.016
  )
}

#' @noRd
resilience_empirical_means <- function() {
  data.frame(
    true_phase      = c("Recovery", "Shock", "Stable", "Turbulent", "Volatile"),
    vsi             = c(24.65, 75.61, 1455.42, 2857.7, 15764.89),
    arch_lm         = c(18.05, 17.97, 9.47, 16.9, 2.03),
    cv              = c(0.13, 0.16, 0.07, 0.16, 0.16),
    recovery_time   = c(6.57, 6.35, 6.19, 6.17, 3.92),
    recovery_slope  = c(0.61, -0.83, -0.58, 0.51, 0.28),
    sample_entropy  = c(1.75, 0.93, 1.53, 1.57, 1.6),
    dfa_alpha       = c(0.97, 0.86, 1.3, 1.45, 1.04),
    ac_ratio        = c(1.39, 1.33, 2.8, 2.29, 1.79),
    stringsAsFactors = FALSE
  )
}

#' @noRd
scale_directional <- function(values, metric, method) {
  # Handle Inf in sample_entropy (match original behavior)
  if (metric == "sample_entropy") {
    inf_mask <- is.infinite(values)
    if (any(inf_mask)) {
      finite_vals <- values[!inf_mask & !is.na(values)]
      if (length(finite_vals) > 0L) {
        values[inf_mask] <- max(finite_vals) * 1.1
      }
    }
  }
  valid <- values[!is.na(values) & is.finite(values)]
  if (length(valid) < 2L) return(rep(0.5, length(values)))
  lower <- is_lower_better(metric)
  scaled <- switch(
    method,
    empirical_state_aware = {
      empirical <- resilience_empirical_means()
      good_states  <- c("Stable", "Recovery")
      bad_states   <- c("Shock")
      good_vals <- empirical[[metric]][empirical$true_phase %in% good_states]
      bad_vals  <- empirical[[metric]][empirical$true_phase %in% bad_states]
      good_vals <- good_vals[!is.na(good_vals) & is.finite(good_vals)]
      bad_vals  <- bad_vals[!is.na(bad_vals) & is.finite(bad_vals)]
      if (length(good_vals) < 1L || length(bad_vals) < 1L) {
        # Fall back to percentile
        sc <- scale_directional(values, metric, "empirical_percentile")
        return(sc)
      }
      if (lower) {
        excellent_th <- stats::quantile(good_vals, 0.75, na.rm = TRUE)
        troubled_th  <- stats::quantile(bad_vals, 0.25, na.rm = TRUE)
      } else {
        excellent_th <- stats::quantile(good_vals, 0.25, na.rm = TRUE)
        troubled_th  <- stats::quantile(bad_vals, 0.75, na.rm = TRUE)
      }
      sc <- numeric(length(values))
      for (i in seq_along(values)) {
        v <- values[i]
        if (is.na(v) || !is.finite(v)) { sc[i] <- NA_real_; next }
        if (lower) {
          if (v <= excellent_th) {
            sc[i] <- 0
          } else if (v >= troubled_th) {
            sc[i] <- 1
          } else {
            rng <- troubled_th - excellent_th
            sc[i] <- if (rng > 0) (v - excellent_th) / rng else 0.5
          }
        } else {
          if (v >= excellent_th) {
            sc[i] <- 0
          } else if (v <= troubled_th) {
            sc[i] <- 1
          } else {
            rng <- excellent_th - troubled_th
            sc[i] <- if (rng > 0) 1 - (v - troubled_th) / rng else 0.5
          }
        }
      }
      pmax(0, pmin(1, sc))
    },
    empirical_percentile = {
      pctiles <- stats::quantile(
        valid, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE
      )
      sc <- numeric(length(values))
      for (i in seq_along(values)) {
        v <- values[i]
        if (is.na(v) || !is.finite(v)) { sc[i] <- NA_real_; next }
        if (lower) {
          if      (v <= pctiles[1L]) sc[i] <- 0.0
          else if (v <= pctiles[2L]) sc[i] <- 0.2
          else if (v <= pctiles[3L]) sc[i] <- 0.4
          else if (v <= pctiles[4L]) sc[i] <- 0.6
          else if (v <= pctiles[5L]) sc[i] <- 0.8
          else                       sc[i] <- 1.0
        } else {
          if      (v >= pctiles[5L]) sc[i] <- 0.0
          else if (v >= pctiles[4L]) sc[i] <- 0.2
          else if (v >= pctiles[3L]) sc[i] <- 0.4
          else if (v >= pctiles[2L]) sc[i] <- 0.6
          else if (v >= pctiles[1L]) sc[i] <- 0.8
          else                       sc[i] <- 1.0
        }
      }
      sc
    },
    percentile = {
      sc <- stats::ecdf(valid)(values)
      if (!lower) sc <- 1 - sc
      sc
    },
    minmax = {
      mn <- min(valid); mx <- max(valid); rng <- mx - mn
      sc <- if (rng == 0) rep(0.5, length(values)) else (values - mn) / rng
      if (lower) sc else 1 - sc
    },
    zscore = {
      mu <- mean(valid); s <- stats::sd(valid)
      sc <- if (s == 0) rep(0, length(values)) else (values - mu) / s
      if (lower) sc else -sc
    },
    divide_max = {
      mx <- max(abs(valid))
      sc <- if (mx == 0) rep(0, length(values)) else values / mx
      if (lower) sc else 1 - sc
    }
  )
  scaled
}

#' @noRd
resilience_smooth_ <- function(values, window_size, method) {
  n <- length(values)
  if (window_size <= 1L || n < window_size) return(values)
  smoothed <- numeric(n)
  half <- floor(window_size / 2)
  smooth_fn <- if (method == "median") stats::median else mean
  for (i in seq_len(n)) {
    s <- max(1L, i - half)
    e <- min(n, i + half)
    smoothed[i] <- smooth_fn(values[s:e], na.rm = TRUE)
  }
  smoothed
}

# --------------------------------------------------------------------------
# S3 methods
# --------------------------------------------------------------------------

#' Plot Resilience Results
#'
#' Visualizes resilience metric trajectories alongside the original time
#' series, either as a combined ribbon plot or as faceted line plots.
#'
#' @details
#' Two complementary visualizations are available:
#'
#' **`type = "ribbons"` (default).** Produces a single-panel plot combining
#' the original time series (black line, upper portion) with a set of
#' gradient-colored horizontal bands below it. Each band corresponds to one
#' resilience metric (plus the composite score). Color at each time point
#' encodes the directional score from [classify_resilience()], ranging from
#' green ("Excellent") through yellow ("Fair") to red ("Troubled"). This
#' layout answers the question *"when and where does the system shift between
#' resilience states?"* at a glance, because all metrics share the same time
#' axis. Requires classified data; call [classify_resilience()] first.
#'
#' **`type = "lines"`.** Produces a faceted plot with one panel per raw
#' metric, each showing the metric's rolling-window value over time. No
#' classification is needed. This view answers the question *"what is the
#' raw trajectory of each metric?"* and is useful for diagnosing which
#' metrics drive changes visible in the ribbon plot. Panels use free y-axes
#' so that metrics on different scales are individually legible.
#'
#' @export
#' @param x \[`resilience`\]\cr
#'   An object of class `resilience`, optionally augmented by
#'   [classify_resilience()] (required for `type = "ribbons"`).
#' @param type \[`character(1)`: `"ribbons"`\]\cr
#'   Plot type. The available options are:
#'   * `"ribbons"`: Time series with gradient-colored risk ribbon below.
#'     Requires classified data from [classify_resilience()].
#'   * `"lines"`: Faceted line plots of raw metric values.
#' @param ... Additional arguments (currently unused).
#' @return A [ggplot2::ggplot()] object that can be further customized with
#'   standard ggplot2 functions (e.g., `+ theme()`, `+ labs()`).
#'
#' @seealso [resilience()] for computing the metrics;
#'   [classify_resilience()] for directional normalization required by the
#'   ribbon plot.
#' @family resilience
#' @concept time series
#' @concept resilience
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' res <- resilience(x, window = 50L)
#'
#' # Faceted raw metric trajectories (no classification needed)
#' plot(res, type = "lines")
#'
#' # Ribbon plot (requires classification first)
#' cls <- classify_resilience(res)
#' plot(cls, type = "ribbons")
#' }
plot.resilience <- function(x, type = "ribbons", ...) {
  check_missing(x)
  check_class(x, "resilience")
  type <- check_match(type, c("ribbons", "lines"))
  if (type == "ribbons") {
    plot_resilience_ribbons_(x, ...)
  } else {
    plot_resilience_lines_(x, ...)
  }
}

#' Print a Resilience Object
#'
#' @describeIn resilience Print method for `"resilience"` objects. Dispatches
#'   to the tibble print method, displaying the time series data and computed
#'   resilience metrics.
#'
#' @export
#' @param x \[`resilience`\]\cr
#'   A `resilience` object.
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
print.resilience <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

# --------------------------------------------------------------------------
# Internal plot helpers (@noRd)
# --------------------------------------------------------------------------

#' @noRd
resilience_display_names <- function() {
  c(
    vsi            = "Stability",
    arch_lm        = "Volatility",
    cv             = "Variability",
    recovery_time  = "Recovery time",
    recovery_slope = "Recovery rate",
    sample_entropy = "Complexity",
    dfa_alpha      = "Memory",
    ac_ratio       = "Correlation",
    composite      = "OVERALL"
  )
}

#' @noRd
plot_resilience_ribbons_ <- function(x, show_metrics = NULL, ...) {
  has_scores <- any(grepl("_score$", names(x)))
  if (!has_scores) {
    stop_(
      "Ribbon plot requires classified data.
      Run {.fn classify_resilience} first."
    )
  }
  display <- resilience_display_names()
  all_metrics <- attr(x, "metrics")

  # Default: exclude recovery_time and recovery_slope (matches original)
  if (is.null(show_metrics)) {
    plot_metrics <- all_metrics[!all_metrics %in%
                                  c("recovery_time", "recovery_slope")]
  } else if (length(show_metrics) == 1L && show_metrics == "all") {
    plot_metrics <- all_metrics
  } else {
    plot_metrics <- show_metrics
  }

  # Order: composite first, then individual metrics
  ordered_metrics <- c("composite", plot_metrics)

  time_vals <- x[["time"]]
  ts_vals   <- x[["value"]]

  ts_range <- diff(range(ts_vals, na.rm = TRUE))
  ts_min   <- min(ts_vals, na.rm = TRUE)
  ts_max   <- max(ts_vals, na.rm = TRUE)

  # Compress time series into upper 60% of its range

  ts_compressed_range <- ts_range * 0.6
  ts_new_max <- ts_max
  ts_new_min <- ts_new_max - ts_compressed_range

  # Ribbon sizing
  ribbon_height  <- ts_range * 0.08
  ribbon_spacing <- ribbon_height * 1.3
  ribbon_start_y <- ts_new_min - ts_range * 0.15

  ts_df <- data.frame(time = time_vals, value = ts_vals)

  p <- ggplot2::ggplot(ts_df, ggplot2::aes(
    x = !!rlang::sym("time"), y = !!rlang::sym("value")
  )) +
    ggplot2::geom_line(color = "black", linewidth = 0.5) +
    ggplot2::labs(
      title = "System Resilience Analysis",
      x     = "Time",
      y     = "Value"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(size = 10),
      plot.title       = ggplot2::element_text(size = 14, face = "bold"),
      legend.text      = ggplot2::element_text(size = 8),
      axis.title       = ggplot2::element_text(color = "black", face = "bold"),
      axis.text        = ggplot2::element_text(color = "black")
    )

  time_diff <- if (length(time_vals) > 1L) {
    stats::median(diff(time_vals), na.rm = TRUE)
  } else 1
  t_min <- min(time_vals, na.rm = TRUE)
  t_max <- max(time_vals, na.rm = TRUE)

  for (i in seq_along(ordered_metrics)) {
    metric <- ordered_metrics[i]
    score_col <- if (metric == "composite") "composite_score" else
      paste0(metric, "_score")
    if (!score_col %in% names(x)) next

    ribbon_y_bottom <- ribbon_start_y - (i * ribbon_spacing)
    ribbon_y_center <- ribbon_y_bottom + ribbon_height / 2

    scores <- x[[score_col]]
    scores[is.na(scores)] <- -1

    tile_df <- data.frame(
      x        = time_vals,
      y        = ribbon_y_center,
      severity = scores
    )

    p <- p + ggplot2::geom_tile(
      data   = tile_df,
      ggplot2::aes(x = !!rlang::sym("x"),
                   y = !!rlang::sym("y"),
                   fill = !!rlang::sym("severity")),
      width  = time_diff,
      height = ribbon_height,
      alpha  = 0.8
    )

    # Label on the right side
    label_color <- if (metric == "composite") "darkblue" else "black"
    label_size  <- if (metric == "composite") 3.5 else 3
    disp_name   <- if (metric == "composite") "OVERALL" else
      display[metric]

    p <- p + ggplot2::annotate(
      "text",
      x        = t_max + (t_max - t_min) / 70,
      y        = ribbon_y_center,
      label    = disp_name,
      hjust    = 0,
      size     = label_size,
      fontface = "bold",
      color    = label_color
    )
  }

  # Extend x-axis to make room for labels; set y limits
  n_metrics <- length(ordered_metrics)
  y_bottom  <- ribbon_start_y - (n_metrics * ribbon_spacing) - ribbon_height
  p <- p +
    ggplot2::coord_cartesian(
      xlim = c(t_min, t_max + (t_max - t_min) / 6),
      ylim = c(y_bottom, ts_new_max + ts_range * 0.1),
      clip = "off"
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c("#cccccc", "#2d7d32", "#66bb6a", "#ffeb99",
                 "#ffcc66", "#ff9999", "#cc0000"),
      values = (c(-1, 0, 0.2, 0.4, 0.6, 0.8, 1.0) + 1) / 2,
      name   = "Resilience Level",
      breaks = c(-1, 0, 0.2, 0.4, 0.6, 0.8, 1.0),
      labels = c("Insufficient", "Excellent", "Solid", "Fair",
                  "Vulnerable", "Failing", "Troubled"),
      na.value = "#cccccc",
      guide = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust    = 0.5,
        barwidth       = 24,
        barheight      = 1
      )
    )
  p
}

#' @noRd
plot_resilience_lines_ <- function(x, ...) {
  metrics <- attr(x, "metrics")
  display <- resilience_display_names()
  long <- tidyr::pivot_longer(
    dplyr::select(x, !!rlang::sym("time"), tidyselect::all_of(metrics)),
    cols = tidyselect::all_of(metrics),
    names_to = "metric",
    values_to = "score"
  )
  long <- dplyr::mutate(
    long,
    metric_label = factor(
      display[!!rlang::sym("metric")],
      levels = display[metrics]
    )
  )
  ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = !!rlang::sym("time"),
      y = !!rlang::sym("score"),
      color = !!rlang::sym("metric_label")
    )
  ) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::facet_wrap(
      ggplot2::vars(!!rlang::sym("metric_label")),
      scales = "free_y",
      ncol = 2L
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Time", y = "Value", color = "Metric") +
    ggplot2::theme(legend.position = "none")
}
