# ============================================================================
# hurst.R --Hurst exponent, EWS detection, and visualization
# ============================================================================

# --------------------------------------------------------------------------
# Exported: hurst()
# --------------------------------------------------------------------------

#' Hurst Exponent Analysis for Time Series Data
#'
#' Computes the Hurst exponent using Detrended Fluctuation Analysis (DFA),
#' Rescaled Range (R/S) analysis, or Multifractal DFA (MFDFA). Supports
#' rolling-window computation with state classification for tracking
#' temporal dynamics of long-range dependence.
#'
#' @details
#' The Hurst exponent \eqn{H} quantifies long-range dependence in a time
#' series. Three estimation methods are available:
#'
#' **DFA (Detrended Fluctuation Analysis).**
#' The series is first integrated into a cumulative profile
#' \eqn{Y(k) = \sum_{i=1}^{k}(x_i - \bar{x})}. The profile is divided
#' into non-overlapping segments of length \eqn{s}, each segment is
#' locally detrended by a least-squares linear fit, and the root-mean-square
#' (RMS) fluctuation \eqn{F(s)} is computed. Fitting
#' \eqn{\log F(s)} vs. \eqn{\log s} yields the Hurst exponent as the slope.
#' Both forward and backward segment directions are used to reduce
#' boundary effects. DFA is robust to non-stationarity and polynomial
#' trends in the data.
#'
#' **R/S (Rescaled Range).**
#' The classic Hurst method partitions the series into segments of length
#' \eqn{s}, computes the range of the cumulative deviation from the segment
#' mean, and divides by the segment standard deviation to obtain the
#' rescaled range \eqn{R/S}. The Hurst exponent is the slope of
#' \eqn{\log(R/S)} vs. \eqn{\log s}.
#'
#' **MFDFA (Multifractal DFA).**
#' Generalizes DFA by computing fluctuation functions for a range of moment
#' orders \eqn{q}, yielding a generalized Hurst exponent \eqn{h(q)} for
#' each \eqn{q}. Positive \eqn{q} emphasizes large fluctuations; negative
#' \eqn{q} emphasizes small fluctuations. The multifractal spectrum width
#' \eqn{\Delta h = \max(h(q)) - \min(h(q))} measures the degree of
#' multifractality. The overall Hurst exponent is taken at \eqn{q = 2},
#' which corresponds to standard DFA.
#'
#' **State classification.**
#' Rolling Hurst values are classified into five regimes based on standard
#' thresholds:
#' \itemize{
#'   \item \eqn{H < 0.4}: `strong_antipersistent` -- strong mean-reverting
#'     dynamics.
#'   \item \eqn{0.4 \le H < 0.5}: `antipersistent` -- weak mean-reversion.
#'   \item \eqn{0.5 \le H < 0.6}: `random_walk` -- uncorrelated or weakly
#'     correlated increments.
#'   \item \eqn{0.6 \le H < 0.7}: `persistent` -- positive long-range
#'     correlations.
#'   \item \eqn{H \ge 0.7}: `strong_persistent` -- strong trending
#'     behavior.
#' }
#'
#' **Rolling window alignment.**
#' Windows are centered on each time point. Edge observations that fall
#' outside the rolling coverage are filled by a three-step interpolation
#' strategy: (1) forward-fill from the first valid estimate to the series
#' start, (2) backward-fill from the last valid estimate to the series end,
#' and (3) linear interpolation for any internal gaps.
#'
#' @export
#' @param data \[`ts`, `numeric()`\]\cr
#'   Univariate time series data.
#' @param method \[`character(1)`: `"dfa"`\]\cr
#'   Hurst estimation method. The available options are:
#'   * `"dfa"`: Detrended Fluctuation Analysis -- robust to
#'     non-stationarity and trends.
#'   * `"rs"`: Rescaled Range analysis -- classic method based on
#'     range-over-standard-deviation scaling.
#'   * `"mfdfa"`: Multifractal DFA -- extends DFA to characterize the
#'     full multifractal spectrum via generalized Hurst exponents.
#' @param window \[`integer(1)`: `50L`\]\cr
#'   Rolling window size for local Hurst estimation when `states = TRUE`.
#' @param step \[`integer(1)`: `1L`\]\cr
#'   Step size between consecutive windows.
#' @param scaling \[`character(1)`: `"none"`\]\cr
#'   Preprocessing applied to the data. The available options are:
#'   `"none"`, `"center"`, `"standardize"`, `"minmax"`, and `"iqr"`.
#' @param min_scale \[`integer(1)`: `4L`\]\cr
#'   Minimum box size for the DFA/R/S log-log regression.
#' @param max_scale \[`integer(1)`: `NULL`\]\cr
#'   Maximum box size. If `NULL`, defaults to `floor(n / 4)` for DFA/MFDFA
#'   or `floor(n / 2)` for R/S.
#' @param n_scales \[`integer(1)`: `10L`\]\cr
#'   Number of scales (box sizes) to evaluate between `min_scale` and
#'   `max_scale` on a logarithmic grid.
#' @param q \[`numeric()`: `seq(-5, 5, 0.5)`\]\cr
#'   Vector of moment orders for MFDFA. Only used when `method = "mfdfa"`.
#' @param states \[`logical(1)`: `TRUE`\]\cr
#'   If `TRUE`, compute rolling Hurst values with state classification.
#'   If `FALSE`, compute a single global Hurst estimate.
#' @return
#'   If `states = TRUE`, an object of class `"hurst"` (inheriting from
#'   [tibble::tibble()]) with the following columns:
#'   \itemize{
#'     \item `time` -- integer or Date time index.
#'     \item `value` -- original series values.
#'     \item `hurst` -- rolling Hurst exponent estimate.
#'     \item `r_squared` -- \eqn{R^2} of the log-log regression at each
#'       window.
#'     \item `state` -- character state label (one of
#'       `"strong_antipersistent"`, `"antipersistent"`, `"random_walk"`,
#'       `"persistent"`, `"strong_persistent"`).
#'     \item `transition` -- numeric; nonzero where the state changes from
#'       the previous time point.
#'     \item `mf_width` -- (MFDFA only) multifractal spectrum width
#'       \eqn{\Delta h}.
#'     \item `mf_category` -- (MFDFA only) character classification:
#'       `"monofractal"`, `"weak_multifractal"`,
#'       `"moderate_multifractal"`, or `"strong_multifractal"`.
#'   }
#'   Attributes: `window`, `step`, `method`, `scaling`.
#'
#'   If `states = FALSE`, an object of class `"hurst_global"` (a list)
#'   with the following elements:
#'   \itemize{
#'     \item `hurst` -- single global Hurst exponent.
#'     \item `r_squared` -- \eqn{R^2} of the log-log fit.
#'     \item `method` -- character string naming the method used.
#'     \item `n` -- integer series length.
#'     \item `scales` -- numeric vector of scale sizes used.
#'     \item `fluctuations` -- (DFA) numeric vector of \eqn{F(s)} values.
#'     \item `rs_values` -- (R/S) numeric vector of rescaled-range values.
#'     \item `hq` -- (MFDFA) numeric vector of generalized Hurst
#'       exponents \eqn{h(q)}.
#'     \item `tauq` -- (MFDFA) numeric vector of mass exponents
#'       \eqn{\tau(q)}.
#'     \item `mf_width` -- (MFDFA) multifractal spectrum width.
#'   }
#'
#' @references
#' Peng, C.-K., Buldyrev, S.V., Havlin, S., Simons, M., Stanley, H.E.,
#' & Goldberger, A.L. (1994). Mosaic organization of DNA nucleotides.
#' \emph{Physical Review E}, 49(2), 1685--1689.
#' \doi{10.1103/PhysRevE.49.1685}
#'
#' Hurst, H.E. (1951). Long-term storage capacity of reservoirs.
#' \emph{Transactions of the American Society of Civil Engineers}, 116,
#' 770--799.
#'
#' Kantelhardt, J.W., Zschiegner, S.A., Koscielny-Bunde, E.,
#' Havlin, S., Bunde, A., & Stanley, H.E. (2002). Multifractal
#' detrended fluctuation analysis of nonstationary time series.
#' \emph{Physica A}, 316(1--4), 87--114.
#' \doi{10.1016/S0378-4371(02)01383-3}
#'
#' @seealso [detect_hurst_warnings()] for early warning signal detection
#'   from rolling Hurst output; [plot.hurst()] for visualization.
#' @family hurst
#' @concept time series
#' @concept hurst exponent
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' h <- hurst(x, window = 50L)
#' plot(h)
#' }
hurst <- function(data,
                  method = "dfa",
                  window = 50L,
                  step = 1L,
                  scaling = "none",
                  min_scale = 4L,
                  max_scale = NULL,
                  n_scales = 10L,
                  q = seq(-5, 5, 0.5),
                  states = TRUE) {
  check_missing(data)
  data <- prepare_timeseries_data(data)
  values <- data$values
  time <- data$time
  n <- length(values)
  method <- check_match(method, c("dfa", "rs", "mfdfa"))
  scaling <- check_match(
    scaling,
    c("none", "center", "standardize", "minmax", "iqr")
  )
  min_win <- 2L * min_scale
  check_range(window, type = "integer", min = min_win, max = n)
  check_range(step, type = "integer", min = 1L, max = window)
  check_flag(states)
  values <- hurst_scale(values, scaling)
  # Rolling n_scales reduced to min(n_scales, 5) to match Resilience4
  rolling_n_scales <- min(n_scales, 5L)
  hurst_fun <- switch(
    method,
    dfa = function(x, ns) {
      hurst_dfa(x, min_scale, min(max_scale %||% Inf,
                                   floor(length(x) / 4)), ns)
    },
    rs = function(x, ns) {
      hurst_rs(x, min_scale, min(max_scale %||% Inf,
                                  floor(length(x) / 4)), ns)
    },
    mfdfa = function(x, ns) {
      hurst_mfdfa(x, q, min_scale, min(max_scale %||% Inf,
                                        floor(length(x) / 4)), ns)
    }
  )
  if (!states) {
    result <- hurst_fun(values, n_scales)
    out <- list(
      hurst = result$hurst,
      r_squared = result$r_squared,
      method = method,
      n = n
    )
    if (method == "dfa") {
      out$scales <- result$scales
      out$fluctuations <- result$fluctuations
    } else if (method == "rs") {
      out$scales <- result$scales
      out$rs_values <- result$rs_values
    } else {
      out$scales <- result$scales %||% integer(0)
      out$hq <- result$hq
      out$tauq <- result$tauq
      out$mf_width <- result$mf_width
    }
    return(structure(out, class = c("hurst_global", "list")))
  }
  # Rolling centered window --matches Resilience4 exactly
  half <- floor(window / 2)
  hurst_state <- rep(NA_real_, n)
  r2_state <- rep(NA_real_, n)
  mf_state <- if (method == "mfdfa") rep(NA_real_, n) else NULL
  for (i in seq(half + 1L, n - half, by = step)) {
    start_idx <- max(1L, i - half)
    end_idx <- min(n, i + half - 1L)
    # Maintain consistent window size at boundaries
    actual_ws <- end_idx - start_idx + 1L
    if (actual_ws < window && start_idx > 1L) {
      start_idx <- max(1L, end_idx - window + 1L)
    } else if (actual_ws < window && end_idx < n) {
      end_idx <- min(n, start_idx + window - 1L)
    }
    window_data <- values[start_idx:end_idx]
    # Skip if too many NAs
    if (sum(!is.na(window_data)) < 0.8 * length(window_data)) next
    result <- try_(hurst_fun(window_data, rolling_n_scales))
    if (inherits(result, "try-error")) result <- NULL
    if (!is.null(result)) {
      hurst_state[i] <- result$hurst
      r2_state[i] <- result$r_squared
      if (method == "mfdfa") mf_state[i] <- result$mf_width
    }
  }
  # Three-step interpolation matching Resilience4
  first_valid <- which(!is.na(hurst_state))[1L]
  last_valid <- rev(which(!is.na(hurst_state)))[1L]
  if (!is.na(first_valid) && !is.na(last_valid)) {
    # Forward-fill beginning
    if (first_valid > 1L) {
      hurst_state[seq_len(first_valid - 1L)] <- hurst_state[first_valid]
      r2_state[seq_len(first_valid - 1L)] <- r2_state[first_valid]
      if (method == "mfdfa")
        mf_state[seq_len(first_valid - 1L)] <- mf_state[first_valid]
    }
    # Backward-fill end
    if (last_valid < n) {
      hurst_state[(last_valid + 1L):n] <- hurst_state[last_valid]
      r2_state[(last_valid + 1L):n] <- r2_state[last_valid]
      if (method == "mfdfa")
        mf_state[(last_valid + 1L):n] <- mf_state[last_valid]
    }
    # Linear interpolation for internal gaps
    na_pos <- which(is.na(hurst_state))
    internal_na <- na_pos[na_pos > first_valid & na_pos < last_valid]
    if (length(internal_na) > 0L) {
      valid_pos <- which(!is.na(hurst_state))
      hurst_state <- stats::approx(
        valid_pos, hurst_state[valid_pos], xout = seq_len(n), rule = 2
      )$y
      r2_state <- stats::approx(
        valid_pos, r2_state[valid_pos], xout = seq_len(n), rule = 2
      )$y
      if (method == "mfdfa") {
        mf_state <- stats::approx(
          valid_pos, mf_state[valid_pos], xout = seq_len(n), rule = 2
        )$y
      }
    }
  }
  state <- hurst_classify(hurst_state)
  transition <- c(0, diff(as.numeric(factor(state))))
  out <- tibble::tibble(
    time = time,
    value = data$values,
    hurst = hurst_state,
    r_squared = r2_state,
    state = state,
    transition = transition
  )
  if (method == "mfdfa") {
    out$mf_width <- mf_state
    out$mf_category <- hurst_mf_classify(mf_state)
  }
  structure(
    out,
    window = window,
    step = step,
    method = method,
    scaling = scaling,
    class = c("hurst", "tbl_df", "tbl", "data.frame")
  )
}

# --------------------------------------------------------------------------
# Internal computation helpers (@noRd)
# --------------------------------------------------------------------------

#' @noRd
hurst_scale <- function(x, scaling) {
  switch(
    scaling,
    none = x,
    center = x - mean(x, na.rm = TRUE),
    standardize = {
      s <- stats::sd(x, na.rm = TRUE)
      if (s == 0) x - mean(x, na.rm = TRUE)
      else (x - mean(x, na.rm = TRUE)) / s
    },
    minmax = {
      mn <- min(x, na.rm = TRUE)
      mx <- max(x, na.rm = TRUE)
      rng <- mx - mn
      if (rng == 0) rep(0.5, length(x)) else (x - mn) / rng
    },
    iqr = {
      q <- stats::quantile(x, c(0.25, 0.75), na.rm = TRUE)
      iqr_val <- q[2L] - q[1L]
      if (iqr_val == 0) x - stats::median(x, na.rm = TRUE)
      else (x - stats::median(x, na.rm = TRUE)) / iqr_val
    }
  )
}

#' @noRd
hurst_dfa <- function(x, min_scale, max_scale, n_scales) {
  # Matches Resilience4 .calculate_hurst_dfa() exactly
  if (any(is.na(x))) x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2L * min_scale) {
    return(list(hurst = NA_real_, r_squared = NA_real_,
                scales = integer(0), fluctuations = numeric(0)))
  }
  if (min_scale < 4L) min_scale <- 4L
  if (is.null(max_scale)) {
    max_scale <- floor(n / 4)
  } else if (max_scale > n / 2) {
    max_scale <- floor(n / 4)
  }
  if (max_scale <= min_scale) {
    return(list(hurst = NA_real_, r_squared = NA_real_,
                scales = integer(0), fluctuations = numeric(0)))
  }
  # Logarithmically spaced scales with round() (not floor)
  log_sc <- seq(log(min_scale), log(max_scale), length.out = n_scales)
  scales <- unique(round(exp(log_sc)))
  scales <- scales[scales >= min_scale & scales <= max_scale]
  n_sc <- length(scales)
  # Cumulative profile
  ts_mean <- mean(x)
  profile <- cumsum(x - ts_mean)
  fluct <- numeric(n_sc)
  for (si in seq_len(n_sc)) {
    s <- scales[si]
    n_seg <- floor(n / s)
    if (n_seg < 2L) {
      fluct[si] <- NA_real_
      next
    }
    seg_fluct <- numeric(2L * n_seg)
    x_vals <- seq_len(s)
    x_mean <- mean(x_vals)
    xx_var <- sum((x_vals - x_mean)^2)
    # Forward direction
    for (j in seq_len(n_seg)) {
      start_idx <- (j - 1L) * s + 1L
      end_idx <- j * s
      segment <- profile[start_idx:end_idx]
      if (xx_var > 0) {
        y_mean <- mean(segment)
        xy_cov <- sum((x_vals - x_mean) * (segment - y_mean))
        slope <- xy_cov / xx_var
        intercept <- y_mean - slope * x_mean
        trend <- intercept + slope * x_vals
        detrended <- segment - trend
        seg_fluct[j] <- sqrt(mean(detrended^2))
      } else {
        seg_fluct[j] <- 0
      }
    }
    # Backward direction
    for (j in seq_len(n_seg)) {
      start_idx <- n - j * s + 1L
      end_idx <- n - (j - 1L) * s
      segment <- profile[start_idx:end_idx]
      if (xx_var > 0) {
        y_mean <- mean(segment)
        xy_cov <- sum((x_vals - x_mean) * (segment - y_mean))
        slope <- xy_cov / xx_var
        intercept <- y_mean - slope * x_mean
        trend <- intercept + slope * x_vals
        detrended <- segment - trend
        seg_fluct[n_seg + j] <- sqrt(mean(detrended^2))
      } else {
        seg_fluct[n_seg + j] <- 0
      }
    }
    fluct[si] <- sqrt(mean(seg_fluct^2))
  }
  valid <- !is.na(fluct) & fluct > 0
  if (sum(valid) < 3L) {
    return(list(hurst = NA_real_, r_squared = NA_real_,
                scales = scales, fluctuations = fluct))
  }
  log_s <- log(scales[valid])
  log_f <- log(fluct[valid])
  x_m <- mean(log_s)
  y_m <- mean(log_f)
  xy_cov <- sum((log_s - x_m) * (log_f - y_m))
  xx_var <- sum((log_s - x_m)^2)
  slope <- xy_cov / xx_var
  intercept <- y_m - slope * x_m
  y_pred <- intercept + slope * log_s
  ss_tot <- sum((log_f - y_m)^2)
  ss_res <- sum((log_f - y_pred)^2)
  r_sq <- 1 - (ss_res / ss_tot)
  list(
    hurst = slope,
    r_squared = r_sq,
    scales = scales[valid],
    fluctuations = fluct[valid]
  )
}

#' @noRd
hurst_rs <- function(x, min_scale, max_scale, n_scales) {
  n <- length(x)
  if (is.null(max_scale)) max_scale <- floor(n / 2)
  max_scale <- min(max_scale, floor(n / 2))
  if (max_scale < min_scale) {
    return(list(hurst = NA_real_, r_squared = NA_real_,
                scales = integer(0), rs_values = numeric(0)))
  }
  scales <- unique(floor(
    exp(seq(log(min_scale), log(max_scale), length.out = n_scales))
  ))
  rs_vals <- numeric(length(scales))
  for (si in seq_along(scales)) {
    s <- scales[si]
    n_seg <- floor(n / s)
    if (n_seg < 1L) {
      rs_vals[si] <- NA_real_
      next
    }
    rs_seg <- numeric(n_seg)
    for (v in seq_len(n_seg)) {
      idx <- ((v - 1L) * s + 1L):(v * s)
      seg <- x[idx]
      mu <- mean(seg)
      y <- cumsum(seg - mu)
      r <- max(y) - min(y)
      s_val <- stats::sd(seg)
      rs_seg[v] <- ifelse_(s_val == 0, NA_real_, r / s_val)
    }
    rs_vals[si] <- mean(rs_seg, na.rm = TRUE)
  }
  valid <- !is.na(rs_vals) & rs_vals > 0
  if (sum(valid) < 3L) {
    return(list(hurst = NA_real_, r_squared = NA_real_,
                scales = scales, rs_values = rs_vals))
  }
  log_s <- log(scales[valid])
  log_rs <- log(rs_vals[valid])
  fit <- stats::lm.fit(x = cbind(1, log_s), y = log_rs)
  ss_res <- sum(fit$residuals^2)
  ss_tot <- sum((log_rs - mean(log_rs))^2)
  r_sq <- ifelse_(ss_tot == 0, 1, 1 - ss_res / ss_tot)
  list(
    hurst = fit$coefficients[2L],
    r_squared = r_sq,
    scales = scales,
    rs_values = rs_vals
  )
}

#' @noRd
hurst_mfdfa <- function(x, q, min_scale, max_scale, n_scales) {
  n <- length(x)
  if (is.null(max_scale)) max_scale <- floor(n / 4)
  max_scale <- min(max_scale, floor(n / 4))
  if (max_scale < min_scale) {
    return(list(hurst = NA_real_, r_squared = NA_real_,
                hq = numeric(0), tauq = numeric(0), mf_width = NA_real_))
  }
  scales <- unique(floor(
    exp(seq(log(min_scale), log(max_scale), length.out = n_scales))
  ))
  y <- cumsum(x - mean(x, na.rm = TRUE))
  # Compute fluctuation for each scale and q
  fq <- matrix(NA_real_, nrow = length(q), ncol = length(scales))
  for (si in seq_along(scales)) {
    s <- scales[si]
    n_seg <- floor(n / s)
    if (n_seg < 1L) next
    rms_vals <- numeric(n_seg)
    for (v in seq_len(n_seg)) {
      idx <- ((v - 1L) * s + 1L):(v * s)
      seg <- y[idx]
      x_local <- seq_along(seg)
      fit <- stats::lm.fit(x = cbind(1, x_local), y = seg)
      rms_vals[v] <- mean(fit$residuals^2)
    }
    rms_vals <- rms_vals[rms_vals > 0]
    if (length(rms_vals) == 0L) next
    for (qi in seq_along(q)) {
      qv <- q[qi]
      if (abs(qv) < .Machine$double.eps) {
        fq[qi, si] <- exp(0.5 * mean(log(rms_vals)))
      } else {
        fq[qi, si] <- (mean(rms_vals^(qv / 2)))^(1 / qv)
      }
    }
  }
  # Compute generalized Hurst exponent h(q) for each q
  hq <- numeric(length(q))
  for (qi in seq_along(q)) {
    vals <- fq[qi, ]
    valid <- !is.na(vals) & vals > 0
    if (sum(valid) < 3L) {
      hq[qi] <- NA_real_
      next
    }
    fit <- stats::lm.fit(
      x = cbind(1, log(scales[valid])),
      y = log(vals[valid])
    )
    hq[qi] <- fit$coefficients[2L]
  }
  # Multifractal spectrum
  tauq <- q * hq - 1
  # Overall Hurst at q=2
  q2_idx <- which.min(abs(q - 2))
  h_main <- hq[q2_idx]
  # Multifractal width
  valid_hq <- hq[!is.na(hq)]
  mf_width <- ifelse_(
    length(valid_hq) >= 2L,
    max(valid_hq) - min(valid_hq),
    NA_real_
  )
  # R-squared for q=2
  vals_q2 <- fq[q2_idx, ]
  valid_q2 <- !is.na(vals_q2) & vals_q2 > 0
  r_sq <- NA_real_
  if (sum(valid_q2) >= 3L) {
    log_s <- log(scales[valid_q2])
    log_f <- log(vals_q2[valid_q2])
    fit <- stats::lm.fit(x = cbind(1, log_s), y = log_f)
    ss_res <- sum(fit$residuals^2)
    ss_tot <- sum((log_f - mean(log_f))^2)
    r_sq <- ifelse_(ss_tot == 0, 1, 1 - ss_res / ss_tot)
  }
  list(
    hurst = h_main,
    r_squared = r_sq,
    hq = hq,
    tauq = tauq,
    mf_width = mf_width
  )
}

# --------------------------------------------------------------------------
# Internal state classification helpers (@noRd)
# --------------------------------------------------------------------------

#' @noRd
hurst_classify <- function(h) {
  labels <- c(
    "strong_antipersistent", "antipersistent", "random_walk",
    "persistent", "strong_persistent"
  )
  breaks <- c(-Inf, 0.4, 0.5, 0.6, 0.7, Inf)
  as.character(cut(h, breaks = breaks, labels = labels, right = FALSE))
}

#' @noRd
hurst_mf_classify <- function(width) {
  labels <- c("monofractal", "weak_multifractal",
              "moderate_multifractal", "strong_multifractal")
  breaks <- c(-Inf, 0.1, 0.3, 0.5, Inf)
  as.character(cut(width, breaks = breaks, labels = labels, right = FALSE))
}

#' @noRd
hurst_state_colors <- function() {
  c(
    strong_antipersistent = "#FF6B6B",
    antipersistent        = "#FFA06B",
    random_walk           = "#FFD93D",
    persistent            = "#6BCF7F",
    strong_persistent     = "#4ECDC4"
  )
}

# --------------------------------------------------------------------------
# Exported: detect_hurst_warnings()
# --------------------------------------------------------------------------

#' Detect Early Warning Signals from Hurst Analysis
#'
#' Identifies early warning signals of critical transitions by analyzing
#' temporal patterns in the Hurst exponent, including trends, volatility,
#' flickering, and spectral shifts. Produces a graded warning level for
#' each time point.
#'
#' @details
#' The function computes 10 binary (0/1) indicators at each time point:
#' \enumerate{
#'   \item `extreme_low` -- Hurst falls below the 5th data-driven
#'     percentile.
#'   \item `extreme_high` -- Hurst exceeds the 95th data-driven
#'     percentile.
#'   \item `trend_up` -- significant upward trend in the Hurst exponent
#'     within a rolling window (slope exceeds a dynamic threshold and
#'     \eqn{t}-statistic > 2).
#'   \item `trend_down` -- significant downward trend (symmetric to
#'     `trend_up`).
#'   \item `high_volatility` -- rolling RMSD of Hurst exceeds the 90th
#'     percentile of all rolling volatility values.
#'   \item `flickering` -- rapid state transitions within a rolling
#'     window (> 20\% transitions or \eqn{\ge 3} distinct states).
#'   \item `variance_ratio` -- deviation of the variance ratio at lags
#'     2 and 4 from a random-walk baseline exceeds the 90th percentile,
#'     computed over 60-point windows.
#'   \item `spectral_shift` -- dominant spectral period in the recent
#'     50-point window exceeds 1.5 times the dominant period of the
#'     preceding 50-point window, indicating critical slowing down.
#'   \item `autocorr_increase` -- mean lag-1 to lag-5 autocorrelation
#'     in the recent 50-point window exceeds the preceding window by
#'     more than 0.2, consistent with rising memory before a transition.
#'   \item `state_persistence` -- the current run of identical states
#'     exceeds the 90th percentile of all run lengths, indicating the
#'     system is locked into a regime.
#' }
#'
#' **Weighted scoring.**
#' Each indicator is multiplied by a fixed weight (extreme and flickering
#' indicators: 2.0; trend and volatility indicators: 1.5; autocorrelation
#' increase: 1.5; variance ratio, spectral shift, state persistence: 1.0)
#' and the products are summed to produce a raw warning score per time
#' point.
#'
#' **Proximity enhancement.**
#' Raw scores are boosted when neighbouring time points (within a
#' \eqn{\pm 10}-point window) also carry nonzero scores. The boost factor
#' is \eqn{1 + p}, where \eqn{p} is the proportion of neighbours with
#' \eqn{\mathrm{score} > 0}. This smooths isolated false positives while
#' amplifying clustered warnings.
#'
#' **Dynamic warning levels.**
#' Enhanced scores are assigned to five levels (0--4) using quantile
#' breakpoints at the 40th, 60th, 80th, and 95th percentiles of the
#' nonzero score distribution:
#' \itemize{
#'   \item Level 0 (`"none"`) -- score is zero.
#'   \item Level 1 (`"low"`) -- score \eqn{\le} 40th percentile.
#'   \item Level 2 (`"moderate"`) -- score \eqn{\le} 60th percentile.
#'   \item Level 3 (`"high"`) -- score \eqn{\le} 80th percentile.
#'   \item Level 4 (`"critical"`) -- score > 95th percentile.
#' }
#' Because thresholds are data-driven, the system adapts to each series
#' without requiring tuning.
#'
#' @export
#' @param data \[`hurst`\]\cr
#'   An object of class `hurst` as returned by [hurst()] with
#'   `states = TRUE`.
#' @param trend_window \[`integer(1)`: `30L`\]\cr
#'   Window size for detecting trends in the Hurst exponent.
#' @param volatility_window \[`integer(1)`: `30L`\]\cr
#'   Window size for measuring volatility of the Hurst exponent.
#' @param flicker_window \[`integer(1)`: `20L`\]\cr
#'   Window size for detecting flickering (rapid state transitions).
#' @return An object of class `"hurst_ews"` (inheriting from
#'   [tibble::tibble()]) containing the original `hurst` columns plus:
#'   \itemize{
#'     \item `extreme_low` -- binary indicator (0/1).
#'     \item `extreme_high` -- binary indicator (0/1).
#'     \item `trend_up` -- binary indicator (0/1).
#'     \item `trend_down` -- binary indicator (0/1).
#'     \item `high_volatility` -- binary indicator (0/1).
#'     \item `flickering` -- binary indicator (0/1).
#'     \item `variance_ratio` -- binary indicator (0/1).
#'     \item `spectral_shift` -- binary indicator (0/1).
#'     \item `autocorr_increase` -- binary indicator (0/1).
#'     \item `state_persistence` -- binary indicator (0/1).
#'     \item `warning_score` -- numeric score normalized to \eqn{[0, 1]}.
#'     \item `warning_level` -- integer warning level (0--4).
#'     \item `warning_label` -- character label: `"none"`, `"low"`,
#'       `"moderate"`, `"high"`, or `"critical"`.
#'   }
#'   Attribute: `indicator_weights` (named numeric vector of the 10
#'   indicator weights used in scoring).
#'
#' @references
#' Scheffer, M., Bascompte, J., Brock, W.A., Brovkin, V., Carpenter,
#' S.R., Dakos, V., Held, H., van Nes, E.H., Rietkerk, M., &
#' Sugihara, G. (2009). Early-warning signals for critical transitions.
#' \emph{Nature}, 461, 53--59. \doi{10.1038/nature08227}
#'
#' Peng, C.-K., Buldyrev, S.V., Havlin, S., Simons, M., Stanley, H.E.,
#' & Goldberger, A.L. (1994). Mosaic organization of DNA nucleotides.
#' \emph{Physical Review E}, 49(2), 1685--1689.
#' \doi{10.1103/PhysRevE.49.1685}
#'
#' @seealso [hurst()] for computing rolling Hurst exponents;
#'   [plot.hurst_ews()] for three-panel visualization of warning signals.
#' @family hurst
#' @concept time series
#' @concept early warning signals
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' h <- hurst(x, window = 50L)
#' ews <- detect_hurst_warnings(h)
#' summary(ews)
#' }
detect_hurst_warnings <- function(data,
                                  trend_window = 30L,
                                  volatility_window = 30L,
                                  flicker_window = 20L) {
  check_missing(data)
  check_class(data, "hurst")
  n <- nrow(data)
  check_range(trend_window, type = "integer", min = 5L, max = n)
  check_range(volatility_window, type = "integer", min = 5L, max = n)
  check_range(flicker_window, type = "integer", min = 5L, max = n)
  h <- data$hurst
  states <- data$state
  transitions <- data$transition

  # All indicators are binary (0/1), matching the original design
  # ---------------------------------------------------------------

  # --- Indicator 1-2: Extreme values (data-driven percentiles) ---
  thresh_low  <- stats::quantile(h, 0.05, na.rm = TRUE)
  thresh_high <- stats::quantile(h, 0.95, na.rm = TRUE)
  extreme_low  <- as.numeric(h < thresh_low)
  extreme_high <- as.numeric(h > thresh_high)

  # --- Indicator 3-4: Significant trends (slope + t-test) ---
  trend_up   <- rep(0, n)
  trend_down <- rep(0, n)
  trend_thresh <- stats::quantile(
    abs(h - mean(h, na.rm = TRUE)), 0.90, na.rm = TRUE
  ) / trend_window
  if (n >= trend_window) {
    for (i in trend_window:n) {
      idx <- (i - trend_window + 1L):i
      y <- h[idx]
      if (sum(!is.na(y)) >= trend_window * 0.8) {
        x_local <- seq_len(trend_window)
        fit <- stats::lm(y ~ x_local)
        slope <- stats::coef(fit)[2L]
        se <- summary(fit)$coefficients[2L, 2L]
        t_stat <- if (!is.na(se) && se > 0) abs(slope) / se else 0
        if (!is.na(t_stat) && t_stat > 2) {
          if (slope > trend_thresh) trend_up[i] <- 1
          else if (slope < -trend_thresh) trend_down[i] <- 1
        }
      }
    }
  }

  # --- Indicator 5: High volatility (rolling RMSD > 90th pctl) ---
  high_volatility <- rep(0, n)
  if (n >= volatility_window) {
    rolling_vol <- rep(NA_real_, n)
    for (i in volatility_window:n) {
      idx <- (i - volatility_window + 1L):i
      local_mean <- mean(h[idx], na.rm = TRUE)
      rolling_vol[i] <- sqrt(mean((h[idx] - local_mean)^2, na.rm = TRUE))
    }
    vol_thresh <- stats::quantile(rolling_vol, 0.90, na.rm = TRUE)
    high_volatility <- as.numeric(!is.na(rolling_vol) & rolling_vol > vol_thresh)
  }

  # --- Indicator 6: Flickering (many transitions or states in window) ---
  flickering <- rep(0, n)
  if (n >= flicker_window) {
    state_num <- as.numeric(factor(states))
    for (i in flicker_window:n) {
      idx <- (i - flicker_window + 1L):i
      n_unique <- length(unique(state_num[idx]))
      n_trans  <- sum(diff(state_num[idx]) != 0, na.rm = TRUE)
      flickering[i] <- as.numeric(
        n_trans > flicker_window * 0.2 || n_unique >= 3L
      )
    }
  }

  # --- Indicator 7: Variance ratio (top 10% deviations only) ---
  # Computes rolling variance ratio, then flags only the most
  # extreme deviations (> 90th percentile of all local deviations).
  variance_ratio <- rep(0, n)
  vr_window <- 60L
  if (n >= vr_window) {
    vr_raw <- rep(NA_real_, n)
    for (i in vr_window:n) {
      idx <- (i - vr_window + 1L):i
      seg <- h[idx]
      if (sum(!is.na(seg)) >= 50L) {
        var_2 <- stats::var(diff(seg), na.rm = TRUE)
        if (!is.na(var_2) && var_2 > 0) {
          var_4 <- stats::var(diff(seg, lag = 2), na.rm = TRUE) / 2
          var_8 <- stats::var(diff(seg, lag = 4), na.rm = TRUE) / 4
          vr_2 <- abs(var_4 / var_2 - 1)
          vr_4 <- abs(var_8 / var_2 - 1)
          vr_raw[i] <- max(vr_2, vr_4)
        }
      }
    }
    vr_thresh <- stats::quantile(vr_raw, 0.90, na.rm = TRUE)
    variance_ratio <- as.numeric(!is.na(vr_raw) & vr_raw > vr_thresh)
  }

  # --- Indicator 8: Spectral shift (recent vs past dominant period) ---
  # Compares recent 50-point dominant period to past 50-point period.
  # Fires when recent dominant period is substantially longer.
  spectral_shift <- rep(0, n)
  spec_half <- 50L
  if (n >= 2L * spec_half) {
    for (i in (2L * spec_half):n) {
      recent_idx <- (i - spec_half + 1L):i
      past_idx   <- (i - 2L * spec_half + 1L):(i - spec_half)
      sp_recent <- try_(
        stats::spectrum(h[recent_idx], method = "pgram", plot = FALSE)
      )
      sp_past <- try_(
        stats::spectrum(h[past_idx], method = "pgram", plot = FALSE)
      )
      if (!inherits(sp_recent, "try-error") &&
          !inherits(sp_past, "try-error")) {
        period_recent <- 1 / sp_recent$freq[which.max(sp_recent$spec)]
        period_past   <- 1 / sp_past$freq[which.max(sp_past$spec)]
        # Fire when recent period exceeds past period substantially
        spectral_shift[i] <- as.numeric(period_recent > period_past * 1.5)
      }
    }
  }

  # --- Indicator 9: Autocorrelation increase (recent > past + 0.2) ---
  autocorr_increase <- rep(0, n)
  acf_half <- 50L
  if (n >= 2L * acf_half) {
    for (i in (2L * acf_half):n) {
      recent_idx <- (i - acf_half + 1L):i
      past_idx   <- (i - 2L * acf_half + 1L):(i - acf_half)
      recent_acf <- try_(
        stats::acf(h[recent_idx], lag.max = 5L, plot = FALSE)$acf[2L:6L]
      )
      past_acf <- try_(
        stats::acf(h[past_idx], lag.max = 5L, plot = FALSE)$acf[2L:6L]
      )
      if (!inherits(recent_acf, "try-error") &&
          !inherits(past_acf, "try-error")) {
        autocorr_increase[i] <- as.numeric(
          mean(recent_acf, na.rm = TRUE) >
            mean(past_acf, na.rm = TRUE) + 0.2
        )
      }
    }
  }

  # --- Indicator 10: State persistence (runs > 90th percentile) ---
  state_persistence <- rep(0, n)
  state_runs <- rle(as.character(states))
  if (length(state_runs$lengths) > 1L) {
    run_thresh <- stats::quantile(state_runs$lengths, 0.90)
    run_end   <- cumsum(state_runs$lengths)
    run_start <- c(1L, run_end[-length(run_end)] + 1L)
    for (j in seq_along(state_runs$lengths)) {
      if (state_runs$lengths[j] > run_thresh) {
        state_persistence[run_start[j]:run_end[j]] <- 1
      }
    }
  }

  # --- Weighted scoring (matching original) ---
  # Build indicator matrix (binary 0/1)
  indicators <- cbind(
    extreme_low, extreme_high, trend_up, trend_down,
    high_volatility, flickering, variance_ratio,
    spectral_shift, autocorr_increase, state_persistence
  )
  weights <- c(
    extreme_low       = 2.0,
    extreme_high      = 2.0,
    trend_up          = 1.5,
    trend_down        = 1.5,
    high_volatility   = 1.5,
    flickering        = 2.0,
    variance_ratio    = 1.0,
    spectral_shift    = 1.0,
    autocorr_increase = 1.5,
    state_persistence = 1.0
  )
  raw_scores <- as.numeric(indicators %*% weights)

  # Proximity enhancement: boost score where neighbors also have warnings
  proximity_window <- 10L
  enhanced_scores <- raw_scores
  for (i in (proximity_window + 1L):(n - proximity_window)) {
    idx <- (i - proximity_window):(i + proximity_window)
    neighbor_effect <- mean(raw_scores[idx] > 0)
    enhanced_scores[i] <- raw_scores[i] * (1 + neighbor_effect)
  }

  # Dynamic quantile-based warning levels (adapts to data)
  warning_level <- rep(0L, n)
  nonzero <- enhanced_scores[enhanced_scores > 0]
  if (length(nonzero) > 4L) {
    score_q <- stats::quantile(nonzero, c(0.4, 0.6, 0.8, 0.95),
                               na.rm = TRUE)
    warning_level[enhanced_scores > 0 &
                    enhanced_scores <= score_q[1L]] <- 1L
    warning_level[enhanced_scores > score_q[1L] &
                    enhanced_scores <= score_q[2L]] <- 2L
    warning_level[enhanced_scores > score_q[2L] &
                    enhanced_scores <= score_q[3L]] <- 3L
    warning_level[enhanced_scores > score_q[3L]] <- 4L
  }
  # Normalize enhanced_scores to 0-1 for storage
  max_es <- max(enhanced_scores, na.rm = TRUE)
  warning_score <- if (max_es > 0) enhanced_scores / max_es else
    rep(0, n)
  warning_labels <- c("none", "low", "moderate", "high", "critical")
  out <- data
  out$extreme_low <- extreme_low
  out$extreme_high <- extreme_high
  out$trend_up <- trend_up
  out$trend_down <- trend_down
  out$high_volatility <- high_volatility
  out$flickering <- flickering
  out$variance_ratio <- variance_ratio
  out$spectral_shift <- spectral_shift
  out$autocorr_increase <- autocorr_increase
  out$state_persistence <- state_persistence
  out$warning_score <- warning_score
  out$warning_level <- warning_level
  out$warning_label <- warning_labels[warning_level + 1L]
  structure(
    tibble::as_tibble(out),
    indicator_weights = weights,
    class = c("hurst_ews", "tbl_df", "tbl", "data.frame")
  )
}

# --------------------------------------------------------------------------
# S3: plot.hurst()
# --------------------------------------------------------------------------

#' Plot Hurst Analysis Results
#'
#' Visualizes rolling Hurst exponent results with state-colored backgrounds
#' and threshold bands. Supports displaying the original time series, the
#' Hurst trajectory, or both panels stacked.
#'
#' @details
#' Three plot types are available:
#'
#' **`"series"`.**
#' Draws the original time series as a line plot with the background
#' shaded according to the Hurst state classification at each time point.
#' The five states (`strong_antipersistent`, `antipersistent`,
#' `random_walk`, `persistent`, `strong_persistent`) are mapped to a
#' fixed color palette. This view answers the question: "What was the
#' system doing when it was in each long-range dependence regime?"
#'
#' **`"states"`.**
#' Plots the rolling Hurst exponent trajectory as a line, overlaid on
#' horizontal bands at the classification thresholds (0.4, 0.5, 0.6,
#' 0.7). A dashed reference line at \eqn{H = 0.5} marks the boundary
#' between antipersistent and persistent behavior. This view answers:
#' "How did the Hurst exponent evolve, and how close is it to regime
#' boundaries?"
#'
#' **`"both"`** (default).
#' Stacks the series panel on top and the states panel below using
#' [patchwork::wrap_plots()], with shared time axes. The top panel
#' suppresses its x-axis labels to avoid duplication.
#'
#' @export
#' @param x \[`hurst`\]\cr
#'   An object of class `hurst`.
#' @param type \[`character(1)`: `"both"`\]\cr
#'   Plot type. The available options are:
#'   * `"series"`: Time series with state-colored background.
#'   * `"states"`: Hurst exponent trajectory with threshold bands.
#'   * `"both"`: Both panels stacked vertically.
#' @param ... Additional arguments (currently unused).
#' @return A [ggplot2::ggplot()] object (or a
#'   [patchwork::wrap_plots()] composite when `type = "both"`).
#'
#' @references
#' Peng, C.-K., Buldyrev, S.V., Havlin, S., Simons, M., Stanley, H.E.,
#' & Goldberger, A.L. (1994). Mosaic organization of DNA nucleotides.
#' \emph{Physical Review E}, 49(2), 1685--1689.
#' \doi{10.1103/PhysRevE.49.1685}
#'
#' @seealso [hurst()] for computing the rolling Hurst exponent.
#' @family hurst
#' @concept time series
#' @concept hurst exponent
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' h <- hurst(x, window = 50L)
#' plot(h, type = "series")
#' plot(h, type = "states")
#' plot(h, type = "both")
#' }
plot.hurst <- function(x, type = "both", ...) {
  check_missing(x)
  check_class(x, "hurst")
  type <- check_match(type, c("series", "states", "both"))
  colors <- hurst_state_colors()
  all_states <- c(
    "strong_antipersistent", "antipersistent", "random_walk",
    "persistent", "strong_persistent"
  )
  d <- dplyr::mutate(
    x,
    state_f = factor(
      !!rlang::sym("state"),
      levels = all_states
    )
  )
  # Build state-colored background rectangles
  rects <- build_state_rects_(d)
  if (type == "series" || type == "both") {
    val <- d[["value"]]
    val <- val[is.finite(val)]
    y_rng <- range(val, na.rm = TRUE)
    y_pad <- diff(y_rng) * 0.08
    y_lo  <- y_rng[1L] - y_pad
    y_hi  <- y_rng[2L] + y_pad
    rects[["ymin"]] <- y_lo
    rects[["ymax"]] <- y_hi
    p_series <- ggplot2::ggplot(
      d,
      ggplot2::aes(x = !!rlang::sym("time"), y = !!rlang::sym("value"))
    ) +
      ggplot2::geom_rect(
        data = rects,
        ggplot2::aes(
          xmin = !!rlang::sym("xmin"),
          xmax = !!rlang::sym("xmax"),
          ymin = !!rlang::sym("ymin"),
          ymax = !!rlang::sym("ymax"),
          fill = !!rlang::sym("state_f")
        ),
        alpha = 0.25,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_line(linewidth = 0.4) +
      ggplot2::scale_fill_manual(
        values = colors,
        limits = all_states,
        name = "State",
        drop = FALSE
      ) +
      ggplot2::coord_cartesian(
        ylim = c(y_lo, y_hi),
        clip = "off"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = NULL, y = "Value") +
      ggplot2::theme(
        legend.position = "bottom",
        plot.margin = ggplot2::margin(t = 12, r = 5, b = 5, l = 5)
      )
    if (type == "series") return(p_series)
  }
  if (type == "states" || type == "both") {
    h_vals <- d[["hurst"]]
    h_vals <- h_vals[is.finite(h_vals)]
    h_max  <- max(h_vals, 1.0, na.rm = TRUE)
    h_top  <- h_max + 0.05
    p_states <- ggplot2::ggplot(
      d,
      ggplot2::aes(x = !!rlang::sym("time"), y = !!rlang::sym("hurst"))
    ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 0, ymax = 0.4,
        fill = colors["strong_antipersistent"], alpha = 0.15
      ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 0.4, ymax = 0.5,
        fill = colors["antipersistent"], alpha = 0.15
      ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 0.5, ymax = 0.6,
        fill = colors["random_walk"], alpha = 0.15
      ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 0.6, ymax = 0.7,
        fill = colors["persistent"], alpha = 0.15
      ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 0.7, ymax = h_top,
        fill = colors["strong_persistent"], alpha = 0.15
      ) +
      ggplot2::geom_line(linewidth = 0.5) +
      ggplot2::geom_hline(
        yintercept = 0.5,
        linetype = "dashed",
        color = "grey40"
      ) +
      ggplot2::scale_y_continuous(
        breaks = seq(0, ceiling(h_top * 4) / 4, 0.25)
      ) +
      ggplot2::coord_cartesian(
        ylim = c(-0.02, h_top),
        clip = "off"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Time", y = "Hurst Exponent") +
      ggplot2::theme(
        plot.margin = ggplot2::margin(t = 8, r = 5, b = 5, l = 5)
      )
    if (type == "states") return(p_states)
  }
  p_series <- p_series +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  patchwork::wrap_plots(p_series, p_states, ncol = 1L, heights = c(1, 1))
}

# --------------------------------------------------------------------------
# S3: plot.hurst_ews()
# --------------------------------------------------------------------------

#' Plot Hurst Early Warning Signals
#'
#' Creates a three-panel stacked visualization of early warning signals
#' derived from rolling Hurst exponent analysis.
#'
#' @details
#' The plot consists of three vertically stacked panels:
#'
#' **Panel 1 -- Hurst trajectory with warning backgrounds.**
#' The rolling Hurst exponent is drawn as a line. Behind it,
#' colored rectangles shade each time region according to its warning
#' level (green = none, yellow = moderate, orange = high, red =
#' critical). A dashed line at \eqn{H = 0.5} marks the random-walk
#' boundary. This panel links changes in long-range dependence to the
#' composite warning score.
#'
#' **Panel 2 -- Warning level step plot.**
#' A step function shows the discrete warning level (0--4) over time.
#' Sudden jumps reveal the onset and duration of elevated warning
#' periods. Reading this jointly with Panel 1 shows whether warning
#' escalation coincides with Hurst crossing a regime boundary.
#'
#' **Panel 3 -- Indicator heatmap.**
#' Each of the selected binary indicators is displayed as a row, with
#' time on the x-axis. Active indicators (value 1) are shaded with a
#' warm colour gradient; inactive indicators are light grey. This
#' reveals which specific mechanisms (e.g., flickering, trend, spectral
#' shift) contribute to each warning, enabling targeted interpretation
#' rather than relying solely on the composite score.
#'
#' The `indicators` argument controls which rows appear in Panel 3.
#' Subsetting to a few indicators of interest can simplify the display
#' for presentation.
#'
#' @export
#' @param x \[`hurst_ews`\]\cr
#'   An object of class `hurst_ews` as returned by
#'   [detect_hurst_warnings()].
#' @param indicators \[`character()`\]\cr
#'   Which indicators to display in the heatmap panel. Defaults to all
#'   10 indicators.
#' @param ... Additional arguments (currently unused).
#' @return A [patchwork::wrap_plots()] composite of three
#'   [ggplot2::ggplot()] panels.
#'
#' @references
#' Scheffer, M., Bascompte, J., Brock, W.A., Brovkin, V., Carpenter,
#' S.R., Dakos, V., Held, H., van Nes, E.H., Rietkerk, M., &
#' Sugihara, G. (2009). Early-warning signals for critical transitions.
#' \emph{Nature}, 461, 53--59. \doi{10.1038/nature08227}
#'
#' @seealso [detect_hurst_warnings()] for computing the warning signals
#'   displayed by this plot.
#' @family hurst
#' @concept early warning signals
#' @concept hurst exponent
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' h <- hurst(x, window = 50L)
#' ews <- detect_hurst_warnings(h)
#'
#' # Default: all indicators
#' plot(ews)
#'
#' # Subset to specific indicators
#' plot(ews, indicators = c("flickering", "trend_up", "spectral_shift"))
#' }
plot.hurst_ews <- function(x,
                           indicators = NULL,
                           ...) {
  check_missing(x)
  check_class(x, "hurst_ews")
  all_indicators <- c(
    "extreme_low", "extreme_high", "trend_up", "trend_down",
    "high_volatility", "flickering", "variance_ratio",
    "spectral_shift", "autocorr_increase", "state_persistence"
  )
  if (is.null(indicators)) indicators <- all_indicators
  indicators <- check_match(
    indicators, all_indicators, several.ok = TRUE
  )
  warning_colors <- c(
    "0" = "#2ca02c",
    "1" = "#98df8a",
    "2" = "#ffcc00",
    "3" = "#ff8c00",
    "4" = "#d62728"
  )
  # Panel 1: Hurst line with warning background
  warn_rects <- build_warning_rects_(x)
  p_hurst <- ggplot2::ggplot(
    x,
    ggplot2::aes(x = !!rlang::sym("time"), y = !!rlang::sym("hurst"))
  )
  if (nrow(warn_rects) > 0L) {
    p_hurst <- p_hurst +
      ggplot2::geom_rect(
        data = warn_rects,
        ggplot2::aes(
          xmin = !!rlang::sym("xmin"),
          xmax = !!rlang::sym("xmax"),
          ymin = -Inf, ymax = Inf,
          fill = !!rlang::sym("level")
        ),
        alpha = 0.3,
        inherit.aes = FALSE
      )
  }
  p_hurst <- p_hurst +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_hline(
      yintercept = 0.5, linetype = "dashed", color = "grey40"
    ) +
    ggplot2::scale_fill_manual(
      values = warning_colors,
      name = "Warning",
      labels = c("None", "Low", "Moderate", "High", "Critical"),
      drop = FALSE
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = "Hurst") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "bottom"
    )
  # Panel 2: Warning level step plot
  p_level <- ggplot2::ggplot(
    x,
    ggplot2::aes(
      x = !!rlang::sym("time"),
      y = !!rlang::sym("warning_level")
    )
  ) +
    ggplot2::geom_step(linewidth = 0.5, color = "#d62728") +
    ggplot2::scale_y_continuous(
      breaks = 0:4,
      labels = c("None", "Low", "Mod", "High", "Crit"),
      limits = c(-0.2, 4.2)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = "Warning Level") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  # Panel 3: Indicator heatmap
  indicator_display <- c(
    extreme_low       = "Extreme Low",
    extreme_high      = "Extreme High",
    trend_up          = "Trend Up",
    trend_down        = "Trend Down",
    high_volatility   = "Volatility",
    flickering        = "Flickering",
    variance_ratio    = "Variance Ratio",
    spectral_shift    = "Spectral Shift",
    autocorr_increase = "Autocorrelation",
    state_persistence = "Persistence"
  )
  long <- tidyr::pivot_longer(
    dplyr::select(x, !!rlang::sym("time"),
                  tidyselect::all_of(indicators)),
    cols = tidyselect::all_of(indicators),
    names_to = "indicator",
    values_to = "intensity"
  )
  long <- dplyr::mutate(
    long,
    indicator_label = factor(
      indicator_display[!!rlang::sym("indicator")],
      levels = rev(indicator_display[indicators])
    )
  )
  p_heat <- ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = !!rlang::sym("time"),
      y = !!rlang::sym("indicator_label"),
      fill = !!rlang::sym("intensity")
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colours = c("#f0f0f0", "#ffcc00", "#ff8c00", "#d62728"),
      limits = c(0, 1),
      na.value = "#f0f0f0",
      name = "Intensity"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Time", y = NULL) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank()
    )
  patchwork::wrap_plots(
    p_hurst, p_level, p_heat,
    ncol = 1L,
    heights = c(2, 1, 3)
  )
}

# --------------------------------------------------------------------------
# Internal plot helpers (@noRd)
# --------------------------------------------------------------------------

#' @noRd
build_state_rects_ <- function(d) {
  n <- nrow(d)
  if (n == 0L) {
    return(data.frame(
      xmin = numeric(0), xmax = numeric(0),
      state_f = factor(character(0))
    ))
  }
  states <- d$state
  time <- d$time
  # Find run-length encoding of states
  rle_states <- rle(states)
  n_runs <- length(rle_states$lengths)
  rects <- data.frame(
    xmin = numeric(n_runs),
    xmax = numeric(n_runs),
    state_f = factor(
      rle_states$values,
      levels = c(
        "strong_antipersistent", "antipersistent", "random_walk",
        "persistent", "strong_persistent"
      )
    )
  )
  pos <- 1L
  for (i in seq_len(n_runs)) {
    run_start <- pos
    run_end <- pos + rle_states$lengths[i] - 1L
    rects$xmin[i] <- time[run_start]
    rects$xmax[i] <- time[run_end]
    pos <- run_end + 1L
  }
  rects
}

#' @noRd
build_warning_rects_ <- function(x) {
  n <- nrow(x)
  if (n == 0L) {
    return(data.frame(
      xmin = numeric(0), xmax = numeric(0),
      level = factor(character(0))
    ))
  }
  levels_vec <- x$warning_level
  time <- x$time
  rle_levels <- rle(levels_vec)
  n_runs <- length(rle_levels$lengths)
  rects <- data.frame(
    xmin = numeric(n_runs),
    xmax = numeric(n_runs),
    level = factor(
      as.character(rle_levels$values),
      levels = as.character(0:4)
    )
  )
  pos <- 1L
  for (i in seq_len(n_runs)) {
    run_start <- pos
    run_end <- pos + rle_levels$lengths[i] - 1L
    rects$xmin[i] <- time[run_start]
    rects$xmax[i] <- time[run_end]
    pos <- run_end + 1L
  }
  rects
}

# --------------------------------------------------------------------------
# S3: print / summary methods
# --------------------------------------------------------------------------

#' Print a Hurst Rolling Analysis Object
#'
#' Prints the underlying tibble of rolling Hurst exponent results. Delegates
#' to the default tibble print method.
#'
#' @export
#' @param x \[`hurst`\]\cr
#'   A `hurst` object as returned by [hurst()] with `states = TRUE`.
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
#' @family hurst
print.hurst <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Print a Global Hurst Exponent Result
#'
#' Prints a concise summary of a single global Hurst exponent estimate,
#' including the method, series length, exponent value, \eqn{R^2}, and
#' state classification.
#'
#' @export
#' @param x \[`hurst_global`\]\cr
#'   A `hurst_global` object as returned by [hurst()] with
#'   `states = FALSE`.
#' @param ... Additional arguments (currently unused).
#' @return `x`, invisibly.
#' @family hurst
print.hurst_global <- function(x, ...) {
  cat("Hurst Exponent Analysis\n")
  cat("  Method :", x$method, "\n")
  cat("  N      :", x$n, "\n")
  cat("  Hurst  :", round(x$hurst, 4L), "\n")
  cat("  R\u00b2     :", round(x$r_squared, 4L), "\n")
  cat("  State  :", hurst_classify(x$hurst), "\n")
  invisible(x)
}

#' Print Hurst Early Warning Signal Results
#'
#' Prints the underlying tibble of early warning signal results. Delegates
#' to the default tibble print method.
#'
#' @export
#' @param x \[`hurst_ews`\]\cr
#'   A `hurst_ews` object as returned by [detect_hurst_warnings()].
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
#' @family hurst
print.hurst_ews <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Summarize Hurst Early Warning Signal Results
#'
#' Prints a compact summary of the warning level distribution across all
#' time points, including counts and percentages for each level and the
#' maximum warning level observed.
#'
#' @export
#' @param object \[`hurst_ews`\]\cr
#'   A `hurst_ews` object as returned by [detect_hurst_warnings()].
#' @param ... Additional arguments (currently unused).
#' @return `object`, invisibly.
#' @family hurst
#' @concept early warning signals
summary.hurst_ews <- function(object, ...) {
  labels <- c("none", "low", "moderate", "high", "critical")
  counts <- table(factor(object$warning_label, levels = labels))
  cat("Hurst Early Warning Signal Summary\n")
  cat("  N time points:", nrow(object), "\n")
  cat("  Warning distribution:\n")
  for (i in seq_along(labels)) {
    pct <- round(100 * counts[i] / nrow(object), 1L)
    cat("    ", labels[i], ": ", counts[i], " (", pct, "%)\n", sep = "")
  }
  max_level <- max(object$warning_level, na.rm = TRUE)
  cat("  Max warning level:", labels[max_level + 1L], "\n")
  invisible(object)
}
