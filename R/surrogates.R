# ============================================================================
# surrogates.R -- Surrogate-based significance testing for EWS trends
# Tests whether observed trends in early warning signal metrics are
# statistically distinguishable from trends expected under a null model.
# ============================================================================

# --------------------------------------------------------------------------
# Internal helpers (@noRd)
# --------------------------------------------------------------------------

#' Generate a single surrogate via FFT phase randomization
#'
#' Randomizes the phases of the Fourier transform while preserving the
#' power spectrum and conjugate symmetry for real-valued output.
#'
#' @param x Numeric vector.
#' @return Numeric vector of the same length.
#' @noRd
surrogate_phase_ <- function(x) {
  n <- length(x)
  ft <- stats::fft(x)
  phases <- stats::runif(n, 0, 2 * pi)
  # Preserve conjugate symmetry: freq 0 and Nyquist (if n is even)
  # keep original phase; mirror the rest
  phases[1L] <- 0
  if (n %% 2L == 0L) {
    phases[n / 2L + 1L] <- 0
    idx <- seq(2L, n / 2L)
    phases[n - idx + 2L] <- -phases[idx]
  } else {
    idx <- seq(2L, (n + 1L) / 2L)
    phases[n - idx + 2L] <- -phases[idx]
  }
  randomized <- ft * exp(1i * phases)
  Re(stats::fft(randomized, inverse = TRUE)) / n
}

#' Generate a single surrogate via random permutation
#'
#' @param x Numeric vector.
#' @return Numeric vector of the same length.
#' @noRd
surrogate_shuffle_ <- function(x) {
  sample(x)
}

#' Generate a single surrogate via AAFT
#'
#' Amplitude-Adjusted Fourier Transform: preserves both the amplitude
#' distribution and approximately the power spectrum of the original.
#'
#' @param x Numeric vector.
#' @return Numeric vector of the same length.
#' @noRd
surrogate_aaft_ <- function(x) {
  n <- length(x)
  # Step 1: Create Gaussian series with the same rank order as x
  gaussian <- stats::rnorm(n)
  gaussian_sorted <- sort(gaussian)
  rank_x <- rank(x, ties.method = "random")
  gaussian_ranked <- gaussian_sorted[rank_x]
  # Step 2: Phase-randomize the Gaussian series
  phase_rand <- surrogate_phase_(gaussian_ranked)
  # Step 3: Reorder the original values to match the rank order
  #   of the phase-randomized Gaussian series
  x_sorted <- sort(x)
  rank_phase <- rank(phase_rand, ties.method = "random")
  x_sorted[rank_phase]
}

#' Generate a single surrogate via block bootstrap
#'
#' Samples contiguous blocks of length sqrt(n) with replacement to
#' build a surrogate of the same length.
#'
#' @param x Numeric vector.
#' @return Numeric vector of the same length.
#' @noRd
surrogate_block_ <- function(x) {
  n <- length(x)
  block_len <- max(2L, as.integer(floor(sqrt(n))))
  n_blocks <- ceiling(n / block_len)
  # Sample starting indices for each block
  max_start <- n - block_len + 1L
  starts <- sample.int(max_start, size = n_blocks, replace = TRUE)
  # Build surrogate from concatenated blocks
  surr <- unlist(lapply(starts, function(s) x[s:(s + block_len - 1L)]))
  surr[seq_len(n)]
}

#' Compute a rolling metric over windows of a time series
#'
#' @param x Numeric vector.
#' @param window Integer window size.
#' @param metric Character metric name.
#' @return Numeric vector of rolling metric values (length = length(x) - window + 1).
#' @noRd
surrogate_rolling_metric_ <- function(x, window, metric) {
  n <- length(x)
  n_windows <- n - window + 1L
  if (n_windows < 1L) return(numeric(0L))
  vapply(seq_len(n_windows), function(i) {
    seg <- x[i:(i + window - 1L)]
    surrogate_compute_metric_(seg, metric)
  }, numeric(1L))
}

#' Compute a single metric value on a segment
#'
#' @param seg Numeric vector (window of data).
#' @param metric Character metric name.
#' @return Single numeric value.
#' @noRd
surrogate_compute_metric_ <- function(seg, metric) {
  seg <- seg[!is.na(seg)]
  if (length(seg) < 3L) return(NA_real_)
  switch(
    metric,
    ar1 = {
      acf_val <- try_(
        stats::acf(seg, lag.max = 1L, plot = FALSE)$acf[2L]
      )
      if (inherits(acf_val, "try-error") || is.na(acf_val)) {
        NA_real_
      } else {
        acf_val
      }
    },
    sd = stats::sd(seg),
    skewness = {
      m <- mean(seg)
      s <- stats::sd(seg)
      if (s < .Machine$double.eps) return(NA_real_)
      mean(((seg - m) / s)^3)
    },
    kurtosis = {
      m <- mean(seg)
      s <- stats::sd(seg)
      if (s < .Machine$double.eps) return(NA_real_)
      mean(((seg - m) / s)^4) - 3
    },
    spectral_exponent = {
      sp <- try_(stats::spectrum(seg, method = "pgram", plot = FALSE))
      if (inherits(sp, "try-error")) return(NA_real_)
      valid <- sp$freq > 0 & sp$spec > 0
      if (sum(valid) < 3L) return(NA_real_)
      log_freq <- log(sp$freq[valid])
      log_spec <- log(sp$spec[valid])
      fit <- stats::lm.fit(x = cbind(1, log_freq), y = log_spec)
      fit$coefficients[2L]
    },
    hurst = {
      surrogate_hurst_dfa_(seg)
    },
    NA_real_
  )
}

#' Simplified DFA Hurst exponent for a single segment
#'
#' Uses a compact DFA implementation suitable for short rolling windows.
#'
#' @param x Numeric vector.
#' @return Single numeric Hurst exponent estimate.
#' @noRd
surrogate_hurst_dfa_ <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  min_scale <- 4L
  max_scale <- floor(n / 4)
  if (max_scale <= min_scale || n < 2L * min_scale) return(NA_real_)
  n_scales <- min(5L, max_scale - min_scale + 1L)
  scales <- unique(round(
    exp(seq(log(min_scale), log(max_scale), length.out = n_scales))
  ))
  scales <- scales[scales >= min_scale & scales <= max_scale]
  if (length(scales) < 3L) return(NA_real_)
  profile <- cumsum(x - mean(x))
  fluct <- vapply(scales, function(s) {
    n_seg <- floor(n / s)
    if (n_seg < 2L) return(NA_real_)
    x_vals <- seq_len(s)
    x_mean <- mean(x_vals)
    xx_var <- sum((x_vals - x_mean)^2)
    if (xx_var == 0) return(NA_real_)
    seg_fluct <- vapply(seq_len(n_seg), function(j) {
      start_idx <- (j - 1L) * s + 1L
      end_idx <- j * s
      segment <- profile[start_idx:end_idx]
      y_mean <- mean(segment)
      xy_cov <- sum((x_vals - x_mean) * (segment - y_mean))
      slope <- xy_cov / xx_var
      intercept <- y_mean - slope * x_mean
      trend <- intercept + slope * x_vals
      sqrt(mean((segment - trend)^2))
    }, numeric(1L))
    sqrt(mean(seg_fluct^2))
  }, numeric(1L))
  valid <- !is.na(fluct) & fluct > 0
  if (sum(valid) < 3L) return(NA_real_)
  log_s <- log(scales[valid])
  log_f <- log(fluct[valid])
  fit <- stats::lm.fit(x = cbind(1, log_s), y = log_f)
  fit$coefficients[2L]
}

#' Compute Kendall tau of a metric vector against its time index
#'
#' @param metric_values Numeric vector.
#' @return Single numeric Kendall tau.
#' @noRd
surrogate_kendall_tau_ <- function(metric_values) {
  valid <- !is.na(metric_values)
  if (sum(valid) < 4L) return(NA_real_)
  vals <- metric_values[valid]
  time_idx <- which(valid)
  result <- try_(
    stats::cor.test(time_idx, vals, method = "kendall")$estimate
  )
  if (inherits(result, "try-error")) NA_real_ else result
}

#' Colour palette for surrogate test plots
#'
#' @return Named character vector of colours.
#' @noRd
surrogate_colors_ <- function() {
  c(
    observed   = "#E63946",
    surrogates = "#457B9D",
    envelope   = "#A8DADC",
    fill       = "#1D3557"
  )
}

#' Format a p-value for display, avoiding ugly 0.0000
#' @param p Numeric p-value.
#' @return Character string.
#' @noRd
surrogate_format_p_ <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("p < 0.001")
  sprintf("p = %.3f", p)
}

# --------------------------------------------------------------------------
# Exported: surrogate_test()
# --------------------------------------------------------------------------

#' Surrogate-Based Significance Testing for EWS Trends
#'
#' Tests whether the observed trend in an early warning signal (EWS) metric
#' is statistically significant by comparing it against a null distribution
#' generated from surrogate time series. Four surrogate generation methods
#' are available, covering a range of null hypotheses from uncorrelated
#' noise to nonlinear structure preservation.
#'
#' @details
#' The surrogate testing procedure follows three steps:
#'
#' **1. Surrogate generation.** Multiple surrogate series are created from
#' the original data under a chosen null hypothesis. The four methods differ
#' in what structure they preserve:
#'
#' \describe{
#'   \item{`"phase"` (FFT phase randomization)}{Preserves the power spectrum
#'     (and hence the autocorrelation function) of the original series but
#'     destroys all phase relationships. This tests the null hypothesis that
#'     the observed trend is consistent with a linear Gaussian process with
#'     the same spectral density. Conjugate symmetry of the Fourier
#'     coefficients is maintained to ensure real-valued output.}
#'   \item{`"shuffle"` (random permutation)}{Destroys all temporal structure
#'     (autocorrelation, trend, periodicity). This tests the most permissive
#'     null: that the observed trend could arise from any random rearrangement
#'     of the same values.}
#'   \item{`"aaft"` (Amplitude-Adjusted Fourier Transform)}{Preserves both
#'     the amplitude distribution and approximately the power spectrum of
#'     the original series. This tests the null hypothesis that the data
#'     were generated by a static nonlinear transformation of a linear
#'     Gaussian process (Theiler et al., 1992). Preferable to pure phase
#'     randomization when the original series is non-Gaussian.}
#'   \item{`"block"` (block bootstrap)}{Preserves local temporal structure
#'     within blocks of length \eqn{\sqrt{n}} while shuffling block order.
#'     This tests whether the observed long-range trend is distinguishable
#'     from short-range autocorrelation patterns rearranged randomly.}
#' }
#'
#' **2. Metric and trend computation.** For the original series and each
#' surrogate, a rolling-window metric is computed across the time series.
#' The available metrics are:
#' \describe{
#'   \item{`"ar1"`}{Lag-1 autocorrelation. Rising AR(1) is a hallmark of
#'     critical slowing down.}
#'   \item{`"sd"`}{Standard deviation. Increasing variance signals growing
#'     instability.}
#'   \item{`"skewness"`}{Third standardized moment. Asymmetric fluctuations
#'     may indicate approaching a bifurcation.}
#'   \item{`"kurtosis"`}{Excess kurtosis (fourth standardized moment minus
#'     3). Heavy tails signal intermittent extreme events.}
#'   \item{`"spectral_exponent"`}{Slope of the log-power vs. log-frequency
#'     regression. More negative slopes (reddening) indicate critical
#'     slowing down.}
#'   \item{`"hurst"`}{Hurst exponent via simplified DFA. Rising Hurst
#'     indicates increasing long-range dependence.}
#' }
#'
#' A trend statistic (Kendall's tau of the metric values against time)
#' summarizes each metric trajectory into a single scalar.
#'
#' **3. Significance assessment.** The p-value is the proportion of
#' surrogate Kendall tau values that are greater than or equal to the
#' observed tau (one-sided test). A result is flagged as significant at
#' the 0.05 level.
#'
#' @export
#' @param data \[`ts`, `numeric()`\]\cr
#'   Univariate time series data.
#' @param n_surrogates \[`integer(1)`: `100L`\]\cr
#'   Number of surrogate series to generate.
#' @param method \[`character(1)`: `"phase"`\]\cr
#'   Surrogate generation method. The available options are:
#'   * `"phase"`: FFT phase randomization -- preserves the power spectrum.
#'   * `"shuffle"`: Random permutation -- destroys all temporal structure.
#'   * `"aaft"`: Amplitude-Adjusted Fourier Transform -- preserves both
#'     amplitude distribution and spectral density.
#'   * `"block"`: Block bootstrap with block length sqrt(n) -- preserves
#'     local autocorrelation structure.
#' @param metric \[`character(1)`: `"ar1"`\]\cr
#'   Rolling-window metric to track. The available options are:
#'   `"ar1"`, `"sd"`, `"skewness"`, `"kurtosis"`,
#'   `"spectral_exponent"`, and `"hurst"`.
#' @param window \[`integer(1)`: `50L`\]\cr
#'   Rolling window size for computing the metric.
#' @param test \[`character(1)`: `"trend"`\]\cr
#'   Test statistic. Currently only `"trend"` (Kendall's tau of
#'   the rolling metric values against time) is supported.
#' @return An object of class `"surrogate_test"` (a list) containing:
#'   \describe{
#'     \item{`observed_tau`}{Kendall's tau of the rolling metric for the
#'       original series.}
#'     \item{`observed_metric`}{Numeric vector of rolling metric values
#'       for the original series.}
#'     \item{`surrogate_taus`}{Numeric vector of Kendall's tau for each
#'       surrogate series.}
#'     \item{`surrogate_metrics`}{Matrix with `n_surrogates` rows and
#'       one column per rolling-window position, containing metric values
#'       for all surrogates.}
#'     \item{`p_value`}{Proportion of surrogate taus >= observed tau
#'       (one-sided).}
#'     \item{`significant`}{Logical; `TRUE` if `p_value < 0.05`.}
#'     \item{`metric`}{Character name of the metric used.}
#'     \item{`method`}{Character name of the surrogate method used.}
#'     \item{`n_surrogates`}{Integer number of surrogates generated.}
#'     \item{`window`}{Integer window size used.}
#'     \item{`test`}{Character test statistic name.}
#'   }
#'
#'   All analysis settings are also stored as attributes for retrieval.
#'
#' @references
#' Theiler, J., Eubank, S., Longtin, A., Galdrikian, B., & Farmer, J. D.
#' (1992). Testing for nonlinearity in time series: the method of surrogate
#' data. \emph{Physica D: Nonlinear Phenomena}, 58(1--4), 77--94.
#' \doi{10.1016/0167-2789(92)90102-S}
#'
#' Dakos, V., Carpenter, S. R., Brock, W. A., Ellison, A. M., Guttal, V.,
#' Ives, A. R., ... & Scheffer, M. (2012). Methods for detecting early
#' warnings of critical transitions in time series illustrated using
#' simulated ecological data. \emph{PLoS ONE}, 7(7), e41010.
#' \doi{10.1371/journal.pone.0041010}
#'
#' Lancaster, G., Iatsenko, D., Pidde, A., Ticcinelli, V., &
#' Stefanovska, A. (2018). Surrogate data for hypothesis testing of
#' physical systems. \emph{Physics Reports}, 748, 1--60.
#' \doi{10.1016/j.physrep.2018.06.001}
#'
#' @seealso [hurst()] for Hurst exponent analysis;
#'   [resilience()] for rolling resilience metrics;
#'   [compute_trend()] for trend classification.
#' @family surrogates
#' @concept early warning signals
#' @concept significance testing
#' @concept surrogate data
#' @examples
#' \donttest{
#' set.seed(42)
#' # Series with an increasing trend (should show significant AR1 trend)
#' x <- cumsum(rnorm(300)) + seq(0, 5, length.out = 300)
#' st <- surrogate_test(x, n_surrogates = 50L, metric = "ar1", window = 30L)
#' print(st)
#' summary(st)
#' plot(st)
#' }
surrogate_test <- function(data,
                           n_surrogates = 100L,
                           method = "phase",
                           metric = "ar1",
                           window = 50L,
                           test = "trend") {
  # --- Validation ---
  check_missing(data)
  data <- prepare_timeseries_data(data)
  values <- data$values
  n <- length(values)
  method <- check_match(method, c("phase", "shuffle", "aaft", "block"))
  metric <- check_match(
    metric,
    c("ar1", "sd", "skewness", "kurtosis", "spectral_exponent", "hurst")
  )
  test <- check_match(test, "trend")
  check_range(n_surrogates, type = "integer", min = 10L, max = 100000L)
  check_range(window, type = "integer", min = 2L, max = n - 1L)

  if (n < window + 4L) {
    stop_(
      "Series length ({n}) must be at least {.arg window} + 4 = {window + 4L}
      for meaningful rolling metric computation."
    )
  }

  # --- Select surrogate generator ---
  generate_fn <- switch(
    method,
    phase   = surrogate_phase_,
    shuffle = surrogate_shuffle_,
    aaft    = surrogate_aaft_,
    block   = surrogate_block_
  )

  # --- Compute rolling metric for original series ---
  observed_metric <- surrogate_rolling_metric_(values, window, metric)
  observed_tau <- surrogate_kendall_tau_(observed_metric)

  # --- Generate surrogates and compute their metrics/taus ---
  n_windows <- length(observed_metric)
  surrogate_taus <- numeric(n_surrogates)
  surrogate_metrics <- matrix(NA_real_, nrow = n_surrogates, ncol = n_windows)

  for (i in seq_len(n_surrogates)) {
    surr <- generate_fn(values)
    surr_metric <- surrogate_rolling_metric_(surr, window, metric)
    # Ensure same length (block bootstrap may differ slightly)
    surr_len <- min(length(surr_metric), n_windows)
    if (surr_len > 0L) {
      surrogate_metrics[i, seq_len(surr_len)] <- surr_metric[seq_len(surr_len)]
    }
    surrogate_taus[i] <- surrogate_kendall_tau_(surr_metric)
  }

  # --- Compute p-value (one-sided: proportion of surrogates >= observed) ---
  valid_taus <- surrogate_taus[!is.na(surrogate_taus)]
  if (length(valid_taus) == 0L || is.na(observed_tau)) {
    p_value <- NA_real_
    significant <- NA
  } else {
    p_value <- mean(valid_taus >= observed_tau)
    significant <- p_value < 0.05
  }

  # --- Build return object ---
  out <- list(
    observed_tau      = observed_tau,
    observed_metric   = observed_metric,
    surrogate_taus    = surrogate_taus,
    surrogate_metrics = surrogate_metrics,
    p_value           = p_value,
    significant       = significant,
    metric            = metric,
    method            = method,
    n_surrogates      = as.integer(n_surrogates),
    window            = as.integer(window),
    test              = test
  )

  structure(
    out,
    metric       = metric,
    method       = method,
    n_surrogates = as.integer(n_surrogates),
    window       = as.integer(window),
    test         = test,
    class        = c("surrogate_test", "list")
  )
}

# --------------------------------------------------------------------------
# S3 methods
# --------------------------------------------------------------------------

#' Plot Surrogate Test Results
#'
#' Visualizes the surrogate significance test with a histogram of surrogate
#' Kendall tau values, the observed metric trajectory with a surrogate
#' envelope, or both panels stacked.
#'
#' @details
#' Three plot types are available:
#'
#' **`"histogram"`.**
#' Draws a histogram of Kendall's tau values from all surrogates, with a
#' vertical dashed line marking the observed tau. The shaded region to the
#' right of the observed value represents the p-value visually. This panel
#' answers: "How extreme is the observed trend relative to the null
#' distribution?"
#'
#' **`"envelope"`.**
#' Plots the observed rolling metric trajectory (solid red line) overlaid
#' on a shaded 95\% envelope (2.5th to 97.5th percentiles) computed from
#' all surrogates. Points where the observed trajectory escapes the
#' envelope signal time periods where the trend departs from the null.
#' This panel answers: "When during the series does the metric deviate
#' from what surrogates produce?"
#'
#' **`"both"` (default).**
#' Stacks the histogram panel on top and the envelope panel below using
#' [patchwork::wrap_plots()].
#'
#' @export
#' @param x \[`surrogate_test`\]\cr
#'   An object of class `surrogate_test` as returned by
#'   [surrogate_test()].
#' @param type \[`character(1)`: `"both"`\]\cr
#'   Plot type. The available options are:
#'   * `"histogram"`: Null distribution of surrogate taus with observed
#'     value marked.
#'   * `"envelope"`: Observed metric with 95\% surrogate envelope.
#'   * `"both"`: Both panels stacked vertically.
#' @param ... Additional arguments (currently unused).
#' @return A [ggplot2::ggplot()] object (or a [patchwork::wrap_plots()]
#'   composite when `type = "both"`).
#'
#' @seealso [surrogate_test()] for running the significance test.
#' @family surrogates
#' @concept early warning signals
#' @concept significance testing
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(300))
#' st <- surrogate_test(x, n_surrogates = 50L, metric = "ar1", window = 30L)
#' plot(st, type = "histogram")
#' plot(st, type = "envelope")
#' plot(st, type = "both")
#' }
plot.surrogate_test <- function(x, type = "both", ...) {
  check_missing(x)
  check_class(x, "surrogate_test")
  type <- check_match(type, c("histogram", "envelope", "both"))
  colors <- surrogate_colors_()

  metric_labels <- c(
    ar1                = "AR(1)",
    sd                 = "Std. Dev.",
    skewness           = "Skewness",
    kurtosis           = "Kurtosis",
    spectral_exponent  = "Spectral Exponent",
    hurst              = "Hurst Exponent"
  )
  metric_label <- metric_labels[x$metric]

  # --- Histogram panel ---
  if (type == "histogram" || type == "both") {
    tau_df <- data.frame(tau = x$surrogate_taus[!is.na(x$surrogate_taus)])
    obs_df <- data.frame(obs_tau = x$observed_tau)

    p_hist <- ggplot2::ggplot(tau_df, ggplot2::aes(x = !!rlang::sym("tau"))) +
      ggplot2::geom_histogram(
        fill = colors["surrogates"],
        color = "white",
        bins = 30L,
        alpha = 0.8
      ) +
      ggplot2::geom_vline(
        data = obs_df,
        ggplot2::aes(xintercept = !!rlang::sym("obs_tau")),
        color = colors["observed"],
        linewidth = 1.2,
        linetype = "dashed"
      ) +
      ggplot2::annotate(
        "label",
        x = Inf,
        y = Inf,
        label = sprintf("tau = %.3f\n%s",
                        x$observed_tau, surrogate_format_p_(x$p_value)),
        hjust = 1.05,
        vjust = 1.2,
        color = colors["observed"],
        fill = "white",
        label.size = 0,
        size = 3.8,
        fontface = "bold"
      ) +
      ggplot2::labs(
        title = sprintf(
          "Surrogate Test: %s trend (%s method)",
          metric_label, x$method
        ),
        x = "Kendall's tau",
        y = "Count"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 12, face = "bold")
      )

    if (type == "histogram") return(p_hist)
  }

  # --- Envelope panel ---
  if (type == "envelope" || type == "both") {
    n_windows <- length(x$observed_metric)
    time_idx <- seq_len(n_windows)

    # Compute 2.5th and 97.5th percentiles across surrogates
    surr_mat <- x$surrogate_metrics
    lower_env <- apply(surr_mat, 2L, stats::quantile,
                       probs = 0.025, na.rm = TRUE)
    upper_env <- apply(surr_mat, 2L, stats::quantile,
                       probs = 0.975, na.rm = TRUE)
    median_env <- apply(surr_mat, 2L, stats::median, na.rm = TRUE)

    env_df <- data.frame(
      time = time_idx,
      observed = x$observed_metric,
      lower = lower_env,
      upper = upper_env,
      median = median_env
    )

    p_env <- ggplot2::ggplot(env_df, ggplot2::aes(x = !!rlang::sym("time"))) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = !!rlang::sym("lower"),
          ymax = !!rlang::sym("upper")
        ),
        fill = colors["envelope"],
        alpha = 0.5
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = !!rlang::sym("median")),
        color = colors["surrogates"],
        linewidth = 0.5,
        linetype = "dotted"
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = !!rlang::sym("observed")),
        color = colors["observed"],
        linewidth = 0.8
      ) +
      ggplot2::labs(
        title = sprintf(
          "Rolling %s with 95%% surrogate envelope",
          metric_label
        ),
        subtitle = sprintf(
          "Window = %d | %s surrogates (%s method)",
          x$window, x$n_surrogates, x$method
        ),
        x = "Window position",
        y = metric_label
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 12, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 9, color = "grey40")
      )

    if (type == "envelope") return(p_env)
  }

  # --- Both panels ---
  p_hist <- p_hist +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    )
  patchwork::wrap_plots(p_hist, p_env, ncol = 1L, heights = c(1, 1.5))
}

#' @describeIn surrogate_test Print method for `"surrogate_test"` objects.
#'   Displays the observed tau, p-value, significance, and method settings.
#'
#' @export
#' @param x \[`surrogate_test`\]\cr
#'   A `surrogate_test` object.
#' @param ... Additional arguments (currently unused).
#' @return `x`, invisibly.
print.surrogate_test <- function(x, ...) {
  check_missing(x)
  check_class(x, "surrogate_test")
  cat("Surrogate Significance Test\n")
  cat(strrep("-", 35L), "\n")
  cat("  Method      :", x$method, "\n")
  cat("  Metric      :", x$metric, "\n")
  cat("  Window      :", x$window, "\n")
  cat("  Surrogates  :", x$n_surrogates, "\n")
  cat("  Observed tau:", round(x$observed_tau, 4L), "\n")
  cat("  p-value     :", surrogate_format_p_(x$p_value), "\n")
  sig_label <- ifelse_(
    is.na(x$significant), "NA",
    ifelse_(x$significant, "YES (p < 0.05)", "NO (p >= 0.05)")
  )
  cat("  Significant :", sig_label, "\n")
  invisible(x)
}

#' @describeIn surrogate_test Summary method for `"surrogate_test"` objects.
#'   Provides detailed statistics of the surrogate distribution including
#'   quantiles and the position of the observed value.
#'
#' @export
#' @param object \[`surrogate_test`\]\cr
#'   A `surrogate_test` object.
#' @param ... Additional arguments (currently unused).
#' @return A list with elements `observed_tau`, `p_value`, `significant`,
#'   `surrogate_quantiles` (named numeric vector), `surrogate_mean`,
#'   `surrogate_sd`, `method`, `metric`, `window`, and `n_surrogates`,
#'   returned invisibly.
summary.surrogate_test <- function(object, ...) {
  check_missing(object)
  check_class(object, "surrogate_test")

  valid_taus <- object$surrogate_taus[!is.na(object$surrogate_taus)]
  surr_quantiles <- stats::quantile(
    valid_taus, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1),
    na.rm = TRUE
  )
  surr_mean <- mean(valid_taus, na.rm = TRUE)
  surr_sd <- stats::sd(valid_taus, na.rm = TRUE)

  cat("Surrogate Test Summary\n")
  cat(strrep("=", 45L), "\n")
  cat("  Method         :", object$method, "\n")
  cat("  Metric         :", object$metric, "\n")
  cat("  Window         :", object$window, "\n")
  cat("  Test           :", object$test, "\n")
  cat("  N surrogates   :", object$n_surrogates, "\n")
  cat("  Valid surrogates:", length(valid_taus), "\n\n")

  cat("Observed:\n")
  cat("  Kendall's tau  :", round(object$observed_tau, 4L), "\n")
  cat("  p-value        :", surrogate_format_p_(object$p_value), "\n")
  sig_label <- ifelse_(
    is.na(object$significant), "NA",
    ifelse_(object$significant, "YES (p < 0.05)", "NO (p >= 0.05)")
  )
  cat("  Significant    :", sig_label, "\n\n")

  cat("Surrogate distribution:\n")
  cat("  Mean           :", round(surr_mean, 4L), "\n")
  cat("  SD             :", round(surr_sd, 4L), "\n")
  cat("  Quantiles:\n")
  q_names <- c("Min", "2.5%", "25%", "50%", "75%", "97.5%", "Max")
  for (i in seq_along(surr_quantiles)) {
    cat(sprintf("    %-6s : %8.4f\n", q_names[i], surr_quantiles[i]))
  }

  # Rank of observed tau within surrogates
  if (length(valid_taus) > 0L && !is.na(object$observed_tau)) {
    rank_pct <- mean(valid_taus <= object$observed_tau) * 100
    cat(sprintf(
      "\n  Observed tau ranks at the %.1f%% percentile of surrogates.\n",
      rank_pct
    ))
  }

  invisible(list(
    observed_tau        = object$observed_tau,
    p_value             = object$p_value,
    significant         = object$significant,
    surrogate_quantiles = surr_quantiles,
    surrogate_mean      = surr_mean,
    surrogate_sd        = surr_sd,
    method              = object$method,
    metric              = object$metric,
    window              = object$window,
    n_surrogates        = object$n_surrogates
  ))
}
