# ============================================================================
# spectral.R -- Rolling spectral early warning signals
# Computes rolling power spectral density metrics (spectral exponent, spectral
# ratio) to detect shifts in noise colour that precede critical transitions.
# ============================================================================

# --------------------------------------------------------------------------
# Internal helpers (@noRd)
# --------------------------------------------------------------------------

#' Detrend a numeric segment
#'
#' Applies detrending to a numeric vector prior to spectral estimation.
#' Supports no detrending, linear detrending (OLS residuals), and
#' first-differencing.
#'
#' @param x Numeric vector.
#' @param method One of `"none"`, `"linear"`, or `"diff"`.
#' @return Detrended numeric vector.
#' @noRd
spectral_detrend_ <- function(x, method) {
  switch(
    method,
    none = x,
    linear = {
      n <- length(x)
      t_idx <- seq_len(n)
      t_mean <- mean(t_idx)
      x_mean <- mean(x, na.rm = TRUE)
      xx_var <- sum((t_idx - t_mean)^2)
      if (xx_var > 0) {
        xy_cov <- sum((t_idx - t_mean) * (x - x_mean))
        slope <- xy_cov / xx_var
        intercept <- x_mean - slope * t_mean
        x - (intercept + slope * t_idx)
      } else {
        x - x_mean
      }
    },
    diff = diff(x)
  )
}

#' Compute spectral metrics for a single window
#'
#' Fits a log-log regression to the power spectrum and computes the spectral
#' exponent (beta), spectral ratio (low/high frequency power), and R-squared.
#'
#' @param x Numeric vector (window data, already detrended).
#' @param spec_method One of `"periodogram"` or `"ar"`.
#' @return A named list with elements `beta`, `ratio`, `r_squared`, or `NULL`
#'   on failure.
#' @noRd
spectral_metrics_ <- function(x, spec_method) {
  if (length(x) < 10L || stats::var(x, na.rm = TRUE) < 1e-12) return(NULL)
  spec_result <- if (spec_method == "periodogram") {
    try_(stats::spectrum(x, method = "pgram", plot = FALSE, detrend = FALSE,
                         taper = 0.1))
  } else {
    try_(stats::spectrum(x, method = "ar", plot = FALSE))
  }
  if (inherits(spec_result, "try-error") || is.null(spec_result)) return(NULL)
  freq <- spec_result$freq
  power <- spec_result$spec
  # Remove zero or negative frequencies/power for log transform
  valid <- freq > 0 & power > 0
  if (sum(valid) < 4L) return(NULL)
  freq <- freq[valid]
  power <- power[valid]
  log_freq <- log(freq)
  log_power <- log(power)
  # OLS: log(power) ~ log(freq)
  x_m <- mean(log_freq)
  y_m <- mean(log_power)
  xx_var <- sum((log_freq - x_m)^2)
  if (xx_var < 1e-12) return(NULL)
  xy_cov <- sum((log_freq - x_m) * (log_power - y_m))
  slope <- xy_cov / xx_var
  intercept <- y_m - slope * x_m
  y_pred <- intercept + slope * log_freq
  ss_tot <- sum((log_power - y_m)^2)
  ss_res <- sum((log_power - y_pred)^2)
  r_sq <- ifelse_(ss_tot > 0, 1 - ss_res / ss_tot, NA_real_)
  # Spectral exponent: negate slope so positive beta = reddening
  beta <- -slope
  # Spectral ratio: sum(power at low freq) / sum(power at high freq)
  median_freq <- stats::median(freq)
  low_power <- sum(power[freq <= median_freq])
  high_power <- sum(power[freq > median_freq])
  ratio <- ifelse_(high_power > 0, low_power / high_power, NA_real_)
  list(beta = beta, ratio = ratio, r_squared = r_sq)
}

#' Classify spectral exponent into noise-colour states
#'
#' @param beta Numeric vector of spectral exponents.
#' @return Character vector of state labels.
#' @noRd
spectral_classify_ <- function(beta) {
  labels <- c("white_noise", "pink_noise", "red_noise", "brownian")
  breaks <- c(-Inf, 0.5, 1.5, 2.5, Inf)
  as.character(cut(beta, breaks = breaks, labels = labels, right = FALSE))
}

#' Colour palette for spectral noise states
#' @noRd
spectral_state_colors_ <- function() {
  c(
    white_noise = "#FFD93D",
    pink_noise  = "#FFA06B",
    red_noise   = "#FF6B6B",
    brownian    = "#9B59B6"
  )
}

#' Build state-coloured rectangles for spectral plots
#'
#' @param d Tibble with `time` and `state_f` columns.
#' @return A data.frame with `xmin`, `xmax`, `state_f` columns.
#' @noRd
spectral_build_rects_ <- function(d) {
  n <- nrow(d)
  if (n == 0L) {
    return(data.frame(
      xmin = numeric(0), xmax = numeric(0),
      state_f = factor(character(0))
    ))
  }
  states <- as.character(d$state_f)
  time <- d$time
  rle_states <- rle(states)
  n_runs <- length(rle_states$lengths)
  rects <- data.frame(
    xmin = numeric(n_runs),
    xmax = numeric(n_runs),
    state_f = factor(
      rle_states$values,
      levels = c("white_noise", "pink_noise", "red_noise", "brownian")
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

# --------------------------------------------------------------------------
# Exported: spectral_ews()
# --------------------------------------------------------------------------

#' Rolling Spectral Early Warning Signals
#'
#' Computes rolling-window spectral metrics for a univariate time series,
#' estimating the spectral exponent (beta) and spectral ratio to detect
#' shifts in noise colour that precede critical transitions. Reddening of
#' the power spectrum (increasing beta) is a hallmark of critical slowing
#' down.
#'
#' @details
#' The analysis proceeds in four stages:
#'
#' **1. Window extraction.** A rolling window of size `window` slides
#' across the series with the specified `align`ment. Each window
#' position produces one set of spectral metrics.
#'
#' **2. Detrending.** Before computing the spectrum, each window segment
#' is detrended according to `detrend`:
#' \describe{
#'   \item{`"none"`}{No detrending. Suitable when the data are already
#'     stationary within each window.}
#'   \item{`"linear"` (default)}{Removes a least-squares linear trend
#'     from each window. Appropriate for slowly drifting series.}
#'   \item{`"diff"`}{First-differences the window segment. Effective
#'     for integrated (unit-root) processes but reduces the segment
#'     length by one observation.}
#' }
#'
#' **3. Spectral estimation.** The power spectral density is estimated
#' by one of two methods:
#' \describe{
#'   \item{`"periodogram"` (default)}{Computes the raw periodogram via
#'     [stats::spectrum()] with `method = "pgram"`. Fast and
#'     nonparametric. A 10\% cosine taper is applied by default.}
#'   \item{`"ar"`}{Fits an autoregressive model and derives the spectrum
#'     via [stats::spectrum()] with `method = "ar"`. Smoother estimates,
#'     but assumes the AR model is appropriate.}
#' }
#'
#' **4. Metrics and classification.** For each window:
#' \itemize{
#'   \item The **spectral exponent** (beta) is the negative slope of the
#'     \eqn{\log(\mathrm{power})} vs. \eqn{\log(\mathrm{frequency})}
#'     regression. Positive beta indicates reddened spectra (more low-
#'     frequency power). Rising beta over time signals critical slowing
#'     down.
#'   \item The **spectral ratio** is the ratio of total power at
#'     frequencies below the median to total power above the median.
#'     Values > 1 indicate a red-shifted spectrum.
#'   \item The **R-squared** of the log-log regression quantifies how
#'     well the power-law model fits the spectrum.
#' }
#'
#' **State classification** based on the spectral exponent:
#' \itemize{
#'   \item \eqn{\beta < 0.5}: `white_noise` -- uncorrelated fluctuations
#'     with flat spectrum.
#'   \item \eqn{0.5 \le \beta < 1.5}: `pink_noise` -- \eqn{1/f}-like
#'     dynamics, moderate long-range correlations.
#'   \item \eqn{1.5 \le \beta < 2.5}: `red_noise` -- strong low-frequency
#'     dominance, consistent with critical slowing down.
#'   \item \eqn{\beta \ge 2.5}: `brownian` -- Brownian-motion-like
#'     dynamics with very strong low-frequency power.
#' }
#'
#' **Edge padding.** Window positions at the series boundaries that lack
#' sufficient data are filled by a three-step strategy: (1) forward-fill
#' from the first valid estimate to the series start, (2) backward-fill
#' from the last valid estimate to the series end, and (3) linear
#' interpolation for any internal gaps.
#'
#' @export
#' @param data \[`ts`, `numeric()`\]\cr
#'   Univariate time series data. Accepts a numeric vector or a `ts` object.
#' @param window \[`integer(1)`: `50L`\]\cr
#'   Rolling window size in observations.
#' @param align \[`character(1)`: `"right"`\]\cr
#'   Window alignment. The available options are:
#'   * `"right"`: each output value uses the `window` observations
#'     ending at that time point.
#'   * `"center"`: each output value is centered on the time point.
#'   * `"left"`: each output value uses the `window` observations
#'     starting at that time point.
#' @param method \[`character(1)`: `"periodogram"`\]\cr
#'   Spectral estimation method. The available options are:
#'   * `"periodogram"`: raw periodogram via [stats::spectrum()].
#'   * `"ar"`: autoregressive spectral estimate via [stats::spectrum()].
#' @param detrend \[`character(1)`: `"linear"`\]\cr
#'   Detrending applied to each window segment before spectral estimation.
#'   The available options are:
#'   * `"none"`: no detrending.
#'   * `"linear"`: remove a least-squares linear trend.
#'   * `"diff"`: first-difference the segment.
#' @param states \[`logical(1)`: `TRUE`\]\cr
#'   If `TRUE`, classify each time point into a noise-colour state based
#'   on the spectral exponent. If `FALSE`, return the tibble without the
#'   `state` column.
#' @param min_points \[`integer(1)`: `20L`\]\cr
#'   Minimum number of non-NA observations required in a window for
#'   spectral estimation.
#' @return An object of class `"spectral"` (inheriting from
#'   [tibble::tibble()]) with the following columns:
#'   \itemize{
#'     \item `time` -- integer or Date time index.
#'     \item `value` -- original series values.
#'     \item `spectral_exponent` -- rolling spectral exponent (beta).
#'     \item `spectral_ratio` -- rolling low-to-high frequency power ratio.
#'     \item `r_squared` -- R-squared of the log-log fit at each window.
#'     \item `state` -- (when `states = TRUE`) factor with levels
#'       `"white_noise"`, `"pink_noise"`, `"red_noise"`, `"brownian"`.
#'   }
#'   Attributes: `window`, `align`, `method`, `detrend`.
#'
#' @references
#' Held, H. & Kleinen, T. (2004). Detection of climate system bifurcations
#' by degenerate fingerprinting. \emph{Geophysical Research Letters}, 31,
#' L23207. \doi{10.1029/2004GL020972}
#'
#' Kleinen, T., Held, H., & Petschel-Held, G. (2003). The potential role
#' of spectral properties in detecting thresholds in the Earth system:
#' application to the thermohaline circulation. \emph{Ocean Dynamics}, 53,
#' 53--63. \doi{10.1007/s10236-002-0023-6}
#'
#' Dakos, V., Scheffer, M., van Nes, E.H., Brovkin, V., Petoukhov, V.,
#' & Held, H. (2008). Slowing down as an early warning signal for abrupt
#' climate change. \emph{Proceedings of the National Academy of Sciences},
#' 105(38), 14308--14312. \doi{10.1073/pnas.0802430105}
#'
#' @seealso [hurst()] for long-range dependence analysis;
#'   [compute_trend()] for trend classification;
#'   [detect_multivariate_warnings()] for multivariate EWS.
#' @family spectral
#' @concept time series
#' @concept early warning signals
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' sp <- spectral_ews(x, window = 50L)
#' plot(sp)
#' summary(sp)
#' }
spectral_ews <- function(data,
                         window = 50L,
                         align = "right",
                         method = "periodogram",
                         detrend = "linear",
                         states = TRUE,
                         min_points = 20L) {
  # --- 1. Input validation ---
  check_missing(data)
  data <- prepare_timeseries_data(data)
  values <- data$values
  time <- data$time
  n <- length(values)
  method <- check_match(method, c("periodogram", "ar"))
  detrend <- check_match(detrend, c("none", "linear", "diff"))
  align <- check_match(align, c("right", "center", "left"))
  check_range(window, type = "integer", min = 2L, max = n)
  check_range(min_points, type = "integer", min = 4L, max = window)
  check_flag(states)

  # --- 2. Rolling window computation ---
  beta_vec <- rep(NA_real_, n)
  ratio_vec <- rep(NA_real_, n)
  r2_vec <- rep(NA_real_, n)

  half_left <- (window - 1L) %/% 2L
  half_right <- window - 1L - half_left

  iter_start <- 1L
  iter_end <- n
  if (align == "center") {
    iter_start <- 1L + half_left
    iter_end <- n - half_right
  } else if (align == "right") {
    iter_start <- window
  } else {
    iter_end <- n - window + 1L
  }

  for (k in iter_start:iter_end) {
    w_idx <- if (align == "center") {
      (k - half_left):(k + half_right)
    } else if (align == "right") {
      (k - window + 1L):k
    } else {
      k:(k + window - 1L)
    }
    seg <- values[w_idx]
    n_valid <- sum(!is.na(seg))
    if (n_valid < min_points) next

    # Remove NAs for spectral computation
    seg_clean <- seg[!is.na(seg)]

    # Detrend
    seg_detrended <- spectral_detrend_(seg_clean, detrend)
    if (length(seg_detrended) < 10L) next

    # Compute spectral metrics
    result <- spectral_metrics_(seg_detrended, method)
    if (!is.null(result)) {
      beta_vec[k] <- result$beta
      ratio_vec[k] <- result$ratio
      r2_vec[k] <- result$r_squared
    }
  }

  # --- 3. Edge padding: three-step interpolation ---
  first_valid <- which(!is.na(beta_vec))[1L]
  last_valid <- rev(which(!is.na(beta_vec)))[1L]

  if (!is.na(first_valid) && !is.na(last_valid)) {
    # Forward-fill beginning
    if (first_valid > 1L) {
      beta_vec[seq_len(first_valid - 1L)] <- beta_vec[first_valid]
      ratio_vec[seq_len(first_valid - 1L)] <- ratio_vec[first_valid]
      r2_vec[seq_len(first_valid - 1L)] <- r2_vec[first_valid]
    }
    # Backward-fill end
    if (last_valid < n) {
      beta_vec[(last_valid + 1L):n] <- beta_vec[last_valid]
      ratio_vec[(last_valid + 1L):n] <- ratio_vec[last_valid]
      r2_vec[(last_valid + 1L):n] <- r2_vec[last_valid]
    }
    # Linear interpolation for internal gaps
    na_pos <- which(is.na(beta_vec))
    internal_na <- na_pos[na_pos > first_valid & na_pos < last_valid]
    if (length(internal_na) > 0L) {
      valid_pos <- which(!is.na(beta_vec))
      beta_vec <- stats::approx(
        valid_pos, beta_vec[valid_pos], xout = seq_len(n), rule = 2
      )$y
      ratio_vec <- stats::approx(
        valid_pos, ratio_vec[valid_pos], xout = seq_len(n), rule = 2
      )$y
      r2_vec <- stats::approx(
        valid_pos, r2_vec[valid_pos], xout = seq_len(n), rule = 2
      )$y
    }
  }

  # --- 4. State classification ---
  out <- tibble::tibble(
    time = time,
    value = data$values,
    spectral_exponent = beta_vec,
    spectral_ratio = ratio_vec,
    r_squared = r2_vec
  )
  if (states) {
    state_labels <- spectral_classify_(beta_vec)
    out$state <- factor(
      state_labels,
      levels = c("white_noise", "pink_noise", "red_noise", "brownian")
    )
  }

  structure(
    out,
    window = window,
    align = align,
    method = method,
    detrend = detrend,
    class = c("spectral", "tbl_df", "tbl", "data.frame")
  )
}

# --------------------------------------------------------------------------
# S3 methods
# --------------------------------------------------------------------------

#' Plot Spectral EWS Results
#'
#' Visualizes rolling spectral early warning signal results with
#' state-coloured backgrounds and spectral exponent trajectories.
#' Supports displaying the original time series, the spectral exponent
#' with threshold bands, or both panels stacked.
#'
#' @details
#' Three plot types are available:
#'
#' **`"series"`.**
#' Draws the original time series as a line plot with the background
#' shaded according to the noise-colour classification at each time
#' point. The four states (`white_noise`, `pink_noise`, `red_noise`,
#' `brownian`) are mapped to a fixed colour palette. This view answers:
#' "What was the system doing when the spectral character changed?"
#'
#' **`"states"`.**
#' Plots the rolling spectral exponent (beta) as a line, overlaid on
#' horizontal bands at the classification thresholds (0.5, 1.5, 2.5).
#' A dashed reference line at \eqn{\beta = 1.0} marks the boundary
#' between white/pink and red noise. This view answers: "How did the
#' spectral exponent evolve, and how close is it to regime boundaries?"
#'
#' **`"both"`** (default).
#' Stacks the series panel on top and the states panel below using
#' [patchwork::wrap_plots()], with shared time axes. The top panel
#' suppresses its x-axis labels to avoid duplication.
#'
#' @export
#' @param x \[`spectral`\]\cr
#'   An object of class `spectral`.
#' @param type \[`character(1)`: `"both"`\]\cr
#'   Plot type. The available options are:
#'   * `"series"`: Time series with state-coloured background.
#'   * `"states"`: Spectral exponent trajectory with threshold bands.
#'   * `"both"`: Both panels stacked vertically.
#' @param ... Additional arguments (currently unused).
#' @return A [ggplot2::ggplot()] object (or a
#'   [patchwork::wrap_plots()] composite when `type = "both"`).
#'
#' @references
#' Dakos, V., Scheffer, M., van Nes, E.H., Brovkin, V., Petoukhov, V.,
#' & Held, H. (2008). Slowing down as an early warning signal for abrupt
#' climate change. \emph{Proceedings of the National Academy of Sciences},
#' 105(38), 14308--14312. \doi{10.1073/pnas.0802430105}
#'
#' @seealso [spectral_ews()] for computing the rolling spectral metrics.
#' @family spectral
#' @concept time series
#' @concept early warning signals
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' sp <- spectral_ews(x, window = 50L)
#' plot(sp, type = "series")
#' plot(sp, type = "states")
#' plot(sp, type = "both")
#' }
plot.spectral <- function(x, type = "both", ...) {
  check_missing(x)
  check_class(x, "spectral")
  type <- check_match(type, c("series", "states", "both"))

  colors <- spectral_state_colors_()
  all_states <- c("white_noise", "pink_noise", "red_noise", "brownian")

  # Ensure state_f column for plotting
  if ("state" %in% names(x)) {
    d <- dplyr::mutate(
      x,
      state_f = factor(
        !!rlang::sym("state"),
        levels = all_states
      )
    )
  } else {
    # Derive states from spectral_exponent if not present
    d <- dplyr::mutate(
      x,
      state_f = factor(
        spectral_classify_(!!rlang::sym("spectral_exponent")),
        levels = all_states
      )
    )
  }

  if (type == "series" || type == "both") {
    rects <- spectral_build_rects_(d)
    val <- d[["value"]]
    val <- val[is.finite(val)]
    y_rng <- range(val, na.rm = TRUE)
    y_pad <- diff(y_rng) * 0.08
    y_lo <- y_rng[1L] - y_pad
    y_hi <- y_rng[2L] + y_pad
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
        name = "Spectral State",
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
    beta_vals <- d[["spectral_exponent"]]
    beta_vals <- beta_vals[is.finite(beta_vals)]
    beta_max <- max(beta_vals, 3.0, na.rm = TRUE)
    beta_min <- min(beta_vals, -0.5, na.rm = TRUE)
    beta_top <- beta_max + 0.1
    beta_bot <- beta_min - 0.1

    p_states <- ggplot2::ggplot(
      d,
      ggplot2::aes(
        x = !!rlang::sym("time"),
        y = !!rlang::sym("spectral_exponent")
      )
    ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = beta_bot, ymax = 0.5,
        fill = colors["white_noise"], alpha = 0.15
      ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 0.5, ymax = 1.5,
        fill = colors["pink_noise"], alpha = 0.15
      ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 1.5, ymax = 2.5,
        fill = colors["red_noise"], alpha = 0.15
      ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 2.5, ymax = beta_top,
        fill = colors["brownian"], alpha = 0.15
      ) +
      ggplot2::geom_line(linewidth = 0.5) +
      ggplot2::geom_hline(
        yintercept = 1.0,
        linetype = "dashed",
        color = "grey40"
      ) +
      ggplot2::scale_y_continuous(
        breaks = seq(
          floor(beta_bot * 2) / 2,
          ceiling(beta_top * 2) / 2,
          0.5
        )
      ) +
      ggplot2::coord_cartesian(
        ylim = c(beta_bot, beta_top),
        clip = "off"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Time", y = "Spectral Exponent (beta)") +
      ggplot2::theme(
        plot.margin = ggplot2::margin(t = 8, r = 5, b = 5, l = 5)
      )
    if (type == "states") return(p_states)
  }

  # type == "both"
  p_series <- p_series +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  patchwork::wrap_plots(p_series, p_states, ncol = 1L, heights = c(1, 1))
}

#' @describeIn spectral_ews Print method for `"spectral"` objects. Dispatches
#'   to the tibble print method.
#'
#' @export
#' @param x \[`spectral`\]\cr
#'   A `spectral` object.
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
print.spectral <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' @describeIn spectral_ews Summary method for `"spectral"` objects. Shows
#'   the state distribution, mean spectral exponent, and analysis settings.
#'
#' @export
#' @param object \[`spectral`\]\cr
#'   A `spectral` object.
#' @param ... Additional arguments (currently unused).
#' @return A list with elements `mean_beta`, `mean_ratio`, `mean_r_squared`,
#'   `state_counts` (if states present), `window`, `method`, `detrend`, and
#'   `align`, returned invisibly.
summary.spectral <- function(object, ...) {
  check_missing(object)
  check_class(object, "spectral")

  mean_beta <- mean(object$spectral_exponent, na.rm = TRUE)
  mean_ratio <- mean(object$spectral_ratio, na.rm = TRUE)
  mean_r2 <- mean(object$r_squared, na.rm = TRUE)

  cat("Spectral Early Warning Signal Summary\n")
  cat(strrep("=", 42L), "\n")
  cat("  Method  :", attr(object, "method"), "\n")
  cat("  Detrend :", attr(object, "detrend"), "\n")
  cat("  Window  :", attr(object, "window"), "\n")
  cat("  Align   :", attr(object, "align"), "\n")
  cat("  N       :", nrow(object), "\n\n")

  cat("  Mean spectral exponent (beta):", round(mean_beta, 4L), "\n")
  cat("  Mean spectral ratio          :", round(mean_ratio, 4L), "\n")
  cat("  Mean R-squared               :", round(mean_r2, 4L), "\n")

  state_counts <- NULL
  if ("state" %in% names(object)) {
    state_counts <- table(object$state)
    props <- prop.table(state_counts)
    cat("\n  State distribution:\n")
    for (i in seq_along(state_counts)) {
      if (state_counts[i] > 0L) {
        cat(sprintf("    %-12s %5d  (%5.1f%%)\n",
                    names(state_counts)[i], state_counts[i],
                    props[i] * 100))
      }
    }
  }

  invisible(list(
    mean_beta = mean_beta,
    mean_ratio = mean_ratio,
    mean_r_squared = mean_r2,
    state_counts = state_counts,
    window = attr(object, "window"),
    method = attr(object, "method"),
    detrend = attr(object, "detrend"),
    align = attr(object, "align")
  ))
}
