# ============================================================================
# sensitivity.R -- Parameter sensitivity analysis for early warning signals
# Evaluates how EWS metric trends (Kendall tau) vary across window sizes
# and detrending methods to assess robustness of detected signals.
# ============================================================================

# --------------------------------------------------------------------------
# Internal helpers (@noRd)
# --------------------------------------------------------------------------

#' Detrend a time series
#'
#' Removes a linear trend via OLS residuals. Returns the original series
#' when `method = "none"`.
#'
#' @param values Numeric vector.
#' @param method One of `"none"` or `"linear"`.
#' @return Numeric vector of same length.
#' @noRd
sensitivity_detrend_ <- function(values, method) {
  if (method == "none") return(values)
  n <- length(values)
  time_idx <- seq_len(n)
  valid <- !is.na(values)
  if (sum(valid) < 3L) return(values)
  fit <- try_(stats::lm(values[valid] ~ time_idx[valid]))
  if (inherits(fit, "try-error")) return(values)
  out <- values
  out[valid] <- stats::residuals(fit)
  out
}

#' Compute a single EWS metric on a numeric window
#'
#' @param x Numeric vector (window data, no NAs expected).
#' @param metric One of "ar1", "sd", "variance", "skewness", "kurtosis", "cv".
#' @return Single numeric value.
#' @noRd
sensitivity_metric_ <- function(x, metric) {
  n <- length(x)
  if (n < 3L) return(NA_real_)
  switch(
    metric,
    ar1 = {
      if (stats::var(x, na.rm = TRUE) < 1e-10) return(NA_real_)
      model <- try_(
        stats::ar.ols(x, aic = FALSE, order.max = 1L,
                      demean = FALSE, intercept = TRUE)
      )
      if (inherits(model, "try-error") || is.null(model) ||
          length(model$ar) == 0L) {
        return(NA_real_)
      }
      max(-0.999, min(0.999, model$ar[1L]))
    },
    sd = stats::sd(x, na.rm = TRUE),
    variance = stats::var(x, na.rm = TRUE),
    skewness = {
      mu <- mean(x, na.rm = TRUE)
      s <- stats::sd(x, na.rm = TRUE)
      if (is.na(s) || s < 1e-10) return(NA_real_)
      mean(((x - mu) / s)^3, na.rm = TRUE)
    },
    kurtosis = {
      mu <- mean(x, na.rm = TRUE)
      s <- stats::sd(x, na.rm = TRUE)
      if (is.na(s) || s < 1e-10) return(NA_real_)
      mean(((x - mu) / s)^4, na.rm = TRUE) - 3
    },
    cv = {
      mu <- mean(x, na.rm = TRUE)
      if (abs(mu) < .Machine$double.eps) return(NA_real_)
      stats::sd(x, na.rm = TRUE) / abs(mu)
    }
  )
}

#' Compute rolling metric values for a detrended series
#'
#' @param values Numeric vector (already detrended).
#' @param window Integer window size.
#' @param metric Character metric name.
#' @return List with `time` (integer indices) and `score` (metric values).
#' @noRd
sensitivity_rolling_ <- function(values, window, metric) {
  n <- length(values)
  n_windows <- n - window + 1L
  if (n_windows < 1L) {
    return(list(
      time = integer(0),
      score = numeric(0)
    ))
  }
  scores <- vapply(seq_len(n_windows), function(i) {
    wd <- values[i:(i + window - 1L)]
    clean <- wd[!is.na(wd)]
    if (length(clean) < 3L) return(NA_real_)
    sensitivity_metric_(clean, metric)
  }, numeric(1L))
  # Right-aligned: metric at position window, window+1, ..., n
  list(
    time = seq(window, n),
    score = scores
  )
}

#' Compute Kendall tau of a metric series against time
#'
#' @param scores Numeric vector of metric values.
#' @return Single numeric Kendall tau value, or `NA_real_`.
#' @noRd
sensitivity_tau_ <- function(scores) {
  valid <- !is.na(scores)
  if (sum(valid) < 4L) return(NA_real_)
  time_idx <- seq_along(scores)[valid]
  vals <- scores[valid]
  result <- try_(stats::cor.test(time_idx, vals, method = "kendall"))
  if (inherits(result, "try-error")) return(NA_real_)
  result$estimate
}

#' Color palette for sensitivity heatmap
#' @noRd
sensitivity_colors_ <- function() {
  c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
    "#F7F7F7",
    "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
}

# --------------------------------------------------------------------------
# Exported: sensitivity_ews()
# --------------------------------------------------------------------------

#' Parameter Sensitivity Analysis for Early Warning Signals
#'
#' Evaluates how the Kendall tau trend statistic of an EWS metric varies
#' across different rolling window sizes and detrending methods. This
#' reveals whether a detected signal is robust to parameter choices or an
#' artifact of a specific configuration.
#'
#' @details
#' The analysis proceeds in four steps:
#'
#' **1. Data preparation.** The input time series is parsed via
#' [prepare_timeseries_data()]. If `windows` is `NULL`, a sequence of 8
#' window sizes is generated spanning from `max(20, n/20)` to
#' `min(n/2, 200)` on a linear grid, rounded to integers. This range
#' covers small windows (high temporal resolution but noisy estimates)
#' through large windows (smooth estimates but low resolution).
#'
#' **2. Detrending.** For each detrend method, the full series is
#' preprocessed: `"none"` uses the raw series; `"linear"` removes the
#' best-fit linear trend via OLS residuals. Detrending before computing
#' EWS metrics eliminates confounding from deterministic trends that
#' could inflate autocorrelation or variance artificially.
#'
#' **3. Rolling metric computation.** For each (window, detrend)
#' combination, the specified metric is computed on right-aligned rolling
#' windows. The available metrics are:
#' \describe{
#'   \item{`"ar1"`}{Lag-1 autocorrelation estimated via OLS AR(1).
#'     Rising AR(1) is a classic signature of critical slowing down
#'     (Scheffer et al., 2009).}
#'   \item{`"sd"`}{Standard deviation. Increasing variance is the
#'     second canonical EWS indicator alongside AR(1).}
#'   \item{`"variance"`}{Variance (SD squared). Equivalent information
#'     to SD but on the original scale, which may be preferable for
#'     some comparisons.}
#'   \item{`"skewness"`}{Third standardized moment. Non-zero skewness
#'     indicates asymmetric fluctuations, which can precede transitions
#'     in systems with asymmetric potential wells.}
#'   \item{`"kurtosis"`}{Excess kurtosis (fourth standardized moment
#'     minus 3). Heavy-tailed fluctuations (positive kurtosis) can
#'     signal proximity to a transition.}
#'   \item{`"cv"`}{Coefficient of variation (SD / |mean|). Useful when
#'     the mean is non-zero and relative dispersion matters more than
#'     absolute dispersion.}
#' }
#'
#' **4. Kendall tau trend test.** For each (window, detrend)
#' combination, Kendall's rank correlation between the metric values and
#' their time indices is computed. A consistently positive tau across
#' many parameter combinations provides robust evidence for an upward
#' trend in the EWS metric (and thus for approaching a critical
#' transition). Sensitivity is assessed by examining whether tau remains
#' consistently positive (or negative) regardless of window size and
#' detrending.
#'
#' @export
#' @param data \[`numeric` or `ts`\]\cr
#'   The time series to analyze. Accepts a numeric vector or a `ts` object.
#' @param metric \[`character(1)`: `"ar1"`\]\cr
#'   EWS metric to test. One of `"ar1"`, `"sd"`, `"variance"`,
#'   `"skewness"`, `"kurtosis"`, `"cv"`.
#' @param windows \[`integer()`: `NULL`\]\cr
#'   Vector of window sizes to evaluate. If `NULL`, a default sequence
#'   of 8 sizes is generated spanning `max(20, n/20)` to
#'   `min(n/2, 200)`.
#' @param detrend_methods \[`character()`: `c("none", "linear")`\]\cr
#'   Detrending methods to compare. `"none"` uses the raw series;
#'   `"linear"` removes a linear trend via OLS residuals.
#' @param method \[`character(1)`: `"rolling"`\]\cr
#'   Analysis method. Currently only `"rolling"` is supported.
#' @return An object of class `"sensitivity_ews"` (inheriting from
#'   [tibble::tbl_df]) with the following columns:
#'   \itemize{
#'     \item `window` -- integer window size used.
#'     \item `detrend` -- character detrend method applied.
#'     \item `time` -- time index for each metric value.
#'     \item `score` -- computed metric value at that time point.
#'     \item `tau` -- Kendall's tau for this (window, detrend)
#'       combination.
#'   }
#'
#'   Attributes stored on the object:
#'   \itemize{
#'     \item `metric` -- the EWS metric name.
#'     \item `windows` -- integer vector of window sizes evaluated.
#'     \item `detrend_methods` -- character vector of detrend methods.
#'   }
#'
#' @references
#' Scheffer, M., Bascompte, J., Brock, W.A., Brovkin, V., Carpenter,
#' S.R., Dakos, V., Held, H., van Nes, E.H., Rietkerk, M., &
#' Sugihara, G. (2009). Early-warning signals for critical transitions.
#' \emph{Nature}, 461, 53--59. \doi{10.1038/nature08227}
#'
#' Dakos, V., Carpenter, S.R., Brock, W.A., Ellison, A.M., Guttal, V.,
#' Ives, A.R., Kefi, S., Livina, V., Seekell, D.A., van Nes, E.H., &
#' Scheffer, M. (2012). Methods for detecting early warnings of
#' critical transitions in time series illustrated using simulated
#' ecological data. \emph{PLoS ONE}, 7(7), e41010.
#' \doi{10.1371/journal.pone.0041010}
#'
#' Lenton, T.M. (2011). Early warning of climate tipping points.
#' \emph{Nature Climate Change}, 1(4), 201--209.
#' \doi{10.1038/nclimate1143}
#'
#' @seealso [resilience()] for rolling resilience metric computation;
#'   [compute_trend()] for trend classification; [hurst()] for
#'   Hurst exponent analysis.
#' @family sensitivity
#' @concept early warning signals
#' @concept sensitivity analysis
#' @examples
#' \donttest{
#' set.seed(42)
#' # Generate a series with increasing variance (approaching transition)
#' x <- cumsum(rnorm(300, sd = seq(0.5, 2, length.out = 300)))
#' sa <- sensitivity_ews(x, metric = "ar1")
#' sa
#' summary(sa)
#' plot(sa, type = "heatmap")
#' plot(sa, type = "lines")
#' }
sensitivity_ews <- function(data,
                            metric = "ar1",
                            windows = NULL,
                            detrend_methods = c("none", "linear"),
                            method = "rolling") {

  # --- 1. Input validation ---
  check_missing(data)
  data <- prepare_timeseries_data(data)
  values <- data$values
  time <- data$time
  n <- length(values)

  metric <- check_match(
    metric, c("ar1", "sd", "variance", "skewness", "kurtosis", "cv")
  )
  detrend_methods <- check_match(
    detrend_methods, c("none", "linear"), several.ok = TRUE
  )
  method <- check_match(method, c("rolling"))

  # --- 2. Generate default window sequence if NULL ---
  if (is.null(windows)) {
    w_min <- max(20L, round(n / 20))
    w_max <- min(round(n / 2), 200L)
    if (w_min >= w_max) {
      windows <- unique(as.integer(round(
        seq(w_min, w_max, length.out = 2L)
      )))
    } else {
      windows <- unique(as.integer(round(
        seq(w_min, w_max, length.out = 8L)
      )))
    }
  } else {
    windows <- as.integer(round(windows))
  }
  # Validate each window size
  windows <- windows[windows >= 2L & windows <= n]
  if (length(windows) == 0L) {
    stop_(
      "No valid window sizes. All values must be between 2 and {n}."
    )
  }

  # --- 3. Compute metrics for all (window, detrend) combinations ---
  results <- vector("list", length(windows) * length(detrend_methods))
  idx <- 0L

  for (dm in detrend_methods) {
    detrended <- sensitivity_detrend_(values, dm)
    for (w in windows) {
      idx <- idx + 1L
      rolling <- sensitivity_rolling_(detrended, w, metric)
      if (length(rolling$score) == 0L) {
        results[[idx]] <- tibble::tibble(
          window  = integer(0),
          detrend = character(0),
          time    = numeric(0),
          score   = numeric(0),
          tau     = numeric(0)
        )
        next
      }
      tau_val <- sensitivity_tau_(rolling$score)
      results[[idx]] <- tibble::tibble(
        window  = rep(as.integer(w), length(rolling$score)),
        detrend = rep(dm, length(rolling$score)),
        time    = time[rolling$time],
        score   = rolling$score,
        tau     = rep(tau_val, length(rolling$score))
      )
    }
  }

  out <- do.call(rbind, results)
  out <- tibble::as_tibble(out)

  structure(
    out,
    metric          = metric,
    windows         = windows,
    detrend_methods = detrend_methods,
    class           = c("sensitivity_ews", "tbl_df", "tbl", "data.frame")
  )
}

# --------------------------------------------------------------------------
# S3 methods
# --------------------------------------------------------------------------

#' Plot Sensitivity Analysis Results
#'
#' Visualizes how the Kendall tau trend statistic and metric trajectories
#' vary across parameter combinations (window size and detrending method).
#'
#' @details
#' Three plot types are available:
#'
#' **`type = "heatmap"`.** Displays a tile plot with window size on the
#' y-axis and detrend method on the x-axis. Each tile is colored by the
#' Kendall tau value using a diverging blue-white-red scale (blue =
#' negative tau, white = zero, red = positive tau). Consistently red
#' tiles across all combinations indicate a robust upward trend in the
#' EWS metric, providing confidence that the signal is not an artifact
#' of a specific parameter choice.
#'
#' **`type = "lines"`.** Faceted line plots showing the metric trajectory
#' over time for each (window, detrend) combination. Panels are faceted
#' by window size and detrend method, with free y-axes. This view reveals
#' how the raw metric evolution differs across parameter settings and
#' helps identify whether apparent trends are driven by specific window
#' sizes or detrending choices.
#'
#' **`type = "both"`.** Stacks the heatmap on top and the line plots
#' below using [patchwork::wrap_plots()], providing a complete overview
#' in a single figure.
#'
#' @export
#' @param x \[`sensitivity_ews`\]\cr
#'   An object of class `sensitivity_ews` as returned by
#'   [sensitivity_ews()].
#' @param type \[`character(1)`: `"heatmap"`\]\cr
#'   Plot type. The available options are:
#'   * `"heatmap"`: Tile plot of Kendall tau across parameter
#'     combinations.
#'   * `"lines"`: Faceted line plots of metric trajectories.
#'   * `"both"`: Both panels stacked vertically.
#' @param ... Additional arguments (currently unused).
#' @return A [ggplot2::ggplot()] object (or a
#'   [patchwork::wrap_plots()] composite when `type = "both"`).
#'
#' @seealso [sensitivity_ews()] for computing the sensitivity analysis.
#' @family sensitivity
#' @concept early warning signals
#' @concept sensitivity analysis
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(300, sd = seq(0.5, 2, length.out = 300)))
#' sa <- sensitivity_ews(x, metric = "ar1")
#' plot(sa, type = "heatmap")
#' plot(sa, type = "lines")
#' plot(sa, type = "both")
#' }
plot.sensitivity_ews <- function(x, type = "heatmap", ...) {
  check_missing(x)
  check_class(x, "sensitivity_ews")
  type <- check_match(type, c("heatmap", "lines", "both"))

  metric_name <- attr(x, "metric") %||% "metric"

  if (type == "heatmap" || type == "both") {
    p_heat <- sensitivity_plot_heatmap_(x, metric_name)
  }
  if (type == "lines" || type == "both") {
    p_lines <- sensitivity_plot_lines_(x, metric_name)
  }

  if (type == "heatmap") return(p_heat)
  if (type == "lines") return(p_lines)
  patchwork::wrap_plots(p_heat, p_lines, ncol = 1L, heights = c(1, 2))
}

#' @noRd
sensitivity_plot_heatmap_ <- function(x, metric_name) {
  # Extract one row per (window, detrend) combination
  tau_df <- unique(x[, c("window", "detrend", "tau")])
  tau_df <- dplyr::mutate(
    tau_df,
    window_label = factor(
      !!rlang::sym("window"),
      levels = sort(unique(!!rlang::sym("window")))
    )
  )

  tau_range <- range(tau_df$tau, na.rm = TRUE)
  abs_max <- max(abs(tau_range), na.rm = TRUE)
  if (abs_max < 0.01) abs_max <- 1

  ggplot2::ggplot(
    tau_df,
    ggplot2::aes(
      x    = !!rlang::sym("detrend"),
      y    = !!rlang::sym("window_label"),
      fill = !!rlang::sym("tau")
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", !!rlang::sym("tau"))),
      size = 3.5,
      color = "black"
    ) +
    ggplot2::scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "#F7F7F7",
      high     = "#B2182B",
      midpoint = 0,
      limits   = c(-abs_max, abs_max),
      name     = "Kendall tau"
    ) +
    ggplot2::labs(
      title = sprintf("Sensitivity of %s: Kendall tau", toupper(metric_name)),
      x     = "Detrend Method",
      y     = "Window Size"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.title = ggplot2::element_text(color = "black", face = "bold"),
      axis.text  = ggplot2::element_text(color = "black"),
      panel.grid = ggplot2::element_blank()
    )
}

#' @noRd
sensitivity_plot_lines_ <- function(x, metric_name) {
  plot_df <- dplyr::mutate(
    x,
    facet_label = sprintf("window = %d, detrend = %s",
                          !!rlang::sym("window"),
                          !!rlang::sym("detrend"))
  )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = !!rlang::sym("time"),
      y = !!rlang::sym("score")
    )
  ) +
    ggplot2::geom_line(linewidth = 0.4, color = "#2166AC") +
    ggplot2::facet_wrap(
      ggplot2::vars(!!rlang::sym("facet_label")),
      scales = "free_y"
    ) +
    ggplot2::labs(
      title = sprintf("Rolling %s Across Parameter Combinations",
                      toupper(metric_name)),
      x     = "Time",
      y     = metric_name
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 14, face = "bold"),
      axis.title      = ggplot2::element_text(color = "black", face = "bold"),
      axis.text       = ggplot2::element_text(color = "black"),
      strip.text      = ggplot2::element_text(face = "bold", size = 8),
      legend.position = "none"
    )
}

#' @describeIn sensitivity_ews Print method for `"sensitivity_ews"` objects.
#'   Dispatches to the tibble print method.
#'
#' @export
#' @param x \[`sensitivity_ews`\]\cr
#'   A `sensitivity_ews` object.
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
print.sensitivity_ews <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' @describeIn sensitivity_ews Summary method for `"sensitivity_ews"` objects.
#'   Prints the range of Kendall tau values, identifies the most and least
#'   robust parameter combinations, and assesses overall consistency.
#'
#' @export
#' @param object \[`sensitivity_ews`\]\cr
#'   A `sensitivity_ews` object.
#' @param ... Additional arguments (currently unused).
#' @return A list with elements `metric`, `tau_range`, `tau_mean`,
#'   `n_positive`, `n_negative`, `n_total`, `most_robust` (combination
#'   with highest |tau|), `least_robust` (combination with lowest |tau|),
#'   and `consistent` (logical), returned invisibly.
summary.sensitivity_ews <- function(object, ...) {
  check_missing(object)
  check_class(object, "sensitivity_ews")

  metric_name <- attr(object, "metric") %||% "unknown"
  tau_df <- unique(object[, c("window", "detrend", "tau")])
  tau_vals <- tau_df$tau[!is.na(tau_df$tau)]

  n_total <- length(tau_vals)
  n_positive <- sum(tau_vals > 0)
  n_negative <- sum(tau_vals < 0)
  tau_mean <- mean(tau_vals, na.rm = TRUE)
  tau_range <- range(tau_vals, na.rm = TRUE)

  # Most robust: highest absolute tau
  abs_tau <- abs(tau_df$tau)
  most_idx <- which.max(abs_tau)
  least_idx <- which.min(abs_tau)
  most_robust <- tau_df[most_idx, ]
  least_robust <- tau_df[least_idx, ]

  # Consistency: all tau values have the same sign
  consistent <- n_total > 0L && (n_positive == n_total || n_negative == n_total)

  cat("Sensitivity Analysis Summary\n")
  cat("  Metric        :", metric_name, "\n")
  cat("  Combinations  :", n_total, "\n")
  cat("  Tau range     : [", sprintf("%.4f", tau_range[1L]),
      ",", sprintf("%.4f", tau_range[2L]), "]\n")
  cat("  Tau mean      :", sprintf("%.4f", tau_mean), "\n")
  cat("  Positive tau  :", n_positive, "/", n_total, "\n")
  cat("  Negative tau  :", n_negative, "/", n_total, "\n\n")

  cat("  Most robust   : window =", most_robust$window,
      ", detrend =", most_robust$detrend,
      ", tau =", sprintf("%.4f", most_robust$tau), "\n")
  cat("  Least robust  : window =", least_robust$window,
      ", detrend =", least_robust$detrend,
      ", tau =", sprintf("%.4f", least_robust$tau), "\n\n")

  if (consistent) {
    direction <- ifelse_(tau_mean > 0, "positive", "negative")
    cat("  Consistency   : ALL tau values are", direction, "\n")
    cat("  Assessment    : Signal is ROBUST across parameter choices.\n")
  } else {
    cat("  Consistency   : Mixed signs detected\n")
    cat("  Assessment    : Signal DEPENDS on parameter choices;",
        "interpret with caution.\n")
  }

  invisible(list(
    metric       = metric_name,
    tau_range    = tau_range,
    tau_mean     = tau_mean,
    n_positive   = n_positive,
    n_negative   = n_negative,
    n_total      = n_total,
    most_robust  = most_robust,
    least_robust = least_robust,
    consistent   = consistent
  ))
}
