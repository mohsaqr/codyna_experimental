# ============================================================================
# trend.R --Rolling trend classification with turbulence detection
# Classifies local time series behavior as ascending, descending, flat,
# or turbulent using rolling-window trend metrics and volatility analysis.
# ============================================================================

# --------------------------------------------------------------------------
# Internal helpers (@noRd)
# --------------------------------------------------------------------------

#' AR(1) coefficient for a single window
#'
#' Computes the first-order autoregressive coefficient via OLS, bounded to
#' (-1, 1). Returns `NA_real_` when data is insufficient or variance is
#' near zero.
#'
#' @param window_data Numeric vector.
#' @param min_points Minimum non-NA points required.
#' @return Single numeric, or `NA_real_`.
#' @noRd
trend_ar1_ <- function(window_data, min_points = 3L) {
  clean <- stats::na.omit(as.numeric(window_data))
  if (length(clean) < min_points) return(NA_real_)
  if (stats::var(clean, na.rm = TRUE) < 1e-10) return(NA_real_)
  model <- try_(suppressWarnings(
    stats::ar.ols(
      clean, aic = FALSE, order.max = 1L,
      demean = FALSE, intercept = TRUE
    )
  ))
  if (inherits(model, "try-error")) model <- NULL
  if (is.null(model) || length(model$ar) == 0L) return(NA_real_)
  max(-0.999, min(0.999, model$ar[1L]))
}

#' Slope metric for a single window
#' @noRd
trend_slope_ <- function(x_vals, y_vals, slope_method, min_points) {
  valid <- !is.na(x_vals) & !is.na(y_vals)
  if (sum(valid) < min_points) return(NA_real_)
  xc <- x_vals[valid]
  yc <- y_vals[valid]
  if (length(unique(xc)) < 2L && slope_method != "robust") return(NA_real_)
  if (slope_method == "robust" && length(xc) < 2L) return(NA_real_)
  if (stats::var(xc, na.rm = TRUE) < 1e-10) {
    return(ifelse_(stats::var(yc, na.rm = TRUE) < 1e-10, 0, NA_real_))
  }
  if (stats::var(yc, na.rm = TRUE) < 1e-10) return(0)
  switch(
    slope_method,
    ols = stats::cov(xc, yc) / stats::var(xc),
    robust = {
      np <- length(xc)
      idx <- which(upper.tri(matrix(0, np, np)), arr.ind = TRUE)
      dx <- xc[idx[, 2L]] - xc[idx[, 1L]]
      dy <- yc[idx[, 2L]] - yc[idx[, 1L]]
      keep <- dx != 0
      if (!any(keep)) return(NA_real_)
      stats::median(dy[keep] / dx[keep], na.rm = TRUE)
    },
    spearman = ,
    kendall = {
      cor_val <- try_(stats::cor(xc, yc, method = slope_method,
                                    use = "complete.obs"))
      if (inherits(cor_val, "try-error")) cor_val <- NA_real_
      sd_x <- stats::sd(xc, na.rm = TRUE)
      if (is.na(cor_val) || is.na(sd_x) || sd_x < 1e-10) {
        NA_real_
      } else {
        cor_val * (stats::sd(yc, na.rm = TRUE) / sd_x)
      }
    }
  )
}

#' Turbulence combined volatility score
#' @noRd
trend_volatility_ <- function(seg) {
  m_sd   <- stats::sd(seg, na.rm = TRUE)
  m_mean <- abs(mean(seg, na.rm = TRUE))
  m_rng  <- diff(range(seg, na.rm = TRUE))
  if (is.na(m_sd) || is.na(m_mean) || is.na(m_rng)) return(NA_real_)
  m_sd / (m_mean + 1e-8) + 0.5 * m_rng / (m_mean + 1e-8)
}

#' Trend state colour palette
#' @noRd
trend_state_colors_ <- function() {
  c(
    ascending    = "#4CAF50",
    descending   = "#F44336",
    flat         = "#FFC107",
    turbulent    = "#2196F3",
    Missing_Data = "#CCCCCC",
    Initial      = "#FFFFFF"
  )
}

# --------------------------------------------------------------------------
# Exported: compute_trend()
# --------------------------------------------------------------------------

#' Compute and Classify Rolling Time Series Trends
#'
#' Performs rolling-window trend analysis on a time series, classifying each
#' time point as **ascending**, **descending**, **flat**, or **turbulent**.
#' Turbulence is detected from the volatility of the trend metric itself,
#' not from the raw series variance.
#'
#' @details
#' The analysis proceeds in three stages:
#'
#' **1. Rolling metric calculation.** A trend metric is computed for each
#' window position. The available methods are:
#' \describe{
#'   \item{`"slope"` (default)}{Rate of change estimated by one of four
#'     sub-methods:
#'     \itemize{
#'       \item `"ols"` --Ordinary Least Squares slope. Efficient but
#'         sensitive to outliers.
#'       \item `"robust"` --Theil-Sen estimator (median of all pairwise
#'         slopes). Highly robust to outliers.
#'       \item `"spearman"` / `"kendall"` --Rank-based correlation scaled
#'         by the ratio of standard deviations to produce a slope-like
#'         metric. Robust to non-linear relationships.
#'     }
#'   }
#'   \item{`"ar1_phi1"`}{First-order autoregressive coefficient. Values
#'     near 1 indicate a strong trend; near 0 indicate no trend.}
#'   \item{`"growth_factor"`}{Ratio of the last to the first value in each
#'     window (\eqn{V_{end} / V_{start}}). Suitable for multiplicative or
#'     exponential processes.}
#' }
#'
#' **2. Initial classification.** Each point is classified based on the
#' metric value relative to a neutral value (0 for slope, 1 for ar1_phi1
#' and growth_factor) plus or minus `epsilon`:
#' \itemize{
#'   \item **ascending**: metric > neutral + epsilon
#'   \item **descending**: metric < neutral - epsilon
#'   \item **flat**: metric within the epsilon band
#' }
#'
#' **3. Turbulence re-classification.** A combined volatility score is
#' computed on a rolling sub-window of the metric itself:
#' \deqn{V_{combined} = CV(metric) + 0.5 \times \frac{Range(metric)}{|mean(metric)|}}
#' If this score exceeds `turbulence_threshold`, the point is re-classified
#' as **turbulent**. The `flat_to_turbulent_factor` introduces hysteresis:
#' flat segments require a higher volatility score to be reclassified,
#' preventing minor noise from triggering false turbulence detections.
#'
#' @export
#' @param data \[`numeric` or `ts`\]\cr
#'   The time series to analyze. Accepts a numeric vector or a `ts` object.
#' @param window \[`integer(1)`: `NULL`\]\cr
#'   Rolling window size. If `NULL`, an adaptive size is used:
#'   `max(3, min(n, round(n / 10)))`.
#' @param method \[`character(1)`: `"slope"`\]\cr
#'   Trend quantification method. One of `"slope"`, `"ar1_phi1"`, or
#'   `"growth_factor"`.
#' @param slope_method \[`character(1)`: `"robust"`\]\cr
#'   Slope estimation method when `method = "slope"`. One of `"ols"`,
#'   `"robust"`, `"spearman"`, or `"kendall"`.
#' @param epsilon \[`numeric(1)`: `0.05`\]\cr
#'   Threshold for the flat classification band around the neutral value.
#' @param turbulence_threshold \[`numeric(1)`: `5`\]\cr
#'   Baseline threshold for the combined volatility score.
#' @param flat_to_turbulent_factor \[`numeric(1)`: `1.5`\]\cr
#'   Multiplier applied to `turbulence_threshold` for points already
#'   classified as flat. Must be >= 1.
#' @param min_points \[`integer(1)`: `3L`\]\cr
#'   Minimum non-NA observations required in a window.
#' @param align \[`character(1)`: `"center"`\]\cr
#'   Window alignment: `"center"`, `"right"`, or `"left"`.
#' @return An object of class `"trend"` (a tibble) with columns:
#'   \itemize{
#'     \item `time` --the time index.
#'     \item `value` --the input time series values.
#'     \item `metric` --the rolling trend metric.
#'     \item `state` --factor with levels `"ascending"`,
#'       `"descending"`, `"flat"`, `"turbulent"`, `"Missing_Data"`,
#'       `"Initial"`.
#'   }
#'
#' @seealso [compare_ts()] for visualizing trend classifications against
#'   the original series; [resilience()] for rolling resilience metrics.
#' @family trend
#' @concept time series
#' @examples
#' \donttest{
#' set.seed(123)
#' x <- cumsum(rnorm(200, sd = 2))
#' tr <- compute_trend(x, window = 20)
#' tr
#' table(tr$state)
#' plot(tr)
#' }
compute_trend <- function(data,
                          window = NULL,
                          method = "slope",
                          slope_method = "robust",
                          epsilon = 0.05,
                          turbulence_threshold = 5,
                          flat_to_turbulent_factor = 1.5,
                          min_points = 3L,
                          align = "center") {

  # --- 1. Input validation ---
  check_missing(data)
  data <- prepare_timeseries_data(data)
  values <- data$values
  time   <- data$time
  n      <- length(values)

  method <- check_match(method, c("slope", "ar1_phi1", "growth_factor"))
  if (method == "slope") {
    slope_method <- check_match(
      slope_method, c("ols", "robust", "spearman", "kendall")
    )
  }
  align <- check_match(align, c("center", "right", "left"))
  check_range(epsilon, type = "numeric", scalar = TRUE, min = 0, max = Inf)
  check_range(
    turbulence_threshold, type = "numeric", scalar = TRUE, min = 0, max = Inf
  )
  check_range(
    flat_to_turbulent_factor, type = "numeric", scalar = TRUE, min = 1, max = Inf
  )

  if (is.null(window)) {
    window <- max(3L, min(n, round(n / 10)))
  }
  check_range(window, type = "integer", min = 2L, max = n)
  if (min_points > window) min_points <- window

  # --- 2. Rolling metric calculation ---
  metric <- rep(NA_real_, n)
  half_left  <- (window - 1L) %/% 2L
  half_right <- window - 1L - half_left

  iter_start <- 1L
  iter_end   <- n
  if (align == "center") {
    iter_start <- 1L + half_left
    iter_end   <- n - half_right
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
    wd <- values[w_idx]
    if (sum(!is.na(wd)) < min_points) next

    metric[k] <- switch(
      method,
      slope = trend_slope_(seq_along(wd), wd, slope_method, min_points),
      ar1_phi1 = trend_ar1_(wd, max(2L, min_points)),
      growth_factor = {
        clean <- wd[!is.na(wd)]
        if (length(clean) >= 2L && abs(clean[1L]) > 1e-10) {
          clean[length(clean)] / clean[1L]
        } else {
          NA_real_
        }
      }
    )
  }

  # --- 3. Edge padding for center alignment ---
  if (align == "center" && any(!is.na(metric))) {
    first_valid <- min(which(!is.na(metric)))
    if (first_valid > 1L) {
      metric[seq_len(first_valid - 1L)] <- metric[first_valid]
    }
    last_valid <- max(which(!is.na(metric)))
    if (last_valid < n) {
      metric[(last_valid + 1L):n] <- metric[last_valid]
    }
  }

  # --- 4. Initial trend classification ---
  state <- rep("Initial", n)
  na_data <- is.na(values)
  state[na_data] <- "Missing_Data"
  neutral <- ifelse_(method %in% c("ar1_phi1", "growth_factor"), 1, 0)
  upper <- neutral + epsilon
  lower <- neutral - epsilon

  valid_idx <- !is.na(metric) & !na_data
  state[valid_idx] <- ifelse(
    metric[valid_idx] > upper, "ascending",
    ifelse(metric[valid_idx] < lower, "descending", "flat")
  )

  # --- 5. Turbulence re-classification ---
  if (sum(valid_idx) >= min_points) {
    vol_win <- min(max(3L, window %/% 2L), sum(valid_idx))
    metrics_valid <- metric[valid_idx]
    orig_indices  <- which(valid_idx)
    if (length(metrics_valid) >= vol_win) {
      for (j in vol_win:length(metrics_valid)) {
        seg <- metrics_valid[(j - vol_win + 1L):j]
        if (sum(!is.na(seg)) < max(2L, min_points - 1L)) next
        combined <- trend_volatility_(seg)
        if (is.na(combined)) next

        idx_orig   <- orig_indices[j]
        base_trend <- state[idx_orig]
        eff_thresh <- ifelse_(
          base_trend == "flat",
          turbulence_threshold * flat_to_turbulent_factor,
          turbulence_threshold
        )
        if (base_trend != "Missing_Data" && combined > eff_thresh) {
          state[idx_orig] <- "turbulent"
        }
      }
    }
  }

  # --- 6. Output ---
  state <- factor(
    state,
    levels = c("ascending", "descending", "flat",
               "turbulent", "Missing_Data", "Initial")
  )
  out <- tibble::tibble(
    time   = time,
    value  = values,
    metric = metric,
    state  = state
  )
  structure(
    out,
    window       = window,
    align        = align,
    method       = method,
    slope_method = if (method == "slope") slope_method else NULL,
    epsilon      = epsilon,
    class        = c("trend", "tbl_df", "tbl", "data.frame")
  )
}

# --------------------------------------------------------------------------
# S3 methods
# --------------------------------------------------------------------------

#' Plot Trend Classification Results
#'
#' Visualizes the time series with trend state classification, either as
#' background-shaded regions or as a ribbon bar beneath the series.
#'
#' @details
#' Two plot types are available:
#'
#' **`type = "series"` (default).** The time series is drawn with coloured
#' background tiles indicating the trend state at each time point. This
#' provides an immediate overlay of direction and trajectory.
#'
#' **`type = "ribbons"`.** The time series is drawn in the upper portion,
#' with a coloured ribbon bar below showing the trend state sequence. This
#' matches the layout used by [plot.resilience()] and makes it easy to
#' compare trend and resilience classifications side by side.
#'
#' @export
#' @param x \[`trend`\]\cr
#'   An object of class `trend` as returned by [compute_trend()].
#' @param type \[`character(1)`: `"series"`\]\cr
#'   Plot type: `"series"` for background shading, `"ribbons"` for a
#'   ribbon bar beneath the time series.
#' @param ... Additional arguments (currently unused).
#' @return A [ggplot2::ggplot()] object.
#'
#' @seealso [compute_trend()] for computing trend classifications;
#'   [plot.resilience()] for the analogous resilience ribbon plot.
#' @family trend
#' @concept time series
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(200, sd = 2))
#' tr <- compute_trend(x, window = 20)
#' plot(tr)
#' plot(tr, type = "ribbons")
#' }
plot.trend <- function(x, type = "series", ...) {
  check_missing(x)
  check_class(x, "trend")
  type <- check_match(type, c("series", "ribbons"))

  if (type == "series") {
    plot_trend_series_(x, ...)
  } else {
    plot_trend_ribbons_(x, ...)
  }
}

#' @noRd
plot_trend_series_ <- function(x, ...) {
  colors <- trend_state_colors_()
  present <- levels(x$state)[levels(x$state) %in% unique(as.character(x$state))]

  ggplot2::ggplot(x, ggplot2::aes(
    x = !!rlang::sym("time"), y = !!rlang::sym("value")
  )) +
    ggplot2::geom_tile(
      ggplot2::aes(
        fill   = !!rlang::sym("state"),
        y      = mean(range(!!rlang::sym("value"), na.rm = TRUE)),
        height = diff(range(!!rlang::sym("value"), na.rm = TRUE)) * 1.1
      ),
      alpha = 0.35,
      width = if (nrow(x) > 1L) {
        stats::median(diff(x$time), na.rm = TRUE)
      } else 1
    ) +
    ggplot2::geom_line(linewidth = 0.6, color = "black") +
    ggplot2::scale_fill_manual(
      name   = "Trend",
      values = colors[present],
      drop   = FALSE
    ) +
    ggplot2::labs(
      title = "Time Series Trend Classification",
      x     = "Time",
      y     = "Value"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title      = ggplot2::element_text(size = 14, face = "bold"),
      axis.title      = ggplot2::element_text(color = "black", face = "bold"),
      axis.text       = ggplot2::element_text(color = "black")
    )
}

#' @noRd
plot_trend_ribbons_ <- function(x, ...) {
  colors <- trend_state_colors_()
  state_levels <- c("ascending", "descending", "flat",
                     "turbulent", "Missing_Data", "Initial")
  state_numeric <- match(as.character(x$state), state_levels)

  time_vals <- x$time
  ts_vals   <- x$value
  ts_range  <- diff(range(ts_vals, na.rm = TRUE))
  ts_min    <- min(ts_vals, na.rm = TRUE)
  ts_max    <- max(ts_vals, na.rm = TRUE)

  # Compress time series into upper portion
  ts_compressed_range <- ts_range * 0.7
  ts_new_max <- ts_max
  ts_new_min <- ts_new_max - ts_compressed_range

  # Ribbon sizing
  ribbon_height  <- ts_range * 0.12
  ribbon_y       <- ts_new_min - ts_range * 0.20
  time_diff <- if (length(time_vals) > 1L) {
    stats::median(diff(time_vals), na.rm = TRUE)
  } else 1

  ts_df <- data.frame(time = time_vals, value = ts_vals)
  tile_df <- data.frame(
    x     = time_vals,
    y     = ribbon_y,
    state = x$state
  )

  t_min <- min(time_vals, na.rm = TRUE)
  t_max <- max(time_vals, na.rm = TRUE)
  y_bottom <- ribbon_y - ribbon_height

  present <- state_levels[state_levels %in% unique(as.character(x$state))]

  ggplot2::ggplot(ts_df, ggplot2::aes(
    x = !!rlang::sym("time"), y = !!rlang::sym("value")
  )) +
    ggplot2::geom_line(color = "black", linewidth = 0.5) +
    ggplot2::geom_tile(
      data    = tile_df,
      ggplot2::aes(
        x    = !!rlang::sym("x"),
        y    = !!rlang::sym("y"),
        fill = !!rlang::sym("state")
      ),
      width  = time_diff,
      height = ribbon_height,
      alpha  = 0.85
    ) +
    ggplot2::annotate(
      "text",
      x        = t_max + (t_max - t_min) / 50,
      y        = ribbon_y,
      label    = "TREND",
      hjust    = 0,
      size     = 3.5,
      fontface = "bold",
      color    = "darkblue"
    ) +
    ggplot2::scale_fill_manual(
      name   = "Trend State",
      values = colors[present],
      drop   = FALSE
    ) +
    ggplot2::coord_cartesian(
      xlim = c(t_min, t_max + (t_max - t_min) / 6),
      ylim = c(y_bottom, ts_new_max + ts_range * 0.1),
      clip = "off"
    ) +
    ggplot2::labs(
      title = "Time Series with Trend Ribbon",
      x     = "Time",
      y     = "Value"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title      = ggplot2::element_text(size = 14, face = "bold"),
      axis.title      = ggplot2::element_text(color = "black", face = "bold"),
      axis.text       = ggplot2::element_text(color = "black")
    )
}

#' Print a Trend Object
#'
#' @describeIn compute_trend Print method for `"trend"` objects. Dispatches
#'   to the tibble print method.
#'
#' @export
#' @param x \[`trend`\]\cr
#'   A `trend` object.
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
print.trend <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Summarize Trend Classification Results
#'
#' @describeIn compute_trend Summary method for `"trend"` objects. Shows
#'   the distribution of trend classifications and the analysis settings.
#'
#' @export
#' @param object \[`trend`\]\cr
#'   A `trend` object.
#' @param ... Additional arguments (currently unused).
#' @return A list with elements `counts` (frequency table of states),
#'   `proportions` (relative frequencies), `window`, `method`, and `align`,
#'   returned invisibly.
summary.trend <- function(object, ...) {
  check_missing(object)
  check_class(object, "trend")
  counts <- table(object$state)
  props  <- prop.table(counts)

  cat("Trend Classification Summary\n")
  cat("  Method :", attr(object, "method"), "\n")
  cat("  Window :", attr(object, "window"), "\n")
  cat("  Align  :", attr(object, "align"), "\n")
  cat("  N      :", nrow(object), "\n\n")
  cat("  State distribution:\n")
  for (i in seq_along(counts)) {
    if (counts[i] > 0L) {
      cat(sprintf("    %-12s %5d  (%5.1f%%)\n",
                  names(counts)[i], counts[i], props[i] * 100))
    }
  }
  invisible(list(
    counts      = counts,
    proportions = props,
    window      = attr(object, "window"),
    method      = attr(object, "method"),
    align       = attr(object, "align")
  ))
}
