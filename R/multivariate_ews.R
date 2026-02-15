# ============================================================================
# multivariate_ews.R --Multivariate early warning signals
# Detects EWS in multivariate time series using Min/Max Autocorrelation
# Factors (MAF), PCA, and covariance-based metrics with rolling or
# expanding window analysis.
# ============================================================================

# --------------------------------------------------------------------------
# Exported: detect_multivariate_warnings()
# --------------------------------------------------------------------------

#' Detect Multivariate Early Warning Signals
#'
#' Performs multivariate early warning signal (EWS) analysis on a system
#' represented by multiple time series. The function computes a suite of
#' metrics --- based on Min/Max Autocorrelation Factors (MAF), Principal
#' Component Analysis (PCA), and covariance structure --- over rolling or
#' expanding windows to detect rising instability before critical transitions.
#'
#' @details
#' Multivariate EWS extend univariate indicators to systems where several
#' variables interact. Individual-variable summaries (mean/max of standard
#' deviations and AR(1) coefficients) are complemented by dimension-reduction
#' indicators that capture coordinated changes across variables.
#'
#' The eleven available metrics are:
#'
#' \describe{
#'   \item{meanSD}{Mean standard deviation across all variables in the
#'     window. Rising values indicate broadening fluctuations system-wide.}
#'   \item{maxSD}{Maximum standard deviation among variables. Flags when any
#'     single variable becomes highly volatile.}
#'   \item{meanAR}{Mean lag-1 autoregressive coefficient across variables.
#'     Increasing values signal system-wide critical slowing down.}
#'   \item{maxAR}{Maximum AR(1) coefficient among variables. Detects
#'     critical slowing down concentrated in a single variable.}
#'   \item{eigenMAF}{Ratio of the minimum eigenvalue of the MAF
#'     decomposition to the total. Declining ratios (increasingly negative
#'     after sign inversion) indicate rising spatial correlation
#'     (Dakos et al., 2012).}
#'   \item{mafAR}{AR(1) coefficient of the first MAF component. The MAF
#'     transformation orders components by temporal smoothness, so the first
#'     component captures the slowest mode. Rising mafAR indicates critical
#'     slowing down in the most persistent system dimension.}
#'   \item{mafSD}{Standard deviation of the first MAF component. Sign-inverted
#'     so that rising values indicate increasing fluctuation in the slowest
#'     mode.}
#'   \item{pcaAR}{AR(1) coefficient of the first principal component. PCA
#'     captures the dominant variance direction; rising pcaAR signals that
#'     this direction is becoming more persistent.}
#'   \item{pcaSD}{Standard deviation of the first principal component.
#'     Rising values indicate that the dominant variance mode is
#'     amplifying.}
#'   \item{eigenCOV}{Dominant eigenvalue of the covariance matrix of
#'     standardized data. Increasing values indicate growing variance
#'     concentration along one axis.}
#'   \item{maxCOV}{Maximum off-diagonal element of the covariance matrix.
#'     Rising values signal strengthening cross-variable coupling.}
#' }
#'
#' **Rolling window** (`method = "rolling"`): a fixed-width window slides
#' across the series. Each metric is standardized (z-scored) across all
#' window positions for comparability. Kendall's tau trend statistic
#' is computed for each metric against time; metrics with
#' \eqn{\tau > 0.7} are flagged as showing significant upward trends.
#'
#' **Expanding window** (`method = "expanding"`): the window starts at
#' `min_window` and grows one observation at a time. A two-stage
#' standardization is applied (matching the EWSmethods implementation):
#' (1) each metric's raw value is z-scored against its own expanding
#' history up to that point; (2) the resulting z-scores are z-scored
#' again using an expanding cumulative mean and SD. This double
#' standardization puts all metrics on a common dimensionless scale
#' regardless of their original units. A warning is flagged when the
#' doubly-standardized strength exceeds `threshold` after the
#' `burn_in` period. Classification into system states (Stable,
#' Vulnerable, Warning, Critical, Failing) is based on
#' how many metrics simultaneously cross the threshold.
#'
#' **Dimension reduction series**: regardless of method, the function
#' also extracts the first MAF (MAF1) and first principal component
#' (PC1) time series, standardizes them, and flags points where the
#' absolute scaled value exceeds the threshold. These provide a
#' low-dimensional summary of the system's trajectory.
#'
#' @export
#' @param data \[`data.frame`\]\cr
#'   A data frame with a time column and at least two numeric time series
#'   columns. If `time_col` is not specified, the first column is assumed
#'   to be time.
#' @param time_col \[`character(1)` or `NULL`: `NULL`\]\cr
#'   Name of the time column. If `NULL`, the first column is used.
#' @param metrics \[`character()`: `"all"`\]\cr
#'   Metrics to compute. Pass `"all"` for all eleven, or a subset of:
#'   `"meanSD"`, `"maxSD"`, `"meanAR"`, `"maxAR"`, `"eigenMAF"`,
#'   `"mafAR"`, `"mafSD"`, `"pcaAR"`, `"pcaSD"`, `"eigenCOV"`,
#'   `"maxCOV"`.
#' @param method \[`character(1)`: `"rolling"`\]\cr
#'   Analysis approach: `"rolling"` or `"expanding"`.
#' @param window \[`numeric(1)`: `50`\]\cr
#'   For `"rolling"`: window size as a percentage of total series length.
#'   For `"expanding"`: not used directly (the window grows from
#'   `min_window` to `n`).
#' @param burn_in \[`integer(1)`: `10L`\]\cr
#'   Number of initial observations after which warnings begin to be
#'   flagged. Only used when `method = "expanding"`.
#' @param threshold \[`numeric(1)`: `2`\]\cr
#'   Number of standard deviations for warning detection. Only used when
#'   `method = "expanding"`.
#' @param tail_direction \[`character(1)`: `"one.tailed"`\]\cr
#'   Whether to test for warnings using `"one.tailed"` (positive
#'   deviations only) or `"two.tailed"` (both directions). Only used when
#'   `method = "expanding"`.
#'
#' @return An object of class `"multi_ews"`, a tibble (inheriting from
#'   `tbl_df`, `tbl`, `data.frame`) with metadata stored as attributes.
#'
#'   **Rolling** columns: `time`, `metric`, `score`, `std`.
#'   Attributes: `orig_data`, `cor` (Kendall tau per metric),
#'   `method` (`"rolling"`), `dimension_reduction`.
#'
#'   **Expanding** columns: `time`, `metric`, `score`, `z_score`,
#'   `detected` (integer 0/1).
#'   Attributes: `orig_data`, `threshold`, `tail_direction`,
#'   `classification` (tibble with `time`, `count`, `state`),
#'   `method` (`"expanding"`), `dimension_reduction`.
#'
#' @references
#' Dakos, V., Carpenter, S. R., Brock, W. A., Ellison, A. M., Guttal, V.,
#' Ives, A. R., ... & Scheffer, M. (2012). Methods for detecting early
#' warnings of critical transitions in time series illustrated using
#' simulated ecological data. *PLoS ONE*, 7(7), e41010.
#' \doi{10.1371/journal.pone.0041010}
#'
#' Scheffer, M., Bascompte, J., Brock, W. A., Brovkin, V., Carpenter,
#' S. R., Dakos, V., ... & Sugihara, G. (2009). Early-warning signals
#' for critical transitions. *Nature*, 461(7260), 53--59.
#' \doi{10.1038/nature08227}
#'
#' @seealso [generate_tipping_data()] for creating test data with known
#'   tipping points; [plot.multi_ews()] for visualization.
#' @family multivariate EWS
#' @concept early warning signals
#' @concept multivariate
#' @examples
#' \donttest{
#' tip <- generate_tipping_data(n_time = 100, n_vars = 3, tipping_point = 60)
#' ews <- detect_multivariate_warnings(tip, method = "expanding", window = 50)
#' plot(ews)
#'
#' ews_roll <- detect_multivariate_warnings(tip, method = "rolling", window = 50)
#' plot(ews_roll)
#' }
detect_multivariate_warnings <- function(data,
                                         time_col = NULL,
                                         metrics = "all",
                                         method = "rolling",
                                         window = 50,
                                         burn_in = 10L,
                                         threshold = 2,
                                         tail_direction = "one.tailed") {
  check_missing(data)

  if (!is.data.frame(data) || ncol(data) < 3L) {
    stop_(
      "{.arg data} must be a data frame with a time column and at least
      two time series columns."
    )
  }

  # --- Identify the time column ---
  if (is.null(time_col)) {
    message_("No {.arg time_col} specified. Using the first column as time.")
    time_col_index <- 1L
  } else {
    time_col_index <- which(names(data) == time_col)
    if (length(time_col_index) == 0L) {
      stop_("Column {.val {time_col}} not found in {.arg data}.")
    }
  }

  original_data <- data
  data_internal <- data[, c(time_col_index, setdiff(seq_len(ncol(data)),
                                                     time_col_index))]

  if (!all(vapply(data_internal[, -1L], is.numeric, logical(1L)))) {
    stop_("All time series columns must be numeric.")
  }

  if (any(is.na(data_internal[, -1L]))) {
    stop_(
      "{.arg data} contains missing values.
      Interpolation before analysis is recommended."
    )
  }

  # --- Validate arguments ---
  method <- check_match(method, c("rolling", "expanding"))
  tail_direction <- check_match(tail_direction, c("one.tailed", "two.tailed"))

  available_metrics <- c(
    "meanSD", "maxSD", "meanAR", "maxAR", "eigenMAF", "mafAR",
    "mafSD", "pcaAR", "pcaSD", "eigenCOV", "maxCOV"
  )

  if (length(metrics) == 1L && metrics == "all") {
    metrics <- available_metrics
  } else {
    invalid <- setdiff(metrics, available_metrics)
    if (length(invalid) > 0L) {
      stop_(
        "Invalid metrics: {.val {invalid}}.
        Available: {.val {available_metrics}}."
      )
    }
  }

  # --- Run core computation ---
  ews_results <- mews_calculate(
    data    = data_internal,
    metrics = metrics,
    window = window,
    method  = method,
    burn_in = burn_in
  )

  # --- Process results into long format ---
  if (method == "rolling") {
    # Kendall tau per metric
    correlations <- vapply(
      ews_results$raw[, metrics, drop = FALSE],
      function(x) {
        if (sum(!is.na(x)) < 4L) return(NA_real_)
        stats::cor.test(x, seq_along(x), method = "kendall")$estimate
      },
      numeric(1L)
    )

    long_data <- tidyr::pivot_longer(
      ews_results$raw,
      cols      = -!!rlang::sym("time"),
      names_to  = "metric",
      values_to = "score"
    )

    # Standardize scores per metric
    long_data$std <- NA_real_
    for (m in unique(long_data$metric)) {
      idx <- long_data$metric == m
      vals <- long_data$score[idx]
      if (sum(!is.na(vals)) > 1L &&
          stats::sd(vals, na.rm = TRUE) > 1e-10) {
        long_data$std[idx] <- as.numeric(scale(vals))
      } else {
        long_data$std[idx] <- 0
      }
    }

  } else {
    # Expanding: metric.score values are already z-scored in Stage 1
    # (inside mews_expanding_metrics). Here we apply Stage 2: a second
    # expanding z-score on those z-scores, matching EWSmethods::wMAF().
    long_data <- tidyr::pivot_longer(
      ews_results$raw,
      cols      = -!!rlang::sym("time"),
      names_to  = "metric",
      values_to = "score"
    )

    long_data <- long_data[order(long_data$metric, long_data$time), ]

    # Stage 2: expanding cumulative mean/sd of the Stage-1 z-scores
    long_data <- long_data |>
      dplyr::group_by(!!rlang::sym("metric")) |>
      dplyr::arrange(!!rlang::sym("time"), .by_group = TRUE) |>
      dplyr::mutate(
        rolling.mean = mews_rolling_mean(!!rlang::sym("score")),
        rolling.sd   = mews_rolling_sd(!!rlang::sym("score"))
      ) |>
      dplyr::ungroup()

    long_data$z_score <- (long_data$score - long_data$rolling.mean) /
      long_data$rolling.sd

    burn_in_time <- data_internal[[1L]][burn_in]

    if (tail_direction == "two.tailed") {
      crossed <- long_data$z_score > threshold | long_data$z_score < -threshold
    } else {
      crossed <- long_data$z_score > threshold
    }
    crossed[is.na(crossed)] <- FALSE

    long_data$detected <- ifelse(
      long_data$time < burn_in_time, 0L,
      ifelse(crossed, 1L, 0L)
    )

    # Remove intermediate columns
    long_data$rolling.mean <- NULL
    long_data$rolling.sd <- NULL
  }

  # --- Build output tibble ---
  out_tbl <- tibble::as_tibble(long_data)

  # --- Classify warnings (expanding only) ---
  classification <- NULL
  if (method == "expanding") {
    n_metrics_count <- length(unique(out_tbl$metric))
    classification <- mews_classify(out_tbl, n_metrics_count)
  }

  # --- Extract MAF1/PC1 dimension reduction series ---
  dr <- mews_dimension_reduction(
    data      = data_internal,
    window   = if (method == "rolling") window else NULL,
    method    = method,
    threshold = if (method == "expanding") threshold else 2
  )

  # --- Assemble as tibble with attributes ---
  if (method == "rolling") {
    structure(
      out_tbl,
      orig_data           = original_data,
      cor                 = correlations,
      method              = "rolling",
      dimension_reduction = dr,
      class               = c("multi_ews", "tbl_df", "tbl", "data.frame")
    )
  } else {
    structure(
      out_tbl,
      orig_data           = original_data,
      threshold           = threshold,
      tail_direction      = tail_direction,
      classification      = classification,
      method              = "expanding",
      dimension_reduction = dr,
      class               = c("multi_ews", "tbl_df", "tbl", "data.frame")
    )
  }
}


# --------------------------------------------------------------------------
# S3: plot.multi_ews()
# --------------------------------------------------------------------------

#' Plot Multivariate EWS Results
#'
#' Creates a multi-panel visualization of multivariate early warning signal
#' analysis results. For expanding windows, shows the standardized metric
#' strengths with threshold lines and a system-state classification ribbon.
#' For rolling windows, shows faceted metric trends with Kendall's tau
#' annotations. Both methods optionally include a dimension-reduction panel
#' showing MAF1 and PC1 trajectories.
#'
#' @details
#' **Expanding window layout** (top to bottom):
#' \enumerate{
#'   \item Dimension reduction panel (MAF1 + PC1) with warning points.
#'   \item Standardized metric strengths with threshold line(s); points
#'     mark individual metric threshold crossings.
#'   \item System-state classification ribbon (Stable through Failing).
#' }
#'
#' **Rolling window layout**:
#' \enumerate{
#'   \item Dimension reduction panel (MAF1 + PC1) with warning points.
#'   \item Faceted panels per metric showing standardized EWS values;
#'     panel titles include Kendall's tau trend statistic.
#' }
#'
#' @export
#' @param x \[`multi_ews`\]\cr
#'   An object of class `"multi_ews"` produced by
#'   [detect_multivariate_warnings()].
#' @param include_dr \[`logical(1)`: `TRUE`\]\cr
#'   Whether to include the dimension reduction (MAF1/PC1) panel.
#' @param ... Additional arguments (currently unused).
#' @return A [ggplot2::ggplot()] or [patchwork::wrap_plots()] object.
#'
#' @seealso [detect_multivariate_warnings()] for computing the EWS.
#' @family multivariate EWS
#' @concept early warning signals
#' @examples
#' \donttest{
#' tip <- generate_tipping_data(n_time = 100, n_vars = 3, tipping_point = 60)
#' ews <- detect_multivariate_warnings(tip, method = "expanding", window = 50)
#' plot(ews)
#' plot(ews, include_dr = FALSE)
#' }
plot.multi_ews <- function(x, include_dr = TRUE, ...) {
  check_missing(x)
  check_class(x, "multi_ews")
  check_flag(include_dr)

  method <- attr(x, "method")

  if (method == "expanding") {
    mews_plot_expanding_(x, include_dr)
  } else {
    mews_plot_rolling_(x, include_dr)
  }
}


#' Print a Multivariate EWS Object
#'
#' @describeIn detect_multivariate_warnings Print method for `"multi_ews"`
#'   objects. Delegates to tibble's default print.
#'
#' @export
#' @param x \[`multi_ews`\]\cr
#'   A `multi_ews` object.
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
print.multi_ews <- function(x, ...) {
  check_missing(x)
  check_class(x, "multi_ews")
  NextMethod(generic = "print", object = x, ...)
}


#' Summarize a Multivariate EWS Object
#'
#' @describeIn detect_multivariate_warnings Summary method for `"multi_ews"`
#'   objects. Provides a condensed overview including per-metric warning
#'   counts and system state distribution.
#'
#' @export
#' @param object \[`multi_ews`\]\cr
#'   A `multi_ews` object.
#' @param ... Additional arguments (currently unused).
#' @return `object`, invisibly.
summary.multi_ews <- function(object, ...) {
  method <- attr(object, "method")
  cat("Multivariate EWS Summary\n")
  cat(strrep("=", 40L), "\n")
  cat("  Method:", method, "\n")
  cat("  Total time points:", length(unique(object$time)), "\n")
  cat("  Metrics computed:", paste(unique(object$metric), collapse = ", "),
      "\n\n")

  # Per-metric warning counts
  if (method == "expanding") {
    cat("Per-metric warnings:\n")
    warned <- object[object$detected == 1L, ]
    if (nrow(warned) > 0L) {
      tbl <- table(warned$metric)
      for (nm in names(tbl)) {
        cat("  ", nm, ":", tbl[[nm]], "time points\n")
      }
    } else {
      cat("  No warnings detected.\n")
    }

    cls <- attr(object, "classification")
    if (!is.null(cls)) {
      cat("\nSystem state distribution:\n")
      tbl <- table(cls$state)
      total <- sum(tbl)
      for (nm in names(tbl)) {
        pct <- round(100 * tbl[[nm]] / total, 1L)
        cat("  ", nm, ":", tbl[[nm]], "(", pct, "%)\n")
      }
    }
  }

  if (method == "rolling") {
    cor_vals <- attr(object, "cor")
    if (!is.null(cor_vals)) {
      cat("Kendall's tau trend statistics:\n")
      tau_clean <- cor_vals
      names(tau_clean) <- sub("\\.tau$", "", names(tau_clean))
      strong <- names(tau_clean)[tau_clean > 0.7 & !is.na(tau_clean)]
      if (length(strong) > 0L) {
        cat("  Strong upward trends (tau > 0.7):", paste(strong, collapse = ", "),
            "\n")
      } else {
        cat("  No metrics show strong upward trends.\n")
      }
    }
  }

  invisible(object)
}


# ==========================================================================
# Internal: core computation engine
# ==========================================================================

#' @noRd
mews_calculate <- function(data, metrics, window, method, burn_in = 10L) {
  time_vec <- data[, 1L]
  ts_data  <- as.matrix(data[, -1L])

  if (method == "rolling") {
    win_size_pts <- round(nrow(ts_data) * window / 100)
    min_window <- max(3L, ncol(ts_data) + 1L)
    if (win_size_pts < min_window) {
      warning_(
        "Window size too small ({win_size_pts} points).
        Using minimum of {min_window} points."
      )
      win_size_pts <- min_window
    }
    if (win_size_pts > nrow(ts_data)) {
      stop_("Window size ({win_size_pts}) exceeds data length ({nrow(ts_data)}).")
    }

    n_obs <- nrow(ts_data)
    rolling_list <- lapply(seq(win_size_pts, n_obs), function(i) {
      w <- ts_data[(i - win_size_pts + 1L):i, , drop = FALSE]
      mews_window_metrics(w, metrics)
    })
    rolling_metrics <- do.call(rbind, rolling_list)
    results_df <- as.data.frame(rolling_metrics)
    results_df$time <- time_vec[win_size_pts:length(time_vec)]
    list(raw = results_df)

  } else {
    # Expanding window: loop starts from burn_in (matching EWSmethods).
    # min_window is the structural minimum (need more obs than variables).
    min_window <- max(3L, ncol(ts_data) + 1L)
    start_from <- max(min_window, burn_in)
    if (start_from >= nrow(ts_data)) {
      stop_("Not enough data points ({nrow(ts_data)}) for expanding analysis.")
    }

    # Accumulators for incremental z-scoring
    env <- new.env(parent = emptyenv())
    env$roll_eigen     <- list()
    env$roll_ar        <- list()
    env$roll_sd        <- list()
    env$roll_pca_ar    <- list()
    env$roll_pca_sd    <- list()
    env$roll_cov       <- list()
    env$roll_mean_ar   <- list()
    env$roll_max_ar    <- list()
    env$roll_mean_sd   <- list()
    env$roll_max_sd    <- list()
    env$roll_cov_eigen <- list()

    results_list <- lapply(start_from:nrow(ts_data), function(i) {
      window_data <- ts_data[1L:i, , drop = FALSE]
      mews_expanding_metrics(window_data, metrics, i, env)
    })

    results_df <- as.data.frame(do.call(rbind, results_list))
    results_df$time <- time_vec[start_from:nrow(ts_data)]
    list(raw = results_df)
  }
}


# ==========================================================================
# Internal: per-window metric computation (rolling)
# ==========================================================================

#' @noRd
mews_window_metrics <- function(window_data, metrics) {
  if (nrow(window_data) < 3L) return(rep(NA_real_, length(metrics)))

  out <- rep(NA_real_, length(metrics))
  names(out) <- metrics

  for (metric in metrics) {
    if (metric %in% c("meanSD", "maxSD")) {
      sds <- apply(window_data, 2L, stats::sd, na.rm = TRUE)
      if (metric == "meanSD") out[metric] <- mean(sds, na.rm = TRUE)
      if (metric == "maxSD")  out[metric] <- max(sds, na.rm = TRUE)
    }

    if (metric %in% c("meanAR", "maxAR")) {
      ars <- apply(window_data, 2L, mews_ar_robust)
      if (metric == "meanAR") out[metric] <- mean(ars, na.rm = TRUE)
      if (metric == "maxAR")  out[metric] <- max(ars, na.rm = TRUE)
    }

    if (metric %in% c("pcaAR", "pcaSD", "eigenCOV", "maxCOV")) {
      window_scaled <- scale(window_data)
      V_scaled <- stats::cov(window_scaled)

      if (all(!is.na(V_scaled))) {
        if (metric %in% c("pcaAR", "pcaSD")) {
          pca_res <- stats::prcomp(window_scaled, scale. = FALSE,
                                   center = FALSE)
          pca_ts <- pca_res$x[, 1L]
          if (metric == "pcaAR") out[metric] <- mews_ar_robust(pca_ts)
          if (metric == "pcaSD") out[metric] <- stats::sd(pca_ts, na.rm = TRUE)
        }
        if (metric == "eigenCOV") {
          out[metric] <- eigen(V_scaled)$values[1L]
        }
        if (metric == "maxCOV") {
          out[metric] <- max(V_scaled[lower.tri(V_scaled, diag = FALSE)])
        }
      }
    }

    if (metric %in% c("eigenMAF", "mafAR", "mafSD")) {
      maf_res <- mews_maf(window_data)
      if (!is.null(maf_res)) {
        if (metric == "eigenMAF") {
          # Negate so that increasing values indicate rising spatial
          # correlation (approaching transition), matching EWSmethods
          out[metric] <- -min(maf_res$eigen_values) /
            sum(maf_res$eigen_values)
        }
        if (metric %in% c("mafAR", "mafSD")) {
          maf_ts_1 <- maf_res$mafs[, 1L]
          if (metric == "mafAR") out[metric] <- mews_ar_robust(maf_ts_1)
          if (metric == "mafSD") out[metric] <- stats::sd(maf_ts_1,
                                                           na.rm = TRUE)
        }
      }
    }
  }

  out
}


# ==========================================================================
# Internal: per-step metric computation (expanding) with incremental z-score
# ==========================================================================

#' @noRd
mews_expanding_metrics <- function(window_data, metrics, i, env) {
  # Stage 1 z-score: each metric's raw value is z-scored against its own

  # expanding history, matching EWSmethods::wMAF() exactly.
  # All metrics use the same formula; only eigenMAF gets sign-inverted.
  result_vec <- rep(NA_real_, length(metrics))
  names(result_vec) <- metrics

  zscore_ <- function(val, acc_list) {
    vals <- unlist(acc_list)
    if (length(vals) > 1L) {
      (val - mean(vals, na.rm = TRUE)) / stats::sd(vals, na.rm = TRUE)
    } else {
      NA_real_
    }
  }

  for (metric in metrics) {
    if (metric %in% c("meanSD", "maxSD")) {
      sds <- apply(window_data, 2L, stats::sd, na.rm = TRUE)
      if (metric == "meanSD") {
        env$roll_mean_sd[[i]] <- mean(sds, na.rm = TRUE)
        result_vec[metric] <- zscore_(env$roll_mean_sd[[i]],
                                       env$roll_mean_sd)
      }
      if (metric == "maxSD") {
        env$roll_max_sd[[i]] <- max(sds, na.rm = TRUE)
        result_vec[metric] <- zscore_(env$roll_max_sd[[i]],
                                       env$roll_max_sd)
      }
    }

    if (metric %in% c("meanAR", "maxAR")) {
      ars <- apply(window_data, 2L, mews_ar_robust)
      if (metric == "meanAR") {
        env$roll_mean_ar[[i]] <- mean(ars, na.rm = TRUE)
        result_vec[metric] <- zscore_(env$roll_mean_ar[[i]],
                                       env$roll_mean_ar)
      }
      if (metric == "maxAR") {
        env$roll_max_ar[[i]] <- max(ars, na.rm = TRUE)
        result_vec[metric] <- zscore_(env$roll_max_ar[[i]],
                                       env$roll_max_ar)
      }
    }

    if (metric %in% c("eigenMAF", "mafAR", "mafSD")) {
      maf_res <- mews_maf(window_data)
      if (!is.null(maf_res)) {
        if (metric == "eigenMAF") {
          env$roll_eigen[[i]] <- min(maf_res$eigen_values /
                                       sum(maf_res$eigen_values))
          result_vec[metric] <- -1 * zscore_(env$roll_eigen[[i]],
                                              env$roll_eigen)
        }
        if (metric == "mafAR") {
          env$roll_ar[[i]] <- mews_ar_robust(maf_res$mafs[, 1L])
          result_vec[metric] <- zscore_(env$roll_ar[[i]], env$roll_ar)
        }
        if (metric == "mafSD") {
          env$roll_sd[[i]] <- stats::sd(maf_res$mafs[, 1L], na.rm = TRUE)
          result_vec[metric] <- zscore_(env$roll_sd[[i]], env$roll_sd)
        }
      }
    }

    if (metric %in% c("pcaAR", "pcaSD")) {
      pca_res <- stats::prcomp(scale(window_data), scale. = FALSE,
                                center = FALSE)
      if (metric == "pcaAR") {
        env$roll_pca_ar[[i]] <- mews_ar_robust(pca_res$x[, 1L])
        result_vec[metric] <- zscore_(env$roll_pca_ar[[i]], env$roll_pca_ar)
      }
      if (metric == "pcaSD") {
        env$roll_pca_sd[[i]] <- stats::sd(pca_res$x[, 1L], na.rm = TRUE)
        result_vec[metric] <- zscore_(env$roll_pca_sd[[i]], env$roll_pca_sd)
      }
    }

    if (metric %in% c("eigenCOV", "maxCOV")) {
      V_scaled <- stats::cov(scale(window_data))
      if (metric == "eigenCOV") {
        env$roll_cov_eigen[[i]] <- max(eigen(V_scaled)$values)
        result_vec[metric] <- zscore_(env$roll_cov_eigen[[i]],
                                       env$roll_cov_eigen)
      }
      if (metric == "maxCOV") {
        env$roll_cov[[i]] <- max(V_scaled[lower.tri(V_scaled, diag = FALSE)])
        result_vec[metric] <- zscore_(env$roll_cov[[i]], env$roll_cov)
      }
    }
  }

  result_vec
}


# ==========================================================================
# Internal: MAF computation
# ==========================================================================

#' @noRd
mews_maf <- function(window_data) {
  if (nrow(window_data) < 3L || ncol(window_data) < 2L) return(NULL)
  result <- try_({
    p <- ncol(window_data)
    n <- nrow(window_data)

    x <- scale(window_data, center = TRUE, scale = TRUE)
    if (any(is.na(x)) || any(!is.finite(x))) return(NULL)

    V <- stats::cov(x)
    svd_V <- svd(V)
    if (any(svd_V$d < .Machine$double.eps)) return(NULL)

    a <- svd_V$u %*% diag(svd_V$d^(-0.5)) %*% t(svd_V$u)
    y <- x %*% a
    dy <- apply(y, 2L, diff)
    cov_dy <- stats::cov(dy) * (n - 2L)
    svd_cov_dy <- svd(cov_dy / (n - 1L))

    u_dy <- svd_cov_dy$u[, p:1L, drop = FALSE]
    aa <- a %*% u_dy
    aa <- apply(aa, 2L, function(col) col / sqrt(sum(col^2)))

    maf <- x %*% aa

    # Sign correction: ensure positive correlation with time
    neg_time_cor <- diag(
      stats::cov(maf, matrix(rep(seq_len(n), p), n, p))
    ) < 0
    mafs <- t(apply(maf, 1L, function(row) {
      (-neg_time_cor + !neg_time_cor) * row
    }))

    d <- svd_cov_dy$d[p:1L]

    list(mafs = mafs, eigen_values = d)
  })
  if (inherits(result, "try-error")) return(NULL)
  result
}


# ==========================================================================
# Internal: robust AR(1) coefficient
# ==========================================================================

#' @noRd
mews_ar_robust <- function(x) {
  if (length(x) < 3L || all(is.na(x)) ||
      stats::sd(x, na.rm = TRUE) < 1e-10) {
    return(NA_real_)
  }
  result <- try_({
    ar_result <- stats::ar.ols(x, aic = FALSE, order.max = 1L,
                                dmean = FALSE, intercept = FALSE)
    if (length(ar_result$ar) > 0L) ar_result$ar[1L] else NA_real_
  })
  if (inherits(result, "try-error")) return(NA_real_)
  result
}


# ==========================================================================
# Internal: rolling cumulative mean / sd
# ==========================================================================

#' @noRd
mews_rolling_mean <- function(x) {
  k <- length(x)
  result <- numeric(k)
  for (i in seq_len(k)) {
    result[i] <- mean(x[1L:i], na.rm = TRUE)
  }
  result
}

#' @noRd
mews_rolling_sd <- function(x) {
  k <- length(x)
  result <- numeric(k)
  for (i in seq_len(k)) {
    result[i] <- stats::sd(x[1L:i], na.rm = TRUE)
  }
  result
}


# ==========================================================================
# Internal: warning classification (expanding only)
# ==========================================================================

#' @noRd
mews_classify <- function(ews_data, n_metrics) {
  ews_valid <- ews_data[!is.na(ews_data$z_score), ]
  if (nrow(ews_valid) == 0L) return(NULL)

  summary_by_time <- stats::aggregate(
    detected ~ time,
    data = ews_valid,
    FUN  = sum
  )
  names(summary_by_time) <- c("time", "count")

  if (n_metrics == 1L) {
    summary_by_time$state <- cut(
      summary_by_time$count,
      breaks = c(-Inf, 0, 1, Inf),
      labels = c("Stable", "Warning", "Critical"),
      right  = TRUE
    )
  } else {
    summary_by_time$state <- cut(
      summary_by_time$count,
      breaks = c(-Inf, 0, 1, 2, 0.5 * n_metrics, Inf),
      labels = c("Stable", "Vulnerable", "Warning",
                  "Critical", "Failing"),
      right  = TRUE
    )
  }

  tibble::as_tibble(summary_by_time)
}


# ==========================================================================
# Internal: dimension reduction (MAF1 + PC1) extraction
# ==========================================================================

#' @noRd
mews_dimension_reduction <- function(data, window, method, threshold) {
  time_vec <- data[, 1L]
  ts_data  <- as.matrix(data[, -1L])

  maf_pc_for_window <- function(window_data) {
    if (nrow(window_data) < 3L || ncol(window_data) < 2L) return(NULL)
    result <- try_({
      p <- ncol(window_data)
      n <- nrow(window_data)

      x <- scale(window_data, center = TRUE, scale = TRUE)
      if (any(is.na(x)) || any(!is.finite(x))) return(NULL)

      V <- stats::cov(x)
      svd_V <- svd(V)
      if (any(svd_V$d < .Machine$double.eps)) return(NULL)

      a <- svd_V$u %*% diag(svd_V$d^(-0.5)) %*% t(svd_V$u)
      y <- x %*% a
      dy <- apply(y, 2L, diff)
      cov_dy <- stats::cov(dy) * (n - 2L)
      svd_cov_dy <- svd(cov_dy / (n - 1L))
      u_dy <- svd_cov_dy$u[, p:1L, drop = FALSE]
      aa <- a %*% u_dy
      aa <- apply(aa, 2L, function(col) col / sqrt(sum(col^2)))
      maf <- x %*% aa
      neg_time_cor <- diag(
        stats::cov(maf, matrix(rep(seq_len(n), p), n, p))
      ) < 0
      mafs <- t(apply(maf, 1L, function(row) {
        (-neg_time_cor + !neg_time_cor) * row
      }))

      pca_res <- stats::prcomp(x, scale. = FALSE, center = FALSE)
      pc1_ts <- pca_res$x[, 1L]

      list(maf1 = mafs[, 1L], pc1 = pc1_ts)
    })
    if (inherits(result, "try-error")) return(NULL)
    result
  }

  if (method == "rolling") {
    win_size_pts <- round(nrow(ts_data) * window / 100)
    min_window <- max(3L, ncol(ts_data) + 1L)
    if (win_size_pts < min_window) win_size_pts <- min_window

    n_obs_dr <- nrow(ts_data)
    maf_pc_list <- lapply(seq(win_size_pts, n_obs_dr), function(i) {
      w <- ts_data[(i - win_size_pts + 1L):i, , drop = FALSE]
      res <- maf_pc_for_window(w)
      if (!is.null(res)) {
        c(maf1 = res$maf1[nrow(w)], pc1 = res$pc1[nrow(w)])
      } else {
        c(maf1 = NA_real_, pc1 = NA_real_)
      }
    })
    maf_pc_results <- do.call(rbind, maf_pc_list)

    results_df <- as.data.frame(maf_pc_results)
    results_df$time <- time_vec[win_size_pts:length(time_vec)]

  } else {
    min_window <- max(3L, ncol(ts_data) + 1L)

    results_list <- lapply(min_window:nrow(ts_data), function(i) {
      window_data <- ts_data[1L:i, , drop = FALSE]
      res <- maf_pc_for_window(window_data)
      if (!is.null(res)) {
        c(time = time_vec[i], maf1 = res$maf1[i], pc1 = res$pc1[i])
      } else {
        c(time = time_vec[i], maf1 = NA_real_, pc1 = NA_real_)
      }
    })

    results_df <- as.data.frame(do.call(rbind, results_list))
  }

  # Standardize and detect warnings
  for (col in c("maf1", "pc1")) {
    if (col %in% names(results_df)) {
      scaled_vals <- as.numeric(scale(results_df[[col]]))
      results_df[[paste0(col, "_scaled")]]  <- scaled_vals
      results_df[[paste0(col, "_warning")]] <- abs(scaled_vals) > threshold
    }
  }

  results_df
}


# ==========================================================================
# Internal: plot helpers
# ==========================================================================

#' @noRd
mews_theme <- function(base_size = 11L) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey90",
                                                  linetype = "dashed"),
      panel.background   = ggplot2::element_rect(fill = "white"),
      panel.border       = ggplot2::element_blank(),
      axis.line          = ggplot2::element_line(colour = "black"),
      strip.background   = ggplot2::element_rect(fill = "grey95",
                                                  colour = "grey95",
                                                  linewidth = 0.5),
      strip.text         = ggplot2::element_text(face = "bold",
                                                  margin = ggplot2::margin(
                                                    5, 0, 5, 0)),
      plot.title         = ggplot2::element_text(face = "bold",
                                                  size = ggplot2::rel(1.2),
                                                  hjust = 0),
      plot.subtitle      = ggplot2::element_text(color = "grey40",
                                                  size = ggplot2::rel(0.9)),
      legend.position    = "left",
      legend.key         = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

#' @noRd
mews_state_colors <- function() {
  c(
    "Stable"     = "#440154FF",
    "Vulnerable" = "#3B528BFF",
    "Warning"    = "#21908CFF",
    "Critical"   = "#5DC863FF",
    "Failing"    = "#FDE725FF"
  )
}


#' @noRd
mews_plot_expanding_ <- function(ews_object, include_dr) {
  classification_data <- attr(ews_object, "classification")
  state_colors <- mews_state_colors()

  if (!is.null(classification_data)) {
    classification_data$state <- factor(
      classification_data$state, levels = names(state_colors)
    )
  }

  # Main metric strength plot
  p <- ggplot2::ggplot(
    ews_object,
    ggplot2::aes(
      x     = !!rlang::sym("time"),
      y     = !!rlang::sym("z_score"),
      color = !!rlang::sym("metric")
    )
  ) +
    ggplot2::geom_line(alpha = 0.8) +
    ggplot2::geom_point(
      data  = function(d) d[d$detected == 1L, ],
      ggplot2::aes(color = !!rlang::sym("metric")),
      size  = 3.5,
      shape = 19,
      alpha = 0.8
    ) +
    ggplot2::geom_hline(
      yintercept = attr(ews_object, "threshold"),
      linetype   = "dashed",
      color      = "grey50"
    )

  if (attr(ews_object, "tail_direction") == "two.tailed") {
    p <- p + ggplot2::geom_hline(
      yintercept = -attr(ews_object, "threshold"),
      linetype   = "dashed",
      color      = "grey50"
    )
  }

  p <- p +
    ggplot2::labs(
      title    = "Multivariate Early Warning Signals (Expanding Window)",
      subtitle = "Points indicate metric strength crossing the threshold.",
      y        = "Strength of EWS",
      x        = NULL,
      color    = "EWS Indicator"
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(shape = 19, size = 3)
      )
    ) +
    mews_theme() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  # Dimension reduction panel
  if (include_dr && !is.null(attr(ews_object, "dimension_reduction"))) {
    p_dr <- mews_plot_dr_(ews_object)

    if (!is.null(classification_data)) {
      p_ribbon <- mews_plot_ribbon_(classification_data, state_colors)
      patchwork::wrap_plots(p_dr, p, p_ribbon, ncol = 1L,
                            heights = c(2, 4, 1))
    } else {
      patchwork::wrap_plots(p_dr, p, ncol = 1L, heights = c(2, 4))
    }
  } else {
    if (!is.null(classification_data)) {
      p_ribbon <- mews_plot_ribbon_(classification_data, state_colors)
      patchwork::wrap_plots(p, p_ribbon, ncol = 1L, heights = c(4, 1))
    } else {
      p + ggplot2::labs(x = "Time Point")
    }
  }
}


#' @noRd
mews_plot_rolling_ <- function(ews_object, include_dr) {
  # Build correlation labels
  cor_vals <- attr(ews_object, "cor")
  cor_names <- sub("\\.tau$", "", names(cor_vals))
  cor_data <- data.frame(
    metric = cor_names,
    tau    = round(cor_vals, 2L),
    stringsAsFactors = FALSE
  )
  cor_data$label <- paste0(cor_data$metric, ": tau=", cor_data$tau)

  plot_data <- merge(as.data.frame(ews_object), cor_data, by = "metric")

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x     = !!rlang::sym("time"),
      y     = !!rlang::sym("std"),
      color = !!rlang::sym("metric")
    )
  ) +
    ggplot2::geom_line(alpha = 0.8, linewidth = 0.5) +
    ggplot2::facet_wrap(
      ggplot2::vars(!!rlang::sym("label")),
      scales = "free_y"
    ) +
    ggplot2::labs(
      title    = "Multivariate Early Warning Signals (Rolling Window)",
      subtitle = paste0("Each panel shows the trend of an EWS metric.",
                        " Kendall's tau indicates trend strength."),
      y        = "Standardized EWS",
      x        = if (include_dr) NULL else "Time Point"
    ) +
    mews_theme() +
    ggplot2::theme(
      legend.position  = "none",
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(size = 8)
    )

  if (include_dr) {
    p <- p + ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }

  if (include_dr && !is.null(attr(ews_object, "dimension_reduction"))) {
    p_dr <- mews_plot_dr_(ews_object) +
      ggplot2::theme(legend.position = "right")
    patchwork::wrap_plots(p_dr, p, ncol = 1L, heights = c(2, 3.5))
  } else {
    p
  }
}


#' @noRd
mews_plot_dr_ <- function(ews_object) {
  dr_data <- attr(ews_object, "dimension_reduction")
  if (is.null(dr_data)) {
    stop_("No dimension reduction data found in {.arg ews_object}.")
  }

  plot_df <- tidyr::pivot_longer(
    dr_data,
    cols      = c(
      !!rlang::sym("maf1_scaled"),
      !!rlang::sym("pc1_scaled")
    ),
    names_to  = "dimension",
    values_to = "scaled_value"
  )

  plot_df$dimension <- ifelse(
    plot_df$dimension == "maf1_scaled", "MAF1", "PC1"
  )
  plot_df$warning <- ifelse(
    plot_df$dimension == "MAF1",
    plot_df$maf1_warning,
    plot_df$pc1_warning
  )

  y_vals <- plot_df$scaled_value[!is.na(plot_df$scaled_value)]
  if (length(y_vals) > 0L) {
    y_range  <- range(y_vals, na.rm = TRUE)
    y_buffer <- diff(y_range) * 0.05
    y_min    <- max(y_range[1L] - y_buffer, -6)
    y_max    <- min(y_range[2L] + y_buffer,  6)
  } else {
    y_min <- -4
    y_max <-  4
  }

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x     = !!rlang::sym("time"),
      y     = !!rlang::sym("scaled_value"),
      color = !!rlang::sym("dimension")
    )
  ) +
    ggplot2::geom_line(alpha = 0.8, linewidth = 1) +
    ggplot2::geom_point(
      data  = function(d) d[d$warning == TRUE, ],
      ggplot2::aes(
        x     = !!rlang::sym("time"),
        y     = !!rlang::sym("scaled_value"),
        color = !!rlang::sym("dimension")
      ),
      size  = 3,
      alpha = 0.8
    ) +
    ggplot2::coord_cartesian(ylim = c(y_min, y_max)) +
    ggplot2::labs(
      title    = "Dimension Reduction Time Series",
      subtitle = "MAF1 and PC1 components with warning points highlighted",
      y        = "Scaled Component Value",
      x        = "Time Point",
      color    = "Component"
    ) +
    ggplot2::scale_color_manual(
      values = c("MAF1" = "#3B82F6", "PC1" = "#F59E0B")
    ) +
    mews_theme() +
    ggplot2::theme(
      legend.position    = "left",
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    )
}


#' @noRd
mews_plot_ribbon_ <- function(classification_data, state_colors) {
  ggplot2::ggplot(
    classification_data,
    ggplot2::aes(x = !!rlang::sym("time"), y = 1)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = !!rlang::sym("state")),
      height = 1
    ) +
    ggplot2::scale_fill_manual(
      values = state_colors,
      name   = "System State",
      drop   = FALSE
    ) +
    ggplot2::labs(x = "Time Point", y = "") +
    mews_theme() +
    ggplot2::theme(
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )
}
