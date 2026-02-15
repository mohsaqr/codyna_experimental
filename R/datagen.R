# ============================================================================
# datagen.R --Synthetic data generation and time series comparison
# Provides generate_ts_data(), generate_tipping_data(), and compare_ts()
# for testing and validating resilience / early-warning analyses.
# ============================================================================

# --------------------------------------------------------------------------
# Exported: generate_ts_data()
# --------------------------------------------------------------------------

#' Generate Synthetic Time Series Data
#'
#' Creates synthetic time series with known ground-truth phases for
#' testing resilience, classification, and early-warning analyses.
#' Supports two modes: `"resilience"` (dynamic phases with shocks,
#' recovery, volatility) and `"clustered"` (distinct value levels).
#'
#' @details
#' Two generation modes are available:
#'
#' **Resilience mode** (`data_type = "resilience"`): produces a time series
#' composed of randomly ordered dynamic phases, with approximately 50\% of
#' the total length allocated to stable segments. The phases are:
#' \describe{
#'   \item{Stable}{An AR(1) process with a low autoregressive coefficient
#'     (0.2) and small noise (sd = 0.5), centered on a given baseline value.
#'     When `n_stable_levels > 1`, multiple stable segments are created at
#'     distinct operating baselines separated by `level_separation` percent
#'     of `mean_stable_value`.}
#'   \item{Shock + Recovery}{A sharp linear drop to 70\% below baseline
#'     (the Shock sub-phase), followed by a linear climb back to baseline
#'     (the Recovery sub-phase). Both include small Gaussian perturbations.}
#'   \item{Volatile}{An AR(1) process with the same autoregressive
#'     coefficient as Stable but much larger noise (sd = 15), producing
#'     high-amplitude fluctuations around the baseline.}
#'   \item{Turbulent}{A random walk (cumulative sum of Gaussian increments,
#'     sd = 8), creating drift without a restoring force.}
#' }
#' The phases are placed in a random order so that the generated series is
#' suitable for testing blind classification algorithms.
#'
#' **Clustered mode** (`data_type = "clustered"`): creates contiguous
#' segments at 2--6 distinct value levels, each with controlled
#' within-cluster Gaussian noise (sd = 1). The function also appends a
#' `quantile_class` column containing a naive quantile-based classification
#' for benchmarking against the true labels.
#'
#' @export
#' @param n_individuals \[`integer(1)`: `1`\]\cr
#'   Number of distinct individuals to generate.
#' @param data_length \[`integer(1)`: `500`\]\cr
#'   Number of data points per individual.
#' @param data_type \[`character(1)`: `"resilience"`\]\cr
#'   One of `"resilience"` or `"clustered"`.
#' @param n_stable_levels \[`integer(1)`: `3`\]\cr
#'   Number of stable levels (resilience mode).
#' @param mean_stable_value \[`numeric(1)`: `100`\]\cr
#'   Central value of the stable phase.
#' @param level_separation \[`numeric(1)`: `20`\]\cr
#'   Percentage separation between stable levels.
#' @param n_levels \[`integer(1)`: `3`\]\cr
#'   Number of clusters (clustered mode, 2--6).
#' @param results_df \[`data.frame` or `NULL`\]\cr
#'   Optional external results for comparison plotting.
#' @param results_col \[`character(1)` or `NULL`\]\cr
#'   Column in `results_df` with classification labels.
#' @param seed \[`integer(1)`: `42`\]\cr
#'   Random seed for reproducibility.
#' @param generate_plot \[`logical(1)`: `TRUE`\]\cr
#'   If `TRUE`, return a list with `$data` and `$plot`.
#'
#' @return Depends on `generate_plot`:
#' \itemize{
#'   \item If `generate_plot = TRUE`, a `list` with:
#'     \describe{
#'       \item{data}{A `data.frame` with columns `time`, `value`,
#'         `true_phase`, and optionally `id` (when `n_individuals > 1`).
#'         Clustered mode also includes `quantile_class`.}
#'       \item{plot}{A `ggplot` or `patchwork` object showing the time
#'         series with classification bands.}
#'     }
#'   \item If `generate_plot = FALSE`, the `data.frame` directly.
#' }
#'
#' @seealso [compare_ts()] for visualising classification results,
#'   [resilience()] for computing resilience metrics,
#'   [hurst()] for long-range dependence analysis.
#' @family data generation
#' @concept synthetic data
#' @concept time series
#'
#' @examples
#' # Resilience mode --returns data.frame only
#' set.seed(1)
#' df <- generate_ts_data(
#'   n_individuals = 1, data_length = 200,
#'   data_type = "resilience", generate_plot = FALSE
#' )
#' head(df)
#' table(df$true_phase)
#'
#' # Clustered mode --returns data.frame only
#' df_cl <- generate_ts_data(
#'   data_length = 300, data_type = "clustered",
#'   n_levels = 4, generate_plot = FALSE
#' )
#' head(df_cl)
#'
#' \donttest{
#' # Resilience mode with plot
#' out <- generate_ts_data(
#'   n_individuals = 2, data_length = 500,
#'   data_type = "resilience", generate_plot = TRUE
#' )
#' out$plot
#' }
generate_ts_data <- function(n_individuals = 1L,
                             data_length = 500L,
                             data_type = "resilience",
                             n_stable_levels = 3L,
                             mean_stable_value = 100,
                             level_separation = 20,
                             n_levels = 3L,
                             results_df = NULL,
                             results_col = NULL,
                             seed = 42L,
                             generate_plot = TRUE) {
  set.seed(seed)
  data_type <- check_match(data_type, c("resilience", "clustered"))
  check_range(n_individuals, type = "integer", min = 1L, max = 100L)
  check_range(data_length, type = "integer", min = 50L, max = 100000L)
  check_flag(generate_plot)
  if (data_type == "clustered" && !n_levels %in% 2L:6L) {
    stop_("{.arg n_levels} must be between 2 and 6 for clustered data.")
  }
  if (!is.null(results_df) && is.null(results_col)) {
    stop_("When {.arg results_df} is provided, {.arg results_col} is required.")
  }

  all_data <- lapply(seq_len(n_individuals), function(i) {
    set.seed(seed + i)
    df <- if (data_type == "clustered") {
      datagen_clustered_(data_length, n_levels, 100, 2, 1,
                         NULL, 0.1, TRUE)
    } else {
      datagen_resilience_(data_length, n_stable_levels,
                          mean_stable_value, level_separation)
    }
    df$id <- i
    df
  })

  final <- do.call(rbind, all_data)
  final$id <- as.factor(final$id)

  if (!generate_plot) return(final)

  message_("Generating plot for {n_individuals} individual(s).")
  fb <- if (n_individuals > 1L) "id" else NULL

  if (!is.null(results_df)) {
    if (!results_col %in% names(results_df)) {
      stop_("Column {.var {results_col}} not found in {.arg results_df}.")
    }
    p <- compare_ts(
      data = final, data2 = results_df,
      ts_col = "value", time_col = "time",
      state1 = "true_phase", state2 = results_col,
      title1 = "Ground Truth Phases",
      title2 = "External Model Classification",
      facet_by = fb
    )
  } else if (data_type == "clustered") {
    p <- compare_ts(
      data = final, ts_col = "value", time_col = "time",
      state1 = "true_phase", state2 = "quantile_class",
      title1 = "True Generated Clusters",
      title2 = "Benchmark: Quantile Classification",
      facet_by = fb
    )
  } else {
    p <- compare_ts(
      data = final, ts_col = "value", time_col = "time",
      state1 = "true_phase",
      title1 = "Generated Resilience Phases",
      facet_by = fb
    )
  }
  list(data = final, plot = p)
}

# --------------------------------------------------------------------------
# Exported: generate_tipping_data()
# --------------------------------------------------------------------------

#' Generate Multivariate Tipping-Point Data
#'
#' Simulates a multivariate system approaching a critical transition.
#' Two simulation models are available: a mean-reverting
#' Ornstein--Uhlenbeck process (`model = "ou"`) with decaying restoring
#' force, and a classic autoregressive model (`model = "ar"`) with
#' increasing AR coefficient and noise.
#'
#' @details
#' **Ornstein--Uhlenbeck model** (`model = "ou"`, default): each variable
#' follows a discrete mean-reverting process centred on a baseline (20):
#'
#' \deqn{x_t = x_{t-1} + \lambda \, (\mu - x_{t-1}) + \sigma \, \varepsilon_t}
#'
#' Before the tipping point, \eqn{\lambda} equals `stability_strength`
#' (default 0.5) and \eqn{\sigma} = 2, giving a stationary process with
#' moderate autocorrelation. After the tipping point, the restoring force
#' \eqn{\lambda} decays toward near-zero while noise ramps up to 4x.
#' As \eqn{\lambda \to 0}, the AR(1) coefficient approaches 1 (critical
#' slowing down) and variance grows as \eqn{\sigma^2 / (2\lambda)}.
#'
#' **Autoregressive model** (`model = "ar"`): each variable follows a
#' standard AR(1) process:
#'
#' \deqn{x_t = \phi \, x_{t-1} + \varepsilon_t}
#'
#' Before the tipping point, \eqn{\phi} = `stability_strength` (default
#' 0.8 for this model) with innovations \eqn{N(20, 5)}. After the
#' tipping point, \eqn{\phi} increases by `forcing_strength` per step
#' and the noise SD grows proportionally. This produces a system that
#' amplifies perturbations increasingly past the tipping point.
#'
#' Both models cap the forcing at `saturation_point` to prevent
#' divergence.
#'
#' @export
#' @param n_time \[`integer(1)`: `200`\]\cr
#'   Total number of time points.
#' @param n_vars \[`integer(1)`: `5`\]\cr
#'   Number of variables in the system.
#' @param tipping_point \[`integer(1)`: `100`\]\cr
#'   Time point at which forcing begins.
#' @param stability_strength \[`numeric(1)`: `0.5`\]\cr
#'   Restoring-force coefficient (OU model) or AR coefficient (AR model)
#'   during the stable period.
#' @param forcing_strength \[`numeric(1)`: `0.02`\]\cr
#'   Rate at which the system destabilises per timestep after the tipping
#'   point. For the OU model this controls restoring-force decay; for the
#'   AR model it controls AR coefficient growth.
#' @param saturation_point \[`integer(1)`: `180`\]\cr
#'   Time after which the forcing stops increasing.
#' @param model \[`character(1)`: `"ou"`\]\cr
#'   Simulation model. `"ou"` for the Ornstein--Uhlenbeck model (clean
#'   critical slowing down signal), `"ar"` for the classic autoregressive
#'   model (stronger amplitude changes).
#'
#' @return A `data.frame` with columns:
#' \itemize{
#'   \item `Time`: integer time index from 1 to `n_time`.
#'   \item `VAR1`, `VAR2`, \ldots, `VARn`: simulated variable values,
#'     where \eqn{n} = `n_vars`.
#' }
#'
#' @seealso [generate_ts_data()] for generating univariate synthetic
#'   series with labelled resilience phases.
#' @family data generation
#' @concept synthetic data
#' @concept time series
#'
#' @examples
#' # OU model (default) -- clean critical slowing down
#' tp_ou <- generate_tipping_data(n_time = 200, n_vars = 3, tipping_point = 100)
#' head(tp_ou)
#'
#' # AR model -- classic autoregressive with growing instability
#' tp_ar <- generate_tipping_data(
#'   n_time = 100, n_vars = 3, tipping_point = 60,
#'   stability_strength = 0.8, forcing_strength = 0.01,
#'   saturation_point = 80, model = "ar"
#' )
#' head(tp_ar)
generate_tipping_data <- function(n_time = 200L,
                                  n_vars = 5L,
                                  tipping_point = 100L,
                                  stability_strength = 0.5,
                                  forcing_strength = 0.02,
                                  saturation_point = 180L,
                                  model = "ou") {
  model <- check_match(model, c("ou", "ar"))
  check_range(n_time, type = "integer", min = 10L, max = 100000L)
  check_range(n_vars, type = "integer", min = 1L, max = 100L)
  check_range(tipping_point, type = "integer", min = 2L, max = n_time - 1L)
  set.seed(123L)

  mat <- matrix(0, nrow = n_time, ncol = n_vars)

  if (model == "ou") {
    # --- Ornstein-Uhlenbeck model ---
    baseline <- 20
    base_noise <- 2
    mat[1L, ] <- baseline + stats::rnorm(n_vars, 0, base_noise)

    # Stable period: strong mean reversion
    for (i in 2L:tipping_point) {
      mat[i, ] <- mat[i - 1L, ] +
        stability_strength * (baseline - mat[i - 1L, ]) +
        stats::rnorm(n_vars, 0, base_noise)
    }

    # Forcing period: restoring force decays, noise ramps up
    for (i in (tipping_point + 1L):n_time) {
      dt <- min(i - tipping_point, saturation_point - tipping_point)
      progress <- dt / (saturation_point - tipping_point)
      lambda <- stability_strength * (1 - 0.95 * progress)
      noise <- base_noise * (1 + 3 * progress)
      mat[i, ] <- mat[i - 1L, ] +
        lambda * (baseline - mat[i - 1L, ]) +
        stats::rnorm(n_vars, 0, noise)
    }

  } else {
    # --- Classic AR model (original) ---
    mat[1L, ] <- stats::rnorm(n_vars, mean = 20, sd = 5)

    # Stable period
    for (i in 2L:tipping_point) {
      mat[i, ] <- mat[i - 1L, ] * stability_strength +
        stats::rnorm(n_vars, mean = 20, sd = 5)
    }

    # Forcing period
    for (i in (tipping_point + 1L):n_time) {
      dt <- if (i > saturation_point) {
        saturation_point - tipping_point
      } else {
        i - tipping_point
      }
      restore <- stability_strength + dt * forcing_strength
      noise   <- 5 + dt * forcing_strength * 20
      mat[i, ] <- mat[i - 1L, ] * restore +
        stats::rnorm(n_vars, mean = 20, sd = noise)
    }
  }

  out <- as.data.frame(mat)
  names(out) <- paste0("VAR", seq_len(n_vars))
  out$Time <- seq_len(n_time)
  out[, c("Time", setdiff(names(out), "Time"))]
}

# --------------------------------------------------------------------------
# Exported: compare_ts()
# --------------------------------------------------------------------------

#' Compare Time Series Classifications
#'
#' Visualises one or two classification results as coloured background
#' bands behind a time series line plot. Useful for comparing ground-truth
#' labels to model predictions.
#'
#' @details
#' Each classification is rendered as a panel containing the time series
#' drawn as a black line with classification labels shown as coloured
#' background rectangles spanning the full y-axis range. Colours are
#' drawn from [RColorBrewer][RColorBrewer::brewer.pal] palettes
#' (`palette1` and `palette2`).
#'
#' When only `state1` is provided the function returns a single `ggplot`
#' panel. When `state2` is also supplied (either from a column already in
#' `data` or from a separate `data2` data frame), two panels are
#' produced and stacked vertically using
#' [patchwork::wrap_plots()], making it straightforward to compare, for
#' example, ground-truth phases against a model's predicted labels.
#'
#' If `facet_by` is set (e.g., to an individual-ID column), each panel
#' is additionally faceted with free y-scales so that multiple time
#' series can be displayed side by side.
#'
#' @export
#' @param data \[`data.frame`\]\cr
#'   Primary data containing the time series and `state1`.
#' @param data2 \[`data.frame` or `NULL`\]\cr
#'   Optional second data frame containing `state2`.
#' @param ts_col \[`character(1)`\]\cr
#'   Column name for the time series values.
#' @param time_col \[`character(1)` or `NULL`\]\cr
#'   Column name for the time axis. If `NULL`, row index is used.
#' @param state1 \[`character(1)`\]\cr
#'   Column name for the first classification.
#' @param state2 \[`character(1)` or `NULL`\]\cr
#'   Column name for the second classification.
#' @param facet_by \[`character(1)` or `NULL`\]\cr
#'   Column name to facet by (e.g., individual ID).
#' @param title1 \[`character(1)`: `"Classification 1"`\]\cr
#'   Title for the first panel.
#' @param title2 \[`character(1)`: `"Classification 2"`\]\cr
#'   Title for the second panel.
#' @param palette1 \[`character(1)`: `"Set2"`\]\cr
#'   RColorBrewer palette for the first panel.
#' @param palette2 \[`character(1)`: `"Set3"`\]\cr
#'   RColorBrewer palette for the second panel.
#'
#' @return
#' \itemize{
#'   \item When `state2` is `NULL`: a single `ggplot` object.
#'   \item When `state2` is provided: a `patchwork` object containing two
#'     vertically stacked `ggplot` panels.
#' }
#'
#' @seealso [generate_ts_data()] for creating synthetic time series with
#'   known ground-truth phases suitable as input to this function.
#' @family data generation
#' @concept synthetic data
#' @concept time series
#'
#' @examples
#' \donttest{
#' # Generate data and compare ground truth to a naive quantile classifier
#' df <- generate_ts_data(
#'   data_length = 300, data_type = "clustered",
#'   n_levels = 3, generate_plot = FALSE
#' )
#' # Single-panel plot
#' p1 <- compare_ts(
#'   data = df, ts_col = "value", time_col = "time",
#'   state1 = "true_phase", title1 = "True Clusters"
#' )
#' print(p1)
#'
#' # Two-panel comparison
#' p2 <- compare_ts(
#'   data = df, ts_col = "value", time_col = "time",
#'   state1 = "true_phase", state2 = "quantile_class",
#'   title1 = "True Clusters",
#'   title2 = "Quantile Classification"
#' )
#' print(p2)
#' }
compare_ts <- function(data,
                       data2 = NULL,
                       ts_col,
                       time_col = NULL,
                       state1,
                       state2 = NULL,
                       facet_by = NULL,
                       title1 = "Classification 1",
                       title2 = "Classification 2",
                       palette1 = "Set2",
                       palette2 = "Set3") {
  check_missing(data)
  check_missing(ts_col)
  check_missing(state1)
  if (!is.data.frame(data)) stop_("{.arg data} must be a data.frame.")

  plot_df <- data

  if (is.null(time_col)) {
    message_("Argument {.arg time_col} not provided. Inferring time index from row order.")
    time_col <- "inferred_time_index"
    plot_df[[time_col]] <- seq_len(nrow(plot_df))
  }

  if (!is.null(state2)) {
    if (!is.null(data2)) {
      if (!is.data.frame(data2)) stop_("{.arg data2} must be a data.frame.")
      if (!state2 %in% names(data2)) {
        stop_("Column {.var {state2}} not found in {.arg data2}.")
      }
      if (nrow(plot_df) != nrow(data2)) {
        stop_("{.arg data} and {.arg data2} must have the same number of rows.")
      }
      plot_df[[state2]] <- data2[[state2]]
    } else {
      if (!state2 %in% names(plot_df)) {
        stop_("Column {.var {state2}} not found in {.arg data}.")
      }
    }
  }

  cols_needed <- c(ts_col, time_col, state1)
  if (!is.null(state2)) cols_needed <- c(cols_needed, state2)
  if (!is.null(facet_by)) cols_needed <- c(cols_needed, facet_by)
  missing <- cols_needed[!cols_needed %in% names(plot_df)]
  if (length(missing) > 0L) {
    stop_("Missing columns: {.var {paste(missing, collapse = ', ')}}.")
  }

  make_panel <- function(df, st_col, ttl, pal) {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = !!rlang::sym(time_col),
        y = !!rlang::sym(ts_col)
      )
    ) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = !!rlang::sym(time_col) - 0.5,
          xmax = !!rlang::sym(time_col) + 0.5,
          ymin = -Inf, ymax = Inf,
          fill = !!rlang::sym(st_col)
        ),
        alpha = 0.6
      ) +
      ggplot2::geom_line(color = "black", linewidth = 0.7) +
      ggplot2::scale_fill_brewer(palette = pal, name = st_col) +
      ggplot2::labs(title = ttl, x = "Time Index", y = "Value") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(legend.position = "right")

    if (!is.null(facet_by)) {
      p <- p + ggplot2::facet_wrap(
        ggplot2::vars(!!rlang::sym(facet_by)),
        scales = "free_y"
      )
    }
    p
  }

  p1 <- make_panel(plot_df, state1, title1, palette1)
  if (is.null(state2)) return(p1)

  p2 <- make_panel(plot_df, state2, title2, palette2)
  patchwork::wrap_plots(p1, p2, ncol = 1L)
}

# --------------------------------------------------------------------------
# Internal helpers (@noRd)
# --------------------------------------------------------------------------

#' Generate random proportions for cluster segments
#' @noRd
datagen_cluster_proportions_ <- function(n_levels, proportions,
                                         min_cluster_size,
                                         randomize) {
  if (!is.null(proportions)) {
    if (length(proportions) != n_levels) {
      stop_("Length of {.arg proportions} must equal {.arg n_levels}.")
    }
    if (abs(sum(proportions) - 1) > 1e-10) {
      stop_("{.arg proportions} must sum to 1.")
    }
    return(proportions)
  }
  if (!randomize) return(rep(1 / n_levels, n_levels))
  for (attempt in seq_len(100L)) {
    w <- stats::runif(n_levels, 0.5, 2)
    w <- w / sum(w)
    if (all(w >= min_cluster_size)) return(w)
  }
  warning_("Could not satisfy constraints after 100 attempts. Using equal proportions.")
  rep(1 / n_levels, n_levels)
}

#' Generate a single resilience time series with dynamic phases
#' @noRd
datagen_resilience_ <- function(data_length, n_stable_levels,
                                mean_stable_value, level_separation) {
  if (n_stable_levels <= 0L) {
    stop_("{.arg n_stable_levels} must be 1 or greater.")
  }
  # --- Stable level means and names ---
  if (n_stable_levels == 1L) {
    stable_levels <- mean_stable_value
    stable_names  <- "Stable"
  } else {
    offset <- mean_stable_value * (level_separation / 100)
    if (n_stable_levels == 2L) {
      stable_levels <- c(mean_stable_value - offset,
                         mean_stable_value + offset)
      stable_names  <- c("Stable Low", "Stable High")
    } else if (n_stable_levels == 3L) {
      stable_levels <- c(mean_stable_value - offset,
                         mean_stable_value,
                         mean_stable_value + offset)
      stable_names  <- c("Stable Low", "Stable", "Stable High")
    } else {
      steps <- seq(-1.5, 1.5, length.out = n_stable_levels)
      stable_levels <- mean_stable_value + steps * offset
      stable_names  <- paste("Stable", seq_len(n_stable_levels))
    }
  }
  stable_levels <- sort(stable_levels)

  dynamic_names <- c("ShockBlock", "Volatile", "Turbulent")
  all_phase_names <- c(stable_names, dynamic_names)

  # --- Phase proportions (~50% stable) ---
  prop_shock    <- 0.10
  prop_volatile <- 0.20
  prop_turb     <- 0.20
  prop_stable   <- 1 - (prop_shock + prop_volatile + prop_turb)
  stable_props  <- rep(prop_stable / n_stable_levels, n_stable_levels)
  phase_props   <- c(stable_props, prop_shock, prop_volatile, prop_turb)
  phase_props   <- phase_props / sum(phase_props)

  # --- Segment lengths (randomised order) ---
  phase_lengths <- round(phase_props * data_length)
  phase_lengths[length(phase_lengths)] <-
    data_length - sum(phase_lengths[-length(phase_lengths)])

  set.seed(as.numeric(Sys.time()))
  idx <- sample(seq_along(all_phase_names))
  ordered_types   <- all_phase_names[idx]
  ordered_lengths <- phase_lengths[idx]

  starts <- cumsum(c(1L, ordered_lengths[-length(ordered_lengths)]))
  ends   <- cumsum(ordered_lengths)

  # --- Generate values per phase ---
  ts_vals <- numeric(data_length)
  phase_vec <- character(data_length)
  baseline <- mean_stable_value

  for (j in seq_along(ordered_types)) {
    s <- starts[j]; e <- ends[j]; len <- e - s + 1L
    ptype <- ordered_types[j]

    if (ptype == "ShockBlock") {
      len_drop  <- round(len * 0.4)
      len_climb <- len - len_drop
      mid <- s + len_drop - 1L
      low_pt <- baseline - baseline * 0.7
      ts_vals[s:mid] <- seq(baseline - 10, low_pt, length.out = len_drop) +
        stats::rnorm(len_drop, 0, 3)
      phase_vec[s:mid] <- "Shock"
      ts_vals[(mid + 1L):e] <- seq(low_pt, baseline, length.out = len_climb) +
        stats::rnorm(len_climb, 0, 3)
      phase_vec[(mid + 1L):e] <- "Recovery"

    } else if (grepl("Stable", ptype)) {
      si <- which(stable_names == ptype)
      mv <- stable_levels[si]
      ts_vals[s:e] <- mv + stats::arima.sim(
        model = list(ar = 0.2), n = len, sd = 0.5
      )
      phase_vec[s:e] <- ptype

    } else if (ptype == "Volatile") {
      ts_vals[s:e] <- baseline + stats::arima.sim(
        model = list(ar = 0.2), n = len, sd = 15
      )
      phase_vec[s:e] <- "Volatile"

    } else if (ptype == "Turbulent") {
      ts_vals[s:e] <- baseline + cumsum(stats::rnorm(len, 0, 8))
      phase_vec[s:e] <- "Turbulent"
    }
  }

  all_levels <- unique(c(stable_names, "Shock", "Recovery",
                         "Volatile", "Turbulent"))
  message_("Generated resilience data with {length(all_levels)} distinct phases.")
  data.frame(
    time       = seq_len(data_length),
    value      = ts_vals,
    true_phase = factor(phase_vec, levels = all_levels)
  )
}

#' Generate a single clustered time series
#' @noRd
datagen_clustered_ <- function(data_length, n_levels, base_value,
                               cluster_sep, within_sd, proportions,
                               min_cluster_size, randomize) {
  configs <- list(
    "2" = list(nm = c("Low", "High"),
               m  = c(-1, 1) * cluster_sep),
    "3" = list(nm = c("Low", "Medium", "High"),
               m  = c(-1.5, 0, 1.5) * cluster_sep),
    "4" = list(nm = c("Very Low", "Low", "High", "Very High"),
               m  = c(-2, -0.7, 0.7, 2) * cluster_sep),
    "5" = list(nm = c("Very Low", "Low", "Medium", "High", "Very High"),
               m  = c(-2, -1, 0, 1, 2) * cluster_sep),
    "6" = list(nm = paste("Level", 1:6),
               m  = c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5) * cluster_sep)
  )
  cfg <- configs[[as.character(n_levels)]]
  props <- datagen_cluster_proportions_(
    n_levels, proportions, min_cluster_size, randomize
  )
  seg_len <- round(props * data_length)
  seg_len[length(seg_len)] <-
    data_length - sum(seg_len[-length(seg_len)])

  ts_vals <- numeric(data_length)
  labels  <- character(data_length)
  pos <- 1L
  for (i in seq_len(n_levels)) {
    ep <- pos + seg_len[i] - 1L
    ts_vals[pos:ep] <- base_value + cfg$m[i] +
      stats::rnorm(seg_len[i], 0, within_sd)
    labels[pos:ep] <- cfg$nm[i]
    pos <- ep + 1L
  }

  out <- data.frame(
    time       = seq_len(data_length),
    value      = ts_vals,
    true_phase = factor(labels, levels = cfg$nm)
  )
  qb <- stats::quantile(out$value,
                         probs = seq(0, 1, length.out = n_levels + 1L))
  out$quantile_class <- cut(out$value, breaks = qb,
                            labels = cfg$nm, include.lowest = TRUE)
  message_("Generated clustered data with {n_levels} clusters.")
  out
}
