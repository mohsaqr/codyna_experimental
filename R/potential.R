# ============================================================================
# potential.R -- Stability landscape and potential analysis
# Estimates the quasi-potential U(x) = -log(P(x)) from kernel density,
# identifies wells (stable attractors) and barriers (unstable equilibria),
# and tracks landscape topology changes over rolling windows.
# ============================================================================

# --------------------------------------------------------------------------
# Internal helpers (@noRd)
# --------------------------------------------------------------------------

#' Find local minima in a numeric vector
#'
#' Identifies indices where the first difference changes from negative to
#' positive (valley). Boundary points are excluded.
#'
#' @param x Numeric vector.
#' @return Integer vector of indices.
#' @noRd
potential_find_minima_ <- function(x) {
  n <- length(x)
  if (n < 3L) return(integer(0))
  d <- diff(x)
  # Sign changes: negative -> positive
  sign_d <- sign(d)
  # Remove exact zeros by forward-filling
  for (i in seq_along(sign_d)) {
    if (sign_d[i] == 0L && i > 1L) sign_d[i] <- sign_d[i - 1L]
  }
  idx <- which(sign_d[-length(sign_d)] < 0 & sign_d[-1L] > 0) + 1L
  idx
}

#' Find local maxima in a numeric vector
#'
#' Identifies indices where the first difference changes from positive to
#' negative (peak). Boundary points are excluded.
#'
#' @param x Numeric vector.
#' @return Integer vector of indices.
#' @noRd
potential_find_maxima_ <- function(x) {
  n <- length(x)
  if (n < 3L) return(integer(0))
  d <- diff(x)
  sign_d <- sign(d)

  for (i in seq_along(sign_d)) {
    if (sign_d[i] == 0L && i > 1L) sign_d[i] <- sign_d[i - 1L]
  }
  idx <- which(sign_d[-length(sign_d)] > 0 & sign_d[-1L] < 0) + 1L
  idx
}

#' Compute potential landscape from a numeric vector
#'
#' Estimates kernel density, derives potential U(x) = -log(density), and
#' locates wells and barriers.
#'
#' @param x Numeric vector (cleaned, no NAs).
#' @param n_bins Integer grid resolution.
#' @param bandwidth Numeric bandwidth or NULL for Silverman default.
#' @return List with landscape, wells, barriers, n_wells.
#' @noRd
potential_landscape_ <- function(x, n_bins, bandwidth) {
  n_grid <- n_bins * 4L
  dens_args <- list(x = x, n = n_grid)
  if (!is.null(bandwidth)) dens_args$bw <- bandwidth
  dens <- do.call(stats::density, dens_args)

  # Potential: U(x) = -log(P(x)), with epsilon to avoid log(0)
  eps <- .Machine$double.eps
  density_vals <- pmax(dens$y, eps)
  potential_vals <- -log(density_vals)

  landscape <- tibble::tibble(
    x = dens$x,
    density = dens$y,
    potential = potential_vals
  )

  # Find wells (local minima of potential) and barriers (local maxima)
  well_idx <- potential_find_minima_(potential_vals)
  barrier_idx <- potential_find_maxima_(potential_vals)

  # Build wells tibble
  if (length(well_idx) == 0L) {
    # If no interior minimum found, use the global minimum
    global_min <- which.min(potential_vals)
    well_idx <- global_min
  }

  # Compute well properties
  well_locations <- dens$x[well_idx]
  well_depths <- numeric(length(well_idx))
  well_widths <- numeric(length(well_idx))

  # For each well, find depth and width relative to nearest barriers

  # Build ordered set of barrier + boundary positions for depth/width
  all_barriers <- sort(unique(c(1L, barrier_idx, length(potential_vals))))

  for (wi in seq_along(well_idx)) {
    w <- well_idx[wi]
    # Find left barrier (closest barrier index < w)
    left_candidates <- all_barriers[all_barriers < w]
    left_b <- ifelse_(
      length(left_candidates) > 0L,
      left_candidates[length(left_candidates)],
      1L
    )
    # Find right barrier (closest barrier index > w)
    right_candidates <- all_barriers[all_barriers > w]
    right_b <- ifelse_(
      length(right_candidates) > 0L,
      right_candidates[1L],
      length(potential_vals)
    )

    # Depth: min of (left_barrier_height - well_bottom, right_barrier_height - well_bottom)
    well_bottom <- potential_vals[w]
    left_height <- potential_vals[left_b]
    right_height <- potential_vals[right_b]
    well_depths[wi] <- min(left_height, right_height) - well_bottom

    # Width: distance in x between left and right barriers
    well_widths[wi] <- dens$x[right_b] - dens$x[left_b]
  }

  wells <- tibble::tibble(
    location = well_locations,
    depth = well_depths,
    width = well_widths
  )

  # Build barriers tibble
  if (length(barrier_idx) > 0L) {
    barriers <- tibble::tibble(
      location = dens$x[barrier_idx],
      height = potential_vals[barrier_idx]
    )
  } else {
    barriers <- tibble::tibble(
      location = numeric(0),
      height = numeric(0)
    )
  }

  list(
    landscape = landscape,
    wells = wells,
    barriers = barriers,
    n_wells = length(well_idx)
  )
}

#' Detrend a numeric vector
#'
#' @param x Numeric vector.
#' @param method One of "none", "linear", "diff".
#' @return Numeric vector (detrended).
#' @noRd
potential_detrend_ <- function(x, method) {
  switch(
    method,
    none = x,
    linear = {
      t_idx <- seq_along(x)
      fit <- stats::lm.fit(x = cbind(1, t_idx), y = x)
      fit$residuals
    },
    diff = {
      diff(x)
    }
  )
}

#' Colour palette for potential landscape plots
#' @noRd
potential_colors_ <- function() {
  c(
    well = "#2196F3",
    barrier = "#F44336",
    potential = "#333333",
    density = "#4CAF50"
  )
}

# --------------------------------------------------------------------------
# Exported: potential_analysis()
# --------------------------------------------------------------------------

#' Stability Landscape and Potential Analysis
#'
#' Estimates the quasi-potential landscape \eqn{U(x) = -\log P(x)} from
#' kernel density estimation, identifying stable attractors (potential
#' wells) and unstable equilibria (barriers). Supports global analysis of
#' the full series or rolling-window analysis to track how the number and
#' location of attractors change over time.
#'
#' @details
#' The potential function is derived from the Boltzmann relation between
#' probability density and energy: \eqn{P(x) \propto \exp(-U(x)/D)},
#' where \eqn{D} is a noise intensity parameter. Setting \eqn{D = 1}
#' gives \eqn{U(x) = -\log P(x)}. The probability density \eqn{P(x)} is
#' estimated via `stats::density()` with a Gaussian kernel.
#'
#' **Wells** (local minima of \eqn{U}) correspond to stable attractors --
#' states the system tends to occupy. **Barriers** (local maxima of
#' \eqn{U}) correspond to unstable equilibria separating attractors.
#'
#' Well depth measures how much "energy" is needed to escape an attractor:
#' deeper wells are more stable. Well width measures the basin of
#' attraction in state space. Together, depth and width characterise the
#' resilience of each attractor -- deep, wide wells indicate robust
#' stability, while shallow, narrow wells are easily escaped by
#' perturbations.
#'
#' **Detrending.** Non-stationary trends can bias the density estimate and
#' create spurious features. The `detrend` parameter offers three options:
#' \describe{
#'   \item{`"none"` (default)}{No detrending. Appropriate for stationary
#'     or nearly stationary data.}
#'   \item{`"linear"`}{Removes a linear trend via OLS before estimating
#'     the landscape. Appropriate when the data has a slow drift.}
#'   \item{`"diff"`}{Uses first differences \eqn{\Delta x_t = x_t -
#'     x_{t-1}}. Appropriate for integrated processes (random walks).
#'     Note: reduces series length by one.}
#' }
#'
#' **Rolling-window mode.** When `window` is specified, the landscape is
#' re-estimated at each window position and the number of wells and the
#' location of the dominant (deepest) well are tracked over time. Changes
#' in the number of wells -- especially transitions from one to two wells
#' or vice versa -- can signal an approaching bifurcation (tipping point).
#'
#' @export
#' @param data \[`numeric` or `ts`\]\cr
#'   The time series to analyze. Accepts a numeric vector or a `ts` object.
#' @param window \[`integer(1)`: `NULL`\]\cr
#'   Rolling window size. If `NULL` (default), a single global landscape
#'   is estimated from the entire series. If specified, the landscape is
#'   re-estimated for each window position.
#' @param n_bins \[`integer(1)`: `50L`\]\cr
#'   Number of bins for the kernel density estimate grid. The actual grid
#'   resolution passed to `stats::density()` is `n_bins * 4` to ensure
#'   smooth interpolation.
#' @param bandwidth \[`numeric(1)`: `NULL`\]\cr
#'   Bandwidth for the kernel density estimate. If `NULL`, Silverman's
#'   rule of thumb is used (the `stats::density()` default).
#' @param detrend \[`character(1)`: `"none"`\]\cr
#'   Detrending method applied before estimating the landscape. One of
#'   `"none"`, `"linear"`, or `"diff"`.
#'
#' @return An object of class `"potential"` (a list) containing:
#'   \describe{
#'     \item{landscape}{\[`tibble`\] Grid of the estimated landscape with
#'       columns `x` (state-space position), `density` (probability
#'       density), and `potential` (\eqn{-\log(\mathrm{density})}).
#'       Present only for global analysis (`window = NULL`).}
#'     \item{wells}{\[`tibble`\] One row per detected well with columns
#'       `location` (x position), `depth` (barrier height minus well
#'       bottom), and `width` (distance between adjacent barriers).
#'       Present only for global analysis.}
#'     \item{barriers}{\[`tibble`\] One row per detected barrier with
#'       columns `location` (x position) and `height` (potential value).
#'       Present only for global analysis.}
#'     \item{n_wells}{\[`integer(1)`\] Number of detected wells.}
#'     \item{values}{\[`numeric()`\] The (possibly detrended) data used
#'       for the analysis.}
#'     \item{time}{\[`numeric()`\] The time index.}
#'     \item{rolling}{\[`tibble`\] Present only for rolling analysis
#'       (`window` specified). Contains columns `time`, `value`,
#'       `n_wells`, and `dominant_well_location`.}
#'   }
#'
#'   Attributes: `n_bins`, `bandwidth`, `detrend`, `window`.
#'
#' @references
#' Livina, V.N., Kwasniok, F., & Lenton, T.M. (2010). Potential analysis
#' reveals changing number of climate states during the last 60 kyr.
#' \emph{Climate of the Past}, 6(1), 77--82.
#' \doi{10.5194/cp-6-77-2010}
#'
#' Scheffer, M., Carpenter, S.R., Lenton, T.M., Bascompte, J.,
#' Brock, W., Dakos, V., van de Koppel, J., van de Leemput, I.A.,
#' Levin, S.A., van Nes, E.H., Pascual, M., & Vandermeer, J. (2012).
#' Anticipating critical transitions.
#' \emph{Science}, 338(6105), 344--348.
#' \doi{10.1126/science.1225244}
#'
#' Strogatz, S.H. (2015). \emph{Nonlinear Dynamics and Chaos} (2nd ed.).
#' Westview Press.
#'
#' @seealso [resilience()] for rolling resilience metrics;
#'   [hurst()] for long-range dependence analysis;
#'   [compute_trend()] for trend classification.
#' @family potential
#' @concept stability landscape
#' @concept potential analysis
#' @concept tipping points
#' @examples
#' \donttest{
#' # Single well: Ornstein-Uhlenbeck process (attractor at 0)
#' set.seed(42)
#' n <- 2000
#' x <- numeric(n)
#' for (i in 2:n) x[i] <- x[i-1] - 0.5 * x[i-1] * 0.01 + 0.3 * rnorm(1)
#' pa <- potential_analysis(x)
#' pa
#' plot(pa)
#'
#' # Two wells: bimodal mixture
#' set.seed(42)
#' bimodal <- c(rnorm(1000, -2, 0.5), rnorm(1000, 2, 0.5))
#' pa2 <- potential_analysis(bimodal)
#' summary(pa2)
#' plot(pa2, type = "both")
#' }
potential_analysis <- function(data,
                               window = NULL,
                               n_bins = 50L,
                               bandwidth = NULL,
                               detrend = "none") {
  # --- 1. Input validation ---
  check_missing(data)
  data <- prepare_timeseries_data(data)
  values <- data$values
  time <- data$time
  n <- length(values)

  detrend <- check_match(detrend, c("none", "linear", "diff"))
  check_range(n_bins, type = "integer", min = 10L, max = 10000L)
  if (!is.null(bandwidth)) {
    check_range(bandwidth, type = "numeric", scalar = TRUE, min = 0)
  }

  # --- 2. Detrending ---
  detrended <- potential_detrend_(values, detrend)
  if (detrend == "diff") {
    # First differencing loses one observation
    time <- time[-1L]
  }
  n_d <- length(detrended)

  # --- 3. Validate window if specified ---
  if (!is.null(window)) {
    check_range(window, type = "integer", min = 30L, max = n_d)
  }

  # --- 4. Analysis ---
  if (is.null(window)) {
    # Global analysis
    result <- potential_landscape_(detrended, n_bins, bandwidth)

    out <- list(
      landscape = result$landscape,
      wells = result$wells,
      barriers = result$barriers,
      n_wells = result$n_wells,
      values = detrended,
      time = time
    )
  } else {
    # Rolling-window analysis
    n_positions <- n_d - window + 1L
    rolling_n_wells <- integer(n_d)
    rolling_dominant <- rep(NA_real_, n_d)

    # Pre-fill with NA for positions without full window
    rolling_n_wells[seq_len(window - 1L)] <- NA_integer_

    for (i in seq_len(n_positions)) {
      end_idx <- i + window - 1L
      w_data <- detrended[i:end_idx]
      w_clean <- w_data[!is.na(w_data)]
      if (length(w_clean) < 20L) {
        rolling_n_wells[end_idx] <- NA_integer_
        next
      }
      w_result <- try_(potential_landscape_(w_clean, n_bins, bandwidth))
      if (inherits(w_result, "try-error") || is.null(w_result)) {
        rolling_n_wells[end_idx] <- NA_integer_
        next
      }
      rolling_n_wells[end_idx] <- w_result$n_wells
      # Dominant well = deepest well
      if (nrow(w_result$wells) > 0L) {
        best <- which.max(w_result$wells$depth)
        rolling_dominant[end_idx] <- w_result$wells$location[best]
      }
    }

    rolling_tbl <- tibble::tibble(
      time = time,
      value = detrended,
      n_wells = rolling_n_wells,
      dominant_well_location = rolling_dominant
    )

    # Also compute a global landscape for reference
    global_result <- potential_landscape_(detrended, n_bins, bandwidth)

    out <- list(
      landscape = global_result$landscape,
      wells = global_result$wells,
      barriers = global_result$barriers,
      n_wells = global_result$n_wells,
      values = detrended,
      time = time,
      rolling = rolling_tbl
    )
  }

  structure(
    out,
    n_bins = n_bins,
    bandwidth = bandwidth,
    detrend = detrend,
    window = window,
    class = c("potential", "list")
  )
}

# --------------------------------------------------------------------------
# S3 methods
# --------------------------------------------------------------------------

#' Plot Potential Analysis Results
#'
#' Visualizes the estimated potential landscape, probability density, or
#' both. Wells are marked as points and barriers as vertical dashed lines.
#'
#' @details
#' Three plot types are available:
#'
#' **`type = "landscape"` (default).** Plots the quasi-potential
#' \eqn{U(x) = -\log P(x)} as a function of state-space position
#' \eqn{x}. Wells (local minima) are shown as blue points and barriers
#' (local maxima) as red vertical dashed lines. This view directly shows
#' the stability landscape: valleys are stable attractors, peaks are
#' unstable equilibria, and the depth of each valley relative to the
#' nearest peak quantifies how resilient each attractor is to
#' perturbations.
#'
#' **`type = "density"`.** Plots the kernel density estimate \eqn{P(x)}.
#' Peaks in density correspond to wells in potential and vice versa. This
#' view is useful for checking the quality of the density estimate and
#' confirming that the identified attractors correspond to genuine modes
#' in the data distribution.
#'
#' **`type = "both"`.** Stacks the density plot on top and the potential
#' landscape below using [patchwork::wrap_plots()]. The shared x-axis
#' makes it easy to see how features in density map to features in
#' potential.
#'
#' When the object contains rolling-window results (i.e., `window` was
#' specified), an additional panel showing the number of wells over time
#' is appended below the landscape panel(s).
#'
#' @export
#' @param x \[`potential`\]\cr
#'   An object of class `potential` as returned by [potential_analysis()].
#' @param type \[`character(1)`: `"landscape"`\]\cr
#'   Plot type: `"landscape"`, `"density"`, or `"both"`.
#' @param ... Additional arguments (currently unused).
#' @return A [ggplot2::ggplot()] object (or a
#'   [patchwork::wrap_plots()] composite when `type = "both"` or when
#'   rolling results are present).
#'
#' @seealso [potential_analysis()] for computing the potential landscape.
#' @family potential
#' @concept stability landscape
#' @examples
#' \donttest{
#' set.seed(42)
#' bimodal <- c(rnorm(1000, -2, 0.5), rnorm(1000, 2, 0.5))
#' pa <- potential_analysis(bimodal)
#' plot(pa, type = "landscape")
#' plot(pa, type = "density")
#' plot(pa, type = "both")
#' }
plot.potential <- function(x, type = "landscape", ...) {
  check_missing(x)
  check_class(x, "potential")
  type <- check_match(type, c("landscape", "density", "both"))

  colors <- potential_colors_()
  landscape <- x$landscape
  wells <- x$wells
  barriers <- x$barriers

  panels <- list()

  # --- Density panel ---
  if (type == "density" || type == "both") {
    p_density <- ggplot2::ggplot(
      landscape,
      ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("density"))
    ) +
      ggplot2::geom_line(
        linewidth = 0.7,
        color = colors["density"]
      ) +
      ggplot2::geom_area(
        alpha = 0.15,
        fill = colors["density"]
      )

    if (nrow(barriers) > 0L) {
      p_density <- p_density +
        ggplot2::geom_vline(
          xintercept = barriers$location,
          linetype = "dashed",
          color = colors["barrier"],
          linewidth = 0.5
        )
    }

    p_density <- p_density +
      ggplot2::labs(
        title = "Probability Density Estimate",
        x = ifelse_(type == "both", NULL, "State (x)"),
        y = "Density P(x)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.title = ggplot2::element_text(color = "black", face = "bold"),
        axis.text = ggplot2::element_text(color = "black")
      )

    if (type == "both") {
      p_density <- p_density +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    }

    if (type == "density") return(p_density)
    panels <- c(panels, list(p_density))
  }

  # --- Potential landscape panel ---
  if (type == "landscape" || type == "both") {
    p_potential <- ggplot2::ggplot(
      landscape,
      ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("potential"))
    ) +
      ggplot2::geom_line(
        linewidth = 0.7,
        color = colors["potential"]
      )

    # Mark wells
    if (nrow(wells) > 0L) {
      # Find potential values at well locations via interpolation
      well_potential <- stats::approx(
        landscape$x, landscape$potential,
        xout = wells$location, rule = 2
      )$y

      well_df <- data.frame(
        x = wells$location,
        y = well_potential
      )

      p_potential <- p_potential +
        ggplot2::geom_point(
          data = well_df,
          ggplot2::aes(
            x = !!rlang::sym("x"),
            y = !!rlang::sym("y")
          ),
          color = colors["well"],
          size = 3,
          shape = 19,
          inherit.aes = FALSE
        )
    }

    # Mark barriers
    if (nrow(barriers) > 0L) {
      p_potential <- p_potential +
        ggplot2::geom_vline(
          xintercept = barriers$location,
          linetype = "dashed",
          color = colors["barrier"],
          linewidth = 0.5
        )
    }

    p_potential <- p_potential +
      ggplot2::labs(
        title = ifelse_(
          type == "both",
          "Potential Landscape",
          "Stability Landscape: Potential U(x) = -log P(x)"
        ),
        x = "State (x)",
        y = "Potential U(x)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.title = ggplot2::element_text(color = "black", face = "bold"),
        axis.text = ggplot2::element_text(color = "black")
      )

    if (type == "landscape" && is.null(x$rolling)) return(p_potential)
    panels <- c(panels, list(p_potential))
  }

  # --- Rolling panel (if present) ---
  if (!is.null(x$rolling)) {
    rolling <- x$rolling
    p_rolling <- ggplot2::ggplot(
      rolling,
      ggplot2::aes(
        x = !!rlang::sym("time"),
        y = !!rlang::sym("n_wells")
      )
    ) +
      ggplot2::geom_step(linewidth = 0.6, color = "#2196F3") +
      ggplot2::scale_y_continuous(
        breaks = function(lims) {
          seq(floor(lims[1L]), ceiling(lims[2L]), by = 1L)
        }
      ) +
      ggplot2::labs(
        title = "Number of Wells Over Time",
        x = "Time",
        y = "Number of Wells"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.title = ggplot2::element_text(color = "black", face = "bold"),
        axis.text = ggplot2::element_text(color = "black")
      )
    panels <- c(panels, list(p_rolling))
  }

  if (length(panels) == 1L) return(panels[[1L]])

  patchwork::wrap_plots(panels, ncol = 1L)
}

#' Print a Potential Analysis Object
#'
#' Prints a concise summary of the potential analysis: the number of
#' wells and barriers, well locations, and analysis settings.
#'
#' @describeIn potential_analysis Print method for `"potential"` objects.
#'
#' @export
#' @param x \[`potential`\]\cr
#'   A `potential` object.
#' @param ... Additional arguments (currently unused).
#' @return `x`, invisibly.
print.potential <- function(x, ...) {
  cat("Potential Analysis\n")
  cat("  Detrend  :", attr(x, "detrend"), "\n")
  cat("  N points :", length(x$values), "\n")
  cat("  N bins   :", attr(x, "n_bins"), "\n")
  bw_str <- ifelse_(
    is.null(attr(x, "bandwidth")),
    "auto (Silverman)",
    as.character(round(attr(x, "bandwidth"), 4))
  )
  cat("  Bandwidth:", bw_str, "\n")
  win_str <- ifelse_(
    is.null(attr(x, "window")),
    "global",
    as.character(attr(x, "window"))
  )
  cat("  Window   :", win_str, "\n")
  cat("  Wells    :", x$n_wells, "\n")
  cat("  Barriers :", nrow(x$barriers), "\n")

  if (x$n_wells > 0L) {
    cat("\n  Well locations:\n")
    for (i in seq_len(nrow(x$wells))) {
      cat(sprintf(
        "    Well %d: x = %.4f  (depth = %.4f, width = %.4f)\n",
        i, x$wells$location[i], x$wells$depth[i], x$wells$width[i]
      ))
    }
  }

  if (!is.null(x$rolling)) {
    n_wells_range <- range(x$rolling$n_wells, na.rm = TRUE)
    cat(sprintf(
      "\n  Rolling: n_wells range = [%d, %d]\n",
      n_wells_range[1L], n_wells_range[2L]
    ))
  }

  invisible(x)
}

#' Summarize Potential Analysis Results
#'
#' Prints detailed information about each well and barrier in the
#' estimated potential landscape, including location, depth, width,
#' and height.
#'
#' @describeIn potential_analysis Summary method for `"potential"` objects.
#'
#' @export
#' @param object \[`potential`\]\cr
#'   A `potential` object.
#' @param ... Additional arguments (currently unused).
#' @return A list with elements `wells`, `barriers`, `n_wells`,
#'   `detrend`, `n_bins`, `bandwidth`, and `window`, returned invisibly.
summary.potential <- function(object, ...) {
  check_missing(object)
  check_class(object, "potential")

  cat("Potential Analysis Summary\n")
  cat(paste0(rep("=", 50), collapse = ""), "\n")
  cat("  Series length  :", length(object$values), "\n")
  cat("  Detrending     :", attr(object, "detrend"), "\n")
  cat("  Grid resolution:", attr(object, "n_bins") * 4L, "points\n")
  bw_str <- ifelse_(
    is.null(attr(object, "bandwidth")),
    "auto (Silverman)",
    as.character(round(attr(object, "bandwidth"), 4))
  )
  cat("  Bandwidth      :", bw_str, "\n")

  cat("\n  Landscape topology\n")
  cat(paste0("  ", rep("-", 40), collapse = ""), "\n")
  cat("  Number of wells   :", object$n_wells, "\n")
  cat("  Number of barriers:", nrow(object$barriers), "\n")

  if (object$n_wells > 0L) {
    cat("\n  Wells:\n")
    for (i in seq_len(nrow(object$wells))) {
      cat(sprintf(
        "    [%d] location = %8.4f | depth = %8.4f | width = %8.4f\n",
        i, object$wells$location[i],
        object$wells$depth[i], object$wells$width[i]
      ))
    }
  }

  if (nrow(object$barriers) > 0L) {
    cat("\n  Barriers:\n")
    for (i in seq_len(nrow(object$barriers))) {
      cat(sprintf(
        "    [%d] location = %8.4f | height = %8.4f\n",
        i, object$barriers$location[i], object$barriers$height[i]
      ))
    }
  }

  if (!is.null(object$rolling)) {
    rolling <- object$rolling
    valid_wells <- rolling$n_wells[!is.na(rolling$n_wells)]
    cat("\n  Rolling-window analysis\n")
    cat(paste0("  ", rep("-", 40), collapse = ""), "\n")
    cat("  Window size        :", attr(object, "window"), "\n")
    cat("  Min wells observed :", min(valid_wells), "\n")
    cat("  Max wells observed :", max(valid_wells), "\n")
    cat("  Mean wells         :", round(mean(valid_wells), 2), "\n")

    # Count transitions in well number
    well_diff <- diff(valid_wells)
    n_transitions <- sum(well_diff != 0L, na.rm = TRUE)
    cat("  Well-count changes :", n_transitions, "\n")
  }

  invisible(list(
    wells = object$wells,
    barriers = object$barriers,
    n_wells = object$n_wells,
    detrend = attr(object, "detrend"),
    n_bins = attr(object, "n_bins"),
    bandwidth = attr(object, "bandwidth"),
    window = attr(object, "window")
  ))
}
