# ============================================================================
# to_tna.R -- Extract state ribbons as wide sequence data
# Converts categorical state trajectories from analysis objects into
# wide-format tibbles: one row per sequence, columns T1, T2, T3, ...
# ============================================================================

# --------------------------------------------------------------------------
# Internal helper
# --------------------------------------------------------------------------

#' Convert a state vector to a single-row wide tibble
#' @param states Character vector of state labels (NAs removed).
#' @param label  Context label for error messages.
#' @return A single-row tibble with columns T1, T2, ..., Tn.
#' @noRd
seq_to_wide_ <- function(states, label = "column") {
  states <- states[!is.na(states)]
  if (length(states) == 0L) {
    stop_("The {.val {label}} contains only NA values -- no states to extract.")
  }
  names(states) <- paste0("T", seq_along(states))
  tibble::as_tibble_row(states)
}

# --------------------------------------------------------------------------
# Generic
# --------------------------------------------------------------------------

#' Extract State Sequences as Wide-Format Data
#'
#' Extracts the categorical state trajectory from an analysis object and
#' returns it as a wide-format tibble with one row per sequence and columns
#' `T1, T2, T3, ...` containing the state at each time point. The output
#' is directly compatible with `tna::prepare_data()` for transition network
#' analysis, `TraMineR::seqdef()` for sequence analysis, or any method
#' that expects wide-format categorical sequences.
#'
#' @param x An analysis object (resilience, hurst, hurst_ews, multi_ews, trend,
#'   changepoint, spectral, or potential).
#' @param ... Additional arguments passed to methods.
#' @return A tibble with one row per sequence. Column names are `T1`, `T2`,
#'   etc. All values are character.
#'
#' @details
#' Each analysis in the pipeline produces a categorical state ribbon over time:
#'
#' \describe{
#'   \item{resilience}{Six-level resilience categories (Excellent to Troubled)
#'     from `classify_resilience()`. Available for the composite score and
#'     each individual metric.}
#'   \item{hurst}{Five Hurst memory states (strong_antipersistent to
#'     strong_persistent), or four multifractal categories when using MFDFA.}
#'   \item{hurst_ews}{Five warning levels (none to critical) from
#'     `detect_hurst_warnings()`.}
#'   \item{multi_ews}{System states (Stable to Failing) from expanding-window
#'     multivariate EWS. Rolling-window objects have no state classification
#'     and are not supported.}
#'   \item{trend}{Four trend states (ascending, descending, flat, turbulent)
#'     from `compute_trend()`.}
#'   \item{changepoint}{Regime state labels combining level and changepoint
#'     type (e.g., `high_change`, `low_return`, `medium_initial`) from
#'     `detect_changepoints()`. The `state` argument selects which column
#'     to extract: `"state"` (default), `"level"`, `"direction"`,
#'     `"regime"`, or `"segment"`.}
#'   \item{spectral}{Four noise-colour states (white_noise, pink_noise,
#'     red_noise, brownian) from `spectral_ews()` with `states = TRUE`.}
#'   \item{potential}{Stability states (unimodal, bimodal, multimodal) derived
#'     from the rolling well count in `potential_analysis()`. Requires
#'     rolling-window mode (`window` specified).}
#' }
#'
#' @examples
#' \donttest{
#' # Resilience state sequences
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' res <- resilience(x, window = 50L)
#' cls <- classify_resilience(res)
#' to_tna(cls)                    # composite categories
#' to_tna(cls, state = "vsi")     # per-metric categories
#'
#' # Hurst memory state sequences
#' h <- hurst(x, method = "dfa", window = 50L, step = 10L)
#' to_tna(h)
#'
#' # EWS warning level sequences
#' ews <- detect_hurst_warnings(h)
#' to_tna(ews)
#' to_tna(ews, state = "state")   # underlying Hurst states
#'
#' # Changepoint segment sequences
#' cp <- detect_changepoints(x, method = "pelt")
#' to_tna(cp)
#'
#' # Spectral noise-colour sequences
#' sp <- spectral_ews(x, window = 50L, states = TRUE)
#' to_tna(sp)
#'
#' # Potential stability state sequences (rolling mode)
#' pot <- potential_analysis(x, window = 100L)
#' to_tna(pot)
#' }
#' @export
to_tna <- function(x, ...) {
  UseMethod("to_tna")
}

# --------------------------------------------------------------------------
# Default: informative error for unsupported objects
# --------------------------------------------------------------------------

#' @rdname to_tna
#' @export
to_tna.default <- function(x, ...) {
  cls <- paste(class(x), collapse = "/")
  supported <- paste(
    "resilience (classified), hurst, hurst_ews, multi_ews (expanding),",
    "trend, changepoint, spectral, potential (rolling)"
  )
  stop_(
    "{.fun to_tna} does not support objects of class {.cls {cls}}.
    Supported: {supported}."
  )
}

# --------------------------------------------------------------------------
# Method: hurst_global --explicit rejection
# --------------------------------------------------------------------------

#' @rdname to_tna
#' @export
to_tna.hurst_global <- function(x, ...) {
  stop_(
    "{.fun to_tna} requires a rolling {.cls hurst} object with states.
    This is a global (single-value) estimate --
    rerun {.fun hurst} with {.code states = TRUE} to get a state sequence."
  )
}

# --------------------------------------------------------------------------
# Method: resilience
# --------------------------------------------------------------------------

#' @rdname to_tna
#' @param state Which state column to extract. `"composite"` (default)
#'   extracts `composite_category`; any other value extracts
#'   `{state}_category` (e.g., `"vsi"` extracts `vsi_category`).
#' @export
to_tna.resilience <- function(x, state = "composite", ...) {
  stopifnot_(
    is.character(state) && length(state) == 1L,
    "{.arg state} must be a single {.cls character} value."
  )

  # Must have been classified first
  if (!any(grepl("_category$", names(x)))) {
    stop_(
      "No state categories found in this {.cls resilience} object.
      Run {.fun classify_resilience} first to create state labels."
    )
  }

  col_name <- paste0(state, "_category")
  if (!col_name %in% names(x)) {
    avail <- grep("_category$", names(x), value = TRUE)
    avail <- sub("_category$", "", avail)
    stop_(
      "State {.val {state}} not found (looked for column {.val {col_name}}).
      Available: {.or {.val {avail}}}."
    )
  }

  states <- as.character(x[[col_name]])
  seq_to_wide_(states, label = col_name)
}

# --------------------------------------------------------------------------
# Method: hurst
# --------------------------------------------------------------------------

#' @rdname to_tna
#' @param state Which state column to extract. `"state"` (default) extracts
#'   the Hurst state classification; `"mf_category"` extracts the
#'   multifractal category (MFDFA only).
#' @export
to_tna.hurst <- function(x, state = "state", ...) {
  stopifnot_(
    is.character(state) && length(state) == 1L,
    "{.arg state} must be a single {.cls character} value."
  )

  valid <- intersect(c("state", "mf_category"), names(x))
  if (length(valid) == 0L) {
    stop_(
      "This {.cls hurst} object has no state columns.
      Rerun {.fun hurst} with {.code states = TRUE}."
    )
  }

  if (!state %in% names(x)) {
    stop_(
      "State {.val {state}} not found. Available: {.or {.val {valid}}}."
    )
  }

  states <- as.character(x[[state]])
  seq_to_wide_(states, label = state)
}

# --------------------------------------------------------------------------
# Method: hurst_ews
# --------------------------------------------------------------------------

#' @rdname to_tna
#' @param state Which state column to extract. `"warning_label"` (default)
#'   extracts the EWS warning level labels; `"state"` extracts the
#'   underlying Hurst state.
#' @export
to_tna.hurst_ews <- function(x, state = "warning_label", ...) {
  stopifnot_(
    is.character(state) && length(state) == 1L,
    "{.arg state} must be a single {.cls character} value."
  )

  valid <- intersect(c("warning_label", "state"), names(x))
  if (length(valid) == 0L) {
    stop_("This {.cls hurst_ews} object has no state columns.")
  }

  if (!state %in% names(x)) {
    stop_(
      "State {.val {state}} not found. Available: {.or {.val {valid}}}."
    )
  }

  states <- as.character(x[[state]])
  seq_to_wide_(states, label = state)
}

# --------------------------------------------------------------------------
# Method: multi_ews (expanding only)
# --------------------------------------------------------------------------

#' @rdname to_tna
#' @export
to_tna.multi_ews <- function(x, ...) {
  method <- attr(x, "method")

  if (is.null(method) || method != "expanding") {
    stop_(
      "{.fun to_tna} requires an expanding-window {.cls multi_ews} object.
      Rolling-window objects do not produce state classifications --
      rerun {.fun detect_multivariate_warnings} with
      {.code method = \"expanding\"}."
    )
  }

  cls <- attr(x, "classification")
  if (is.null(cls) || !"state" %in% names(cls)) {
    stop_(
      "No classification found in this expanding {.cls multi_ews} object.
      The analysis produced no state labels."
    )
  }

  states <- as.character(cls$state)
  seq_to_wide_(states, label = "classification$state")
}

# --------------------------------------------------------------------------
# Method: trend
# --------------------------------------------------------------------------

#' @rdname to_tna
#' @export
to_tna.trend <- function(x, ...) {
  if (!"state" %in% names(x)) {
    stop_("This {.cls trend} object has no {.val state} column.")
  }
  states <- as.character(x$state)
  seq_to_wide_(states, label = "state")
}

# --------------------------------------------------------------------------
# Method: changepoint
# --------------------------------------------------------------------------

#' @rdname to_tna
#' @param state Which column to extract as the TNA sequence. `"state"`
#'   (default) extracts the combined `level_changepoint_type` label.
#'   Also accepts `"level"`, `"direction"`, `"regime"`, or `"segment"`.
#' @export
to_tna.changepoint <- function(x, state = "state", ...) {
  stopifnot_(
    is.character(state) && length(state) == 1L,
    "{.arg state} must be a single {.cls character} value."
  )

  # Map "regime" and "segment" to formatted string columns
  if (state == "regime") {
    states <- paste0("regime_", x$regime)
  } else if (state == "segment") {
    states <- paste0("segment_", x$segment)
  } else {
    valid <- c("state", "level", "direction", "regime", "segment")
    if (!state %in% names(x)) {
      stop_("Column {.val {state}} not found. Available: {.or {.val {valid}}}.")
    }
    states <- as.character(x[[state]])
  }
  seq_to_wide_(states, label = state)
}

# --------------------------------------------------------------------------
# Method: spectral
# --------------------------------------------------------------------------

#' @rdname to_tna
#' @export
to_tna.spectral <- function(x, ...) {
  if (!"state" %in% names(x)) {
    stop_(
      "This {.cls spectral} object has no {.val state} column.
      Rerun {.fun spectral_ews} with {.code states = TRUE}."
    )
  }
  states <- as.character(x$state)
  seq_to_wide_(states, label = "state")
}

# --------------------------------------------------------------------------
# Method: potential (rolling-window only)
# --------------------------------------------------------------------------

#' @rdname to_tna
#' @export
to_tna.potential <- function(x, ...) {
  if (is.null(x$rolling)) {
    stop_(
      "{.fun to_tna} requires a rolling-window {.cls potential} object.
      This is a global (single-landscape) analysis --
      rerun {.fun potential_analysis} with a {.arg window} argument
      to get a stability state sequence over time."
    )
  }

  n_wells <- x$rolling$n_wells
  # Classify well counts into stability states
  states <- vapply(n_wells, function(nw) {
    if (is.na(nw)) return(NA_character_)
    if (nw == 0L) return("flat")
    if (nw == 1L) return("unimodal")
    if (nw == 2L) return("bimodal")
    "multimodal"
  }, character(1L))

  seq_to_wide_(states, label = "n_wells")
}

# --------------------------------------------------------------------------
# Method: surrogate_test -- explicit rejection
# --------------------------------------------------------------------------

#' @rdname to_tna
#' @export
to_tna.surrogate_test <- function(x, ...) {
  stop_(
    "{.fun to_tna} does not support {.cls surrogate_test} objects.
    Surrogate tests produce significance statistics,
    not state sequences over time."
  )
}

# --------------------------------------------------------------------------
# Method: sensitivity_ews -- explicit rejection
# --------------------------------------------------------------------------

#' @rdname to_tna
#' @export
to_tna.sensitivity_ews <- function(x, ...) {
  stop_(
    "{.fun to_tna} does not support {.cls sensitivity_ews} objects.
    Sensitivity analyses produce robustness metrics across
    parameter combinations, not state sequences over time."
  )
}
