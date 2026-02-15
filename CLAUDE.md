# codyna R Package Conventions

All R files in `claude/R/` follow these conventions. Any new function added to this package MUST conform to every rule below.

## Input Handling

- The first parameter for univariate functions is always `data` (never `original_series`, `x`, `series`, etc.)
- Parse it immediately with `prepare_timeseries_data(data)` which returns `list(values, time)` and handles both numeric vectors and `ts` objects
- For multivariate functions taking a data.frame, the first parameter is still `data`

## Parameter Naming

| Concept | Name | NOT |
|---|---|---|
| Window size | `window` | `window_size`, `win_size`, `winsize` |
| Window alignment | `align` | `window_alignment`, `alignment` |
| Minimum observations | `min_points` | `min_data_points`, `min_obs` |
| Analysis method | `method` | `type`, `approach` |
| State output | `states` (logical flag) | `classify`, `labels` |

## Validation

Always validate in this order at the top of every exported function:

```r
check_missing(data)
data <- prepare_timeseries_data(data)
method <- check_match(method, c("option_a", "option_b"))
align <- check_match(align, c("center", "right", "left"))
check_range(window, type = "integer", min = 2L, max = n)
check_flag(states)
```

- Use `check_missing()`, `check_match()`, `check_range()`, `check_flag()`, `check_class()`
- Use `stop_()`, `stopifnot_()`, `warning_()`, `message_()` for messages (cli-formatted)
- NEVER use base `stop()`, `stopifnot()`, `warning()`, `message()` in exported functions
- Use `ifelse_()` instead of `ifelse()` for scalar conditions
- Use `try_()` instead of `tryCatch()` for simple error suppression

## Return Objects

- Always return a **tibble** via `tibble::tibble()` or `tibble::as_tibble()`
- Wrap with `structure()` to set a custom S3 class and store settings as attributes:

```r
structure(
  out,
  window = window,
  align  = align,
  method = method,
  class  = c("trend", "tbl_df", "tbl", "data.frame")
)
```

- Class vector: `c("{class_name}", "tbl_df", "tbl", "data.frame")`
- Store ALL analysis settings as attributes so they can be retrieved later

## Output Columns

- First column: `time` (the time index)
- Second column: `value` (the original series values)
- Remaining columns: computed metrics/results
- Classification column: `state` (factor with ordered levels)
- NEVER use `time_index`, `original_data`, `trend_codes`, `metric_values`

## S3 Methods

Every new S3 class MUST have these three methods:

```r
# print — delegates to tibble
print.xxx <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

# plot — returns a ggplot2 object
plot.xxx <- function(x, ...) {
  check_missing(x)
  check_class(x, "xxx")
  # ... build ggplot ...
}

# summary — prints key info, returns list invisibly
summary.xxx <- function(object, ...) {
  check_missing(object)
  check_class(object, "xxx")
  # ... cat() summary ...
  invisible(list(...))
}
```

## Internal Helpers

- Name pattern: `{topic}_{action}_()` with trailing underscore
- Examples: `trend_ar1_()`, `trend_slope_()`, `hurst_classify_()`, `mews_theme_()`
- Always mark with `@noRd`
- NEVER use `.dotPrefix` naming (not `.calculate_ar1_phi1_single_window`)

## File Structure

```r
# ============================================================================
# {filename}.R — {one-line description}
# {optional context line}
# ============================================================================

# --------------------------------------------------------------------------
# Internal helpers (@noRd)
# --------------------------------------------------------------------------

#' {brief title}
#' @noRd
helper_function_ <- function(...) { }

# --------------------------------------------------------------------------
# Exported: {function_name}()
# --------------------------------------------------------------------------

#' {Title}
#' @export
exported_function <- function(...) { }

# --------------------------------------------------------------------------
# S3 methods
# --------------------------------------------------------------------------

#' @export
plot.class <- function(x, ...) { }

#' @describeIn exported_function Print method.
#' @export
print.class <- function(x, ...) { }

#' @describeIn exported_function Summary method.
#' @export
summary.class <- function(object, ...) { }
```

## Roxygen2 Documentation

- Parameter format: `@param name \[type: default\]\cr` followed by description
  ```
  #' @param window \[`integer(1)`: `50L`\]\cr
  #'   Rolling window size in observations.
  ```
- Always include: `@export`, `@param` (all params), `@return`, `@examples`
- Examples wrap in `\donttest{}` for long-running code
- Use `@family {topic}` to group related functions
- Use `@concept {concept}` for cross-topic discovery
- Use `@seealso` to link to directly related functions
- Use `@describeIn` for S3 methods to group them with the main function
- Use `@noRd` for internal helpers (never `@keywords internal`)

## ggplot2 Plots

- Always use `!!rlang::sym("column_name")` for column references in `aes()`
- Use `ggplot2::` namespace prefix for every ggplot2 function
- Use `patchwork::wrap_plots()` for multi-panel layouts
- Store color palettes in `{topic}_colors_()` internal helpers
- Return the ggplot object (do not print it)

## Dependencies

- Only packages already in codyna DESCRIPTION: checkmate, cli, dplyr, ggplot2, patchwork, rlang, scales, stats, tibble, tidyr, tidyselect
- NEVER add new dependencies without explicit approval
- Always use `pkg::function()` namespace syntax (never `library()` or `@import`)

## Existing Classes

| Class | Created by | Columns |
|---|---|---|
| `resilience` | `resilience()` | time, value, vsi, arch_lm, cv, recovery_time, recovery_slope, sample_entropy, dfa_alpha, ac_ratio |
| `resilience` (classified) | `classify_resilience()` | + {metric}_score, {metric}_category, composite_score, composite_category |
| `hurst` | `hurst()` | time, value, hurst, r_squared, state, transition |
| `hurst_global` | `hurst(states=FALSE)` | list: hurst, r_squared, method, n |
| `hurst_ews` | `detect_hurst_warnings()` | time, value, hurst, r_squared, state, transition, {indicator cols}, warning_score, warning_level, warning_label |
| `multi_ews` (rolling) | `detect_multivariate_warnings()` | time, metric, score, std |
| `multi_ews` (expanding) | `detect_multivariate_warnings()` | time, metric, score, z_score, detected; attr: classification(time, count, state) |
| `trend` | `compute_trend()` | time, value, metric, state |
| `changepoint` | `detect_changepoints()` | time, value, segment, regime, segment_mean, segment_var, level, direction, magnitude, changepoint_type, state, changepoint |
| `spectral` | `spectral_ews()` | time, value, spectral_exponent, spectral_ratio, r_squared, state |
| `potential` | `potential_analysis()` | list: values, time, n_wells, landscape, wells, barriers, rolling |
| `surrogate_test` | `surrogate_test()` | list: observed_tau, surrogate_taus, p_value, significant, ... |
| `sensitivity_ews` | `sensitivity_ews()` | window, detrend, time, score, tau |

## CRAN Conformance Rules

These rules were discovered during `R CMD check --as-cran` and MUST be followed in all code.

### ASCII-Only R Code

- R source files must contain **only ASCII characters** (no Unicode in code or string literals)
- Em dashes: use `--` not `---` (U+2014)
- Greek letters: use the English word (`tau`, `alpha`, `sigma`) not Unicode (`\u03C4`, `\u03B1`, `\u03C3`)
- Superscripts: use `R^2` not `\u00b2`
- Comments are the only place non-ASCII is tolerated, but avoid it everywhere
- Use `tools::showNonASCIIfile("R/file.R")` to audit a file

### Plot Labels

- **Never use Unicode characters in ggplot2 labels, titles, subtitles, or captions**
- `grid.Call(C_textBounds)` crashes on non-ASCII text on some platforms
- Write `"Kendall's tau"` not `"Kendall's \u03C4"`
- This applies to ALL text rendered by the graphics device: `labs()`, `annotate()`, `geom_text()`, facet labels

### Global Function References

- Every function call must be either defined in the package or namespace-qualified (`pkg::fun()`)
- `setNames()` lives in `stats` -- always use `stats::setNames()`
- `Sys.time()`, `Sys.sleep()` etc. are in `base` and don't need qualification
- Use `R CMD check` NOTE about "no visible global function definition" as the definitive list

### Roxygen2 Cross-References

- `@seealso` links like `[function_name()]` only resolve to functions documented in the same package or its dependencies
- Functions in our own package that aren't yet documented will produce warnings during `devtools::document()` -- these resolve once all docs are in place
- To link to functions in Suggests (not Imports), use the full form: `\code{\link[pkg]{function}}`

### Examples

- Wrap long-running examples in `\donttest{}` (not `\dontrun{}`)
- `\donttest{}` examples ARE run by `--as-cran` but not by default `R CMD check`
- Examples must not require user interaction or internet access
- Examples must not write to the user's file system

### DESCRIPTION

- `Authors@R` field must parse correctly via `eval(parse(text=...))` -- avoid non-ASCII in author names when possible
- `R CMD build` expands `Authors@R` into `Author` and `Maintainer` fields; `R CMD check` on a source directory (not tarball) may fail if these are missing -- always check on the built tarball

## Package Integration (codyna)

### Repository

- codyna lives at `github.com/santikka/codyna`
- Our 6 files go into `codyna/R/` alongside the existing 11 files
- Zero new dependencies needed (all packages already in DESCRIPTION Imports)

### Existing codyna Files

| File | Contents |
|------|----------|
| `check.R` | Input validation: `check_class`, `check_flag`, `check_match`, `check_missing`, `check_range`, `check_string`, `check_values`, `check_formula` |
| `utilities.R` | Internal helpers: `roll`, `rollmean`, `stat_mode`, `entropy`, `rmsqd`, `chisq_test`, `get_cols`, `try_`, `ifelse_`, `onlyif`, `warning_`, `stop_`, `stopifnot_`, `message_`, `n_unique`, `%||%`, `%m%` |
| `data.R` | Data loading and `prepare_timeseries_data()` |
| `codyna-package.R` | Package-level documentation |
| `complexity.R` | `complexity()` export |
| `ews.R` | `detect_warnings()` export |
| `indices.R` | `sequence_indices()` export |
| `patterns.R` | `discover_patterns()` export |
| `plot.R` | Plot methods for existing classes |
| `print.R` | Print methods for existing classes |
| `regimes.R` | `detect_regimes()`, `convert()` exports |

### Our 6 Files

| File | Exports | S3 Methods |
|------|---------|------------|
| `datagen.R` | `generate_ts_data`, `generate_tipping_data`, `compare_ts` | -- |
| `hurst.R` | `hurst`, `detect_hurst_warnings` | plot/print/summary for hurst, hurst_global, hurst_ews |
| `multivariate_ews.R` | `detect_multivariate_warnings` | plot/print/summary for multi_ews |
| `resilience.R` | `resilience`, `classify_resilience` | plot/print for resilience |
| `to_tna.R` | `to_tna` (generic) | to_tna methods for resilience, hurst, hurst_ews, hurst_global, multi_ews, trend, changepoint, spectral, potential, surrogate_test, sensitivity_ews |
| `trend.R` | `compute_trend` | plot/print/summary for trend |

### NAMESPACE Entries (auto-generated by roxygen2)

```
# Exports (10 new)
export(classify_resilience)
export(compare_ts)
export(compute_trend)
export(detect_hurst_warnings)
export(detect_multivariate_warnings)
export(generate_tipping_data)
export(generate_ts_data)
export(hurst)
export(resilience)
export(to_tna)

# S3 methods (17 new)
S3method(plot, hurst)
S3method(plot, hurst_ews)
S3method(plot, multi_ews)
S3method(plot, resilience)
S3method(plot, trend)
S3method(print, hurst)
S3method(print, hurst_ews)
S3method(print, hurst_global)
S3method(print, multi_ews)
S3method(print, resilience)
S3method(print, trend)
S3method(summary, hurst_ews)
S3method(summary, multi_ews)
S3method(summary, trend)
S3method(to_tna, default)
S3method(to_tna, hurst)
S3method(to_tna, hurst_ews)
S3method(to_tna, hurst_global)
S3method(to_tna, multi_ews)
S3method(to_tna, resilience)
S3method(to_tna, trend)
S3method(to_tna, changepoint)
S3method(to_tna, spectral)
S3method(to_tna, potential)
S3method(to_tna, surrogate_test)
S3method(to_tna, sensitivity_ews)
```

### Integration Steps

1. Copy the 6 R files into `codyna/R/`
2. Run `devtools::document()` to regenerate NAMESPACE and man/ pages
3. Run `R CMD build .` then `R CMD check --as-cran {tarball}` on the tarball
4. Verify: Status OK with zero errors, warnings, and notes

### TNA Interoperability

- `to_tna()` returns wide-format tibbles (T1, T2, T3, ...) with one row per sequence
- `tna::tna()` accepts this format directly via its `build_model.data.frame` method
- Do NOT use `tna::prepare_data()` with `to_tna()` output -- `prepare_data` expects long-format event data with actor/time/action columns
- Pattern: `tna::tna(to_tna(obj))` builds a transition network model in one step
