# codyna_experimental

Experimental analysis modules for the [codyna](https://github.com/santikka/codyna) R package. These extensions add resilience quantification, early warning signal detection, and transition network analysis capabilities for complex dynamical systems.

## Modules

| Module | Function(s) | Purpose |
|--------|-------------|---------|
| **Resilience** | `resilience()`, `classify_resilience()` | Rolling-window resilience metrics (VSI, ARCH-LM, CV, recovery time/slope, sample entropy, DFA alpha, AC ratio) with six-level classification |
| **Hurst dynamics** | `hurst()`, `detect_hurst_warnings()` | Rolling Hurst exponent via DFA/RS/Higuchi with memory state classification and multi-indicator early warning scores |
| **Trend** | `compute_trend()` | Local trend classification (ascending, descending, flat, turbulent) using slope, variance ratio, and directional consistency |
| **Multivariate EWS** | `detect_multivariate_warnings()` | 11 multivariate early warning metrics (AR, SD, covariance, MAF, PCA) with expanding and rolling window modes |
| **Changepoint** | `detect_changepoints()` | Structural break detection via PELT, BinSeg, and CUSUM with smart regime clustering, level/direction classification, and return detection |
| **Spectral** | `spectral_ews()` | Rolling spectral exponent tracking noise colour transitions (white, pink, red, brownian) as an independent CSD indicator |
| **Potential** | `potential_analysis()` | Stability landscape reconstruction via kernel density with well/barrier detection and rolling modality tracking |
| **Surrogates** | `surrogate_test()` | Statistical significance testing for EWS trends using phase randomization, AAFT, shuffle, and block bootstrap surrogates |
| **Sensitivity** | `sensitivity_ews()` | Parameter robustness analysis across window sizes and detrending methods |
| **Data generation** | `generate_ts_data()`, `generate_tipping_data()` | Synthetic time series with labelled resilience phases; multivariate tipping-point data with Ornstein-Uhlenbeck or AR critical slowing down |
| **TNA bridge** | `to_tna()` | Extracts categorical state sequences from any analysis object as wide-format tibbles compatible with `tna::tna()` and `TraMineR::seqdef()` |

## S3 Classes

Every module returns a structured S3 object with `plot()`, `print()`, and `summary()` methods:

```r
# Full pipeline example
library(codyna)

x <- generate_ts_data(data_length = 500, generate_plot = FALSE)$value

res <- resilience(x, window = 50)          # Rolling resilience metrics
cls <- classify_resilience(res)            # Six-level classification
h   <- hurst(x, method = "dfa", window = 50, step = 10)  # Hurst dynamics
ews <- detect_hurst_warnings(h)            # Early warning signals
tr  <- compute_trend(x, window = 50)       # Trend states
sp  <- spectral_ews(x, window = 50)        # Spectral noise colour
cp  <- detect_changepoints(x)              # Structural breaks
pa  <- potential_analysis(x, window = 100)  # Stability landscape

plot(cls)    # Resilience ribbon
plot(h)      # Hurst trajectory with states
plot(ews)    # Warning dashboard
plot(cp)     # Changepoint segmentation
plot(sp)     # Spectral exponent with noise states

# Export any state sequence to TNA
tna::tna(to_tna(cls))   # Resilience transitions
tna::tna(to_tna(cp))    # Changepoint regime transitions
tna::tna(to_tna(sp))    # Spectral state transitions
```

## Surrogate Significance Testing

```r
tip <- generate_tipping_data(n_time = 200, n_vars = 3, tipping_point = 100)

st <- surrogate_test(tip$VAR1, metric = "sd", n_surrogates = 199, window = 50)
plot(st)     # Histogram + rolling metric with surrogate envelope
summary(st)  # tau, p-value, significance
```

## Tipping-Point Data Generation

Two simulation models for multivariate critical transitions:

```r
# Ornstein-Uhlenbeck (default): clean critical slowing down
tip_ou <- generate_tipping_data(n_time = 200, n_vars = 5, model = "ou")

# Classic AR: stronger amplitude changes
tip_ar <- generate_tipping_data(n_time = 100, n_vars = 5, tipping_point = 60,
                                 stability_strength = 0.8, model = "ar")
```

## Integration with codyna

These files are designed to drop into `codyna/R/` with zero new dependencies. All imports are already in codyna's DESCRIPTION:

```
checkmate, cli, dplyr, ggplot2, patchwork, rlang, scales, stats, tibble, tidyr, tidyselect
```

Integration steps:

1. Copy `R/*.R` into `codyna/R/`
2. Run `devtools::document()` to regenerate NAMESPACE and man pages
3. Run `R CMD build . && R CMD check --as-cran` on the tarball

## Validation

Two test suites verify numerical correctness:

- **Numerical equivalence** (89 tests): all computations match manual base R implementations
- **External equivalence** (800 iterations): validated against `changepoint`, `EWSmethods`, `nonlinearTseries`, and `stats`

```r
source("tests/test_numerical_equivalence.R")   # 89/89 pass
source("tests/test_external_equivalence.R")     # 800/800 pass
```

## Documentation

- `docs/workflow.Rmd` -- Complete tutorial walking through the full analysis pipeline
- `docs/reference.html` -- Function reference with parameter tables and return descriptions
