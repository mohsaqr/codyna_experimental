# ============================================================================
# Numerical Equivalence Tests
# Verifies that spectral.R, potential.R, surrogates.R, changepoint.R,
# sensitivity.R produce numerically correct results by comparing against
# manual base R computations and known analytical cases.
# ============================================================================

cat("=== Loading codyna namespace and sourcing our files ===\n")
library(codyna)
ns <- asNamespace("codyna")
test_env <- new.env(parent = ns)
our_dir <- file.path(
  "/Users/mohammedsaqr/Documents/Documents X/Resileience_Claude/claude/R"
)
files <- list.files(our_dir, pattern = "[.]R$", full.names = TRUE)
for (f in files) source(f, local = test_env)

PASS <- 0L
FAIL <- 0L
check <- function(label, condition) {
  if (isTRUE(condition)) {
    PASS <<- PASS + 1L
    cat("  PASS:", label, "\n")
  } else {
    FAIL <<- FAIL + 1L
    cat("  FAIL:", label, "\n")
  }
}

tol <- 1e-8

# ============================================================================
# 1. SPECTRAL.R -- Numerical equivalence
# ============================================================================
cat("\n=== 1. SPECTRAL.R ===\n")

set.seed(42)
x_spec <- cumsum(rnorm(200))

sp <- test_env$spectral_ews(x_spec, window = 50L, method = "periodogram",
                            detrend = "none", align = "right", states = TRUE)

# --- 1a. Verify spectral exponent via manual stats::spectrum + log-log OLS ---
# Pick window ending at position 100
w_data <- x_spec[51:100]
spec_manual <- stats::spectrum(w_data, method = "pgram", plot = FALSE,
                               detrend = FALSE, taper = 0.1)
freq_m <- spec_manual$freq
power_m <- spec_manual$spec
valid_m <- freq_m > 0 & power_m > 0
lf <- log(freq_m[valid_m])
lp <- log(power_m[valid_m])
xm <- mean(lf)
ym <- mean(lp)
xxv <- sum((lf - xm)^2)
xyc <- sum((lf - xm) * (lp - ym))
slope_manual <- xyc / xxv
beta_manual <- -slope_manual

check("spectral_exponent matches manual (window 51:100)",
      abs(sp$spectral_exponent[100] - beta_manual) < tol)

# --- 1b. Verify spectral ratio ---
med_freq <- stats::median(freq_m[valid_m])
low_p <- sum(power_m[valid_m][freq_m[valid_m] <= med_freq])
high_p <- sum(power_m[valid_m][freq_m[valid_m] > med_freq])
ratio_manual <- low_p / high_p

check("spectral_ratio matches manual",
      abs(sp$spectral_ratio[100] - ratio_manual) < tol)

# --- 1c. Verify R-squared ---
intercept_m <- ym - slope_manual * xm
pred <- intercept_m + slope_manual * lf
ss_tot <- sum((lp - ym)^2)
ss_res <- sum((lp - pred)^2)
r2_manual <- 1 - ss_res / ss_tot

check("r_squared matches manual",
      abs(sp$r_squared[100] - r2_manual) < tol)

# --- 1d. Verify state classification thresholds ---
check("white_noise state: beta < 0.5",
      all(sp$spectral_exponent[sp$state == "white_noise"] < 0.5, na.rm = TRUE))
check("pink_noise state: 0.5 <= beta < 1.5",
      all(sp$spectral_exponent[sp$state == "pink_noise"] >= 0.5 &
          sp$spectral_exponent[sp$state == "pink_noise"] < 1.5, na.rm = TRUE))

# --- 1e. Verify AR method also works ---
sp_ar <- test_env$spectral_ews(x_spec, window = 50L, method = "ar",
                                detrend = "linear", states = TRUE)
check("AR method produces valid output", nrow(sp_ar) == 200L)
check("AR method beta is numeric", all(!is.na(sp_ar$spectral_exponent)))

# --- 1f. Verify linear detrending ---
w_data_dt <- x_spec[51:100]
n_dt <- length(w_data_dt)
t_idx <- seq_len(n_dt)
t_m <- mean(t_idx)
x_m <- mean(w_data_dt)
xx_v <- sum((t_idx - t_m)^2)
xy_c <- sum((t_idx - t_m) * (w_data_dt - x_m))
sl <- xy_c / xx_v
ic <- x_m - sl * t_m
detrended_manual <- w_data_dt - (ic + sl * t_idx)
detrended_fn <- test_env$spectral_detrend_(w_data_dt, "linear")
check("linear detrend matches manual OLS",
      max(abs(detrended_manual - detrended_fn)) < tol)

# --- 1g. Verify diff detrend ---
diff_manual <- diff(w_data_dt)
diff_fn <- test_env$spectral_detrend_(w_data_dt, "diff")
check("diff detrend matches base diff()",
      max(abs(diff_manual - diff_fn)) < tol)

# ============================================================================
# 2. POTENTIAL.R -- Numerical equivalence
# ============================================================================
cat("\n=== 2. POTENTIAL.R ===\n")

# --- 2a. Known case: single well (Ornstein-Uhlenbeck) ---
set.seed(42)
n_ou <- 5000
ou <- numeric(n_ou)
for (i in 2:n_ou) ou[i] <- ou[i - 1] - 0.5 * ou[i - 1] * 0.01 + 0.3 * rnorm(1)

pa <- test_env$potential_analysis(ou, n_bins = 50L)
check("OU process: exactly 1 well detected", pa$n_wells == 1L)
# Well should be near the series mean (mean-reverting process)
check("OU well near series mean (within 1.0)",
      abs(pa$wells$location[1] - mean(ou)) < 1.0)

# --- 2b. Known case: bimodal distribution ---
set.seed(42)
bimodal <- c(rnorm(2000, -3, 0.5), rnorm(2000, 3, 0.5))
pa2 <- test_env$potential_analysis(bimodal, n_bins = 50L)
check("Bimodal: 2 wells detected", pa2$n_wells == 2L)
check("Bimodal: well 1 near -3", abs(pa2$wells$location[1] - (-3)) < 0.5)
check("Bimodal: well 2 near +3", abs(pa2$wells$location[2] - 3) < 0.5)
check("Bimodal: at least 1 barrier", nrow(pa2$barriers) >= 1L)

# --- 2c. Verify U(x) = -log(density) manually ---
dens_manual <- stats::density(bimodal, n = 200)
eps <- .Machine$double.eps
pot_manual <- -log(pmax(dens_manual$y, eps))
# Our function uses n_bins * 4 grid
dens_ours <- stats::density(bimodal, n = 200)
check("potential = -log(density) at grid point 50",
      abs(pa2$landscape$potential[50] -
          (-log(pmax(pa2$landscape$density[50], eps)))) < tol)

# --- 2d. Verify rolling window mode ---
set.seed(42)
x_roll <- cumsum(rnorm(200))
par <- test_env$potential_analysis(x_roll, window = 80L, n_bins = 50L)
check("Rolling: output has 'rolling' element", !is.null(par$rolling))
check("Rolling: rolling tibble has correct length", nrow(par$rolling) == 200L)
check("Rolling: first (window-1) entries are NA",
      all(is.na(par$rolling$n_wells[1:79])))
check("Rolling: position 80+ has values",
      !is.na(par$rolling$n_wells[80]))

# --- 2e. Verify linear detrending ---
set.seed(42)
x_trend <- 1:100 + rnorm(100)
fit <- stats::lm.fit(x = cbind(1, 1:100), y = x_trend)
resid_manual <- fit$residuals
resid_fn <- test_env$potential_detrend_(x_trend, "linear")
check("potential linear detrend matches lm.fit residuals",
      max(abs(resid_manual - resid_fn)) < tol)

# ============================================================================
# 3. SURROGATES.R -- Statistical property verification
# ============================================================================
cat("\n=== 3. SURROGATES.R ===\n")

set.seed(42)
x_surr <- cumsum(rnorm(200))

# --- 3a. FFT phase randomization preserves power spectrum ---
set.seed(42)
surr_phase <- test_env$surrogate_phase_(x_surr)
check("Phase surrogate: same length", length(surr_phase) == length(x_surr))
check("Phase surrogate: output is real", all(!is.complex(surr_phase)))

# Raw FFT magnitude comparison (phase randomization preserves |FFT|^2 exactly,
# but stats::spectrum() applies smoothing/tapering that creates differences)
orig_fft_mag <- Mod(stats::fft(x_surr))^2
surr_fft_mag <- Mod(stats::fft(surr_phase))^2
check("Phase surrogate: raw |FFT|^2 preserved (mean difference < 1e-6)",
      abs(mean(orig_fft_mag) - mean(surr_fft_mag)) < 1e-6)

# --- 3b. Shuffle surrogate preserves amplitude distribution ---
set.seed(42)
surr_shuf <- test_env$surrogate_shuffle_(x_surr)
check("Shuffle: same length", length(surr_shuf) == length(x_surr))
check("Shuffle: same sorted values (identical distribution)",
      all(sort(x_surr) == sort(surr_shuf)))

# --- 3c. AAFT preserves amplitude distribution ---
set.seed(42)
surr_aaft <- test_env$surrogate_aaft_(x_surr)
check("AAFT: same length", length(surr_aaft) == length(x_surr))
check("AAFT: same sorted values (preserves distribution)",
      max(abs(sort(x_surr) - sort(surr_aaft))) < tol)

# --- 3d. Block bootstrap: correct length and real values ---
set.seed(42)
surr_block <- test_env$surrogate_block_(x_surr)
check("Block: same length", length(surr_block) == length(x_surr))
check("Block: all values are from original",
      all(surr_block %in% x_surr))

# --- 3e. Verify AR1 metric computation against acf() ---
set.seed(42)
seg <- rnorm(50)
ar1_ours <- test_env$surrogate_compute_metric_(seg, "ar1")
ar1_manual <- stats::acf(seg, lag.max = 1L, plot = FALSE)$acf[2L]
check("AR1 metric matches stats::acf()",
      abs(ar1_ours - ar1_manual) < tol)

# --- 3f. Verify SD metric ---
sd_ours <- test_env$surrogate_compute_metric_(seg, "sd")
sd_manual <- stats::sd(seg)
check("SD metric matches stats::sd()",
      abs(sd_ours - sd_manual) < tol)

# --- 3g. Verify skewness computation ---
skew_ours <- test_env$surrogate_compute_metric_(seg, "skewness")
mu <- mean(seg); s <- stats::sd(seg)
skew_manual <- mean(((seg - mu) / s)^3)
check("Skewness metric matches manual formula",
      abs(skew_ours - skew_manual) < tol)

# --- 3h. Verify kurtosis computation ---
kurt_ours <- test_env$surrogate_compute_metric_(seg, "kurtosis")
kurt_manual <- mean(((seg - mu) / s)^4) - 3
check("Kurtosis metric matches manual formula",
      abs(kurt_ours - kurt_manual) < tol)

# --- 3i. Verify Kendall tau computation ---
set.seed(42)
metric_vals <- 1:20 + rnorm(20, sd = 0.1)  # Strong upward trend
tau_ours <- test_env$surrogate_kendall_tau_(metric_vals)
tau_manual <- stats::cor.test(1:20, metric_vals, method = "kendall")$estimate
check("Kendall tau matches cor.test()",
      abs(tau_ours - tau_manual) < tol)

# --- 3j. Full surrogate test: trending series should be significant ---
set.seed(42)
x_trending <- cumsum(rnorm(200)) + seq(0, 10, length.out = 200)
st <- test_env$surrogate_test(x_trending, n_surrogates = 99L,
                               metric = "ar1", window = 30L)
check("Surrogate test: p_value is numeric", is.numeric(st$p_value))
check("Surrogate test: p_value in [0,1]",
      st$p_value >= 0 && st$p_value <= 1)
check("Surrogate test: observed_tau is numeric",
      is.numeric(st$observed_tau) && !is.na(st$observed_tau))
check("Surrogate test: correct number of surrogates",
      length(st$surrogate_taus) == 99L)

# --- 3k. DFA Hurst on known series ---
# White noise: H should be ~0.5
set.seed(42)
wn <- rnorm(500)
h_wn <- test_env$surrogate_hurst_dfa_(wn)
check("DFA Hurst on white noise: H near 0.5",
      abs(h_wn - 0.5) < 0.15)

# Random walk: H should be ~1.0 (or at least > 0.7)
set.seed(42)
rw <- cumsum(rnorm(500))
h_rw <- test_env$surrogate_hurst_dfa_(rw)
check("DFA Hurst on random walk: H > 0.7", h_rw > 0.7)

# --- 3l. Verify spectral exponent in surrogates ---
set.seed(42)
seg100 <- rnorm(100)
se_ours <- test_env$surrogate_compute_metric_(seg100, "spectral_exponent")
sp_ref <- stats::spectrum(seg100, method = "pgram", plot = FALSE)
valid_ref <- sp_ref$freq > 0 & sp_ref$spec > 0
lf_ref <- log(sp_ref$freq[valid_ref])
ls_ref <- log(sp_ref$spec[valid_ref])
fit_ref <- stats::lm.fit(x = cbind(1, lf_ref), y = ls_ref)
se_manual <- fit_ref$coefficients[2L]
check("Spectral exponent metric matches manual (note: raw slope, not negated)",
      abs(se_ours - se_manual) < tol)

# ============================================================================
# 4. CHANGEPOINT.R -- Numerical equivalence
# ============================================================================
cat("\n=== 4. CHANGEPOINT.R ===\n")

# --- 4a. Known mean shift: should detect changepoint near 200 ---
set.seed(42)
x_shift <- c(rnorm(200, 0, 1), rnorm(200, 10, 1))
cp <- test_env$detect_changepoints(x_shift, method = "pelt",
                                    penalty = "bic", type = "mean")
cp_locs <- attr(cp, "changepoint_locations")
check("PELT detects at least 1 changepoint", length(cp_locs) >= 1L)
check("PELT changepoint near 200 (within 5)",
      any(abs(cp_locs - 200) <= 5))

# --- 4b. Verify ALL segment means match manual calculation ---
# PELT may detect more than 1 changepoint, so verify every segment
seg_ids <- sort(unique(cp$segment))
all_means_ok <- TRUE
all_vars_ok <- TRUE
for (sid in seg_ids) {
  idx <- which(cp$segment == sid)
  seg_vals <- x_shift[idx]
  manual_mean <- mean(seg_vals, na.rm = TRUE)
  manual_var <- if (length(seg_vals) > 1) stats::var(seg_vals) else 0
  if (abs(cp$segment_mean[idx[1]] - manual_mean) >= tol) all_means_ok <- FALSE
  if (abs(cp$segment_var[idx[1]] - manual_var) >= tol) all_vars_ok <- FALSE
}
check("All segment means match manual calculation", all_means_ok)
check("All segment variances match manual calculation", all_vars_ok)

# --- 4d. Verify cost function for each type ---
cost_fn <- test_env$changepoint_cost_
vals_seg <- rnorm(50)
n_seg <- length(vals_seg)
mu_seg <- mean(vals_seg)

# type="mean": RSS (sum of squared residuals)
cost_rss <- sum((vals_seg - mu_seg)^2)
cost_ours_mean <- cost_fn(vals_seg, "mean")
check("Cost function (mean) matches RSS",
      abs(cost_ours_mean - cost_rss) < tol)

# type="variance": n*log(var_MLE)
v_seg <- sum((vals_seg - mu_seg)^2) / n_seg
cost_nlogv <- n_seg * log(v_seg)
cost_ours_var <- cost_fn(vals_seg, "variance")
check("Cost function (variance) matches n*log(var)",
      abs(cost_ours_var - cost_nlogv) < tol)

# type="both": n*log(var_MLE) + n
cost_both <- n_seg * log(v_seg) + n_seg
cost_ours_both <- cost_fn(vals_seg, "both")
check("Cost function (both) matches n*log(var)+n",
      abs(cost_ours_both - cost_both) < tol)

# --- 4e. Verify BIC penalty = 2*log(n) for mean, 3*log(n) for both ---
pen_bic_mean <- test_env$changepoint_penalty_("bic", NULL, 400L, "mean")
check("BIC penalty (mean) = 2*log(n)", abs(pen_bic_mean - 2 * log(400)) < tol)

pen_bic_both <- test_env$changepoint_penalty_("bic", NULL, 400L, "both")
check("BIC penalty (both) = 3*log(n)", abs(pen_bic_both - 3 * log(400)) < tol)

# --- 4f. Verify AIC penalty = 4 for mean, 6 for both ---
pen_aic_mean <- test_env$changepoint_penalty_("aic", NULL, 400L, "mean")
check("AIC penalty (mean) = 4", abs(pen_aic_mean - 4) < tol)

pen_aic_both <- test_env$changepoint_penalty_("aic", NULL, 400L, "both")
check("AIC penalty (both) = 6", abs(pen_aic_both - 6) < tol)

# --- 4g. Binary segmentation on clear signal ---
cp_bs <- test_env$detect_changepoints(x_shift, method = "binary_segmentation",
                                       penalty = "bic", type = "mean")
cp_locs_bs <- attr(cp_bs, "changepoint_locations")
check("BinSeg detects changepoint near 200",
      any(abs(cp_locs_bs - 200) <= 5))

# --- 4h. CUSUM on clear signal ---
cp_cu <- test_env$detect_changepoints(x_shift, method = "cusum",
                                       penalty = "bic", type = "mean")
cp_locs_cu <- attr(cp_cu, "changepoint_locations")
check("CUSUM detects changepoint near 200",
      any(abs(cp_locs_cu - 200) <= 10))

# --- 4i. Multiple changepoints ---
set.seed(123)
x_multi <- c(rnorm(100, 0), rnorm(100, 5), rnorm(100, -2), rnorm(100, 8))
cp_multi <- test_env$detect_changepoints(x_multi, method = "pelt",
                                          min_segment = 10L)
n_cp <- attr(cp_multi, "n_changepoints")
check("Multiple shifts: at least 3 changepoints", n_cp >= 3L)

# --- 4j. Segment labels are consistent ---
check("Segment IDs are sequential from 1",
      identical(sort(unique(cp_multi$segment)),
                seq_len(max(cp_multi$segment))))

# --- 4k. Changepoint indicator is at right positions ---
cp_indicator_idx <- which(cp_multi$changepoint)
cp_locs_multi <- attr(cp_multi, "changepoint_locations")
expected_indicator <- cp_locs_multi + 1L
check("Changepoint indicators at cps + 1",
      all(cp_indicator_idx %in% expected_indicator))

# --- 4l. CUSUM statistic: verify against manual cumsum ---
vals_test <- c(1, 2, 3, 10, 11, 12)
mu_test <- mean(vals_test)
cusum_manual <- cumsum(vals_test - mu_test)
max_idx <- which.max(abs(cusum_manual))
cusum_fn <- test_env$changepoint_cusum_stat_(vals_test)
check("CUSUM stat index matches manual which.max(abs(cumsum))",
      cusum_fn == max_idx)

# --- 4m. No changepoints in constant series ---
x_const <- rep(5, 100)
cp_const <- test_env$detect_changepoints(x_const, method = "pelt")
check("Constant series: 0 changepoints",
      attr(cp_const, "n_changepoints") == 0L)

# ============================================================================
# 5. SENSITIVITY.R -- Numerical equivalence
# ============================================================================
cat("\n=== 5. SENSITIVITY.R ===\n")

set.seed(42)
x_sens <- cumsum(rnorm(200, sd = seq(0.5, 2, length.out = 200)))

# --- 5a. Verify AR1 metric against stats::ar.ols ---
w <- 50L
seg_ar <- x_sens[1:w]
ar1_sens <- test_env$sensitivity_metric_(seg_ar, "ar1")
ar_model <- stats::ar.ols(seg_ar, aic = FALSE, order.max = 1L,
                           demean = FALSE, intercept = TRUE)
ar1_manual <- max(-0.999, min(0.999, ar_model$ar[1L]))
check("sensitivity AR1 matches ar.ols",
      abs(ar1_sens - ar1_manual) < tol)

# --- 5b. Verify SD metric ---
sd_sens <- test_env$sensitivity_metric_(seg_ar, "sd")
check("sensitivity SD matches stats::sd",
      abs(sd_sens - stats::sd(seg_ar)) < tol)

# --- 5c. Verify variance metric ---
var_sens <- test_env$sensitivity_metric_(seg_ar, "variance")
check("sensitivity variance matches stats::var",
      abs(var_sens - stats::var(seg_ar)) < tol)

# --- 5d. Verify skewness ---
skew_sens <- test_env$sensitivity_metric_(seg_ar, "skewness")
mu_s <- mean(seg_ar); sd_s <- stats::sd(seg_ar)
skew_s <- mean(((seg_ar - mu_s) / sd_s)^3)
check("sensitivity skewness matches manual",
      abs(skew_sens - skew_s) < tol)

# --- 5e. Verify kurtosis ---
kurt_sens <- test_env$sensitivity_metric_(seg_ar, "kurtosis")
kurt_s <- mean(((seg_ar - mu_s) / sd_s)^4) - 3
check("sensitivity kurtosis matches manual",
      abs(kurt_sens - kurt_s) < tol)

# --- 5f. Verify CV ---
cv_sens <- test_env$sensitivity_metric_(seg_ar, "cv")
cv_manual <- stats::sd(seg_ar) / abs(mean(seg_ar))
check("sensitivity CV matches SD/|mean|",
      abs(cv_sens - cv_manual) < tol)

# --- 5g. Verify Kendall tau against cor.test ---
rolling_scores <- vapply(seq_len(200 - w + 1), function(i) {
  test_env$sensitivity_metric_(x_sens[i:(i + w - 1L)], "ar1")
}, numeric(1L))
tau_sens <- test_env$sensitivity_tau_(rolling_scores)
valid_s <- !is.na(rolling_scores)
tau_ref <- stats::cor.test(which(valid_s), rolling_scores[valid_s],
                            method = "kendall")$estimate
check("sensitivity Kendall tau matches cor.test",
      abs(tau_sens - tau_ref) < tol)

# --- 5h. Full sensitivity_ews output structure ---
sa <- test_env$sensitivity_ews(x_sens, metric = "ar1",
                                windows = c(30L, 50L, 80L),
                                detrend_methods = c("none", "linear"))
check("Output has 5 columns",
      all(c("window", "detrend", "time", "score", "tau") %in% names(sa)))
check("All 6 combinations present",
      nrow(unique(sa[, c("window", "detrend")])) == 6L)

# --- 5i. Verify one specific (window, detrend) combination manually ---
# window=50, detrend=none
sa_sub <- sa[sa$window == 50L & sa$detrend == "none", ]
# Manual computation
rolling_manual <- vapply(seq_len(200 - 50 + 1), function(i) {
  test_env$sensitivity_metric_(x_sens[i:(i + 49)], "ar1")
}, numeric(1L))
check("Rolling AR1 values match manual (window=50, detrend=none)",
      max(abs(sa_sub$score - rolling_manual), na.rm = TRUE) < tol)

# --- 5j. Verify linear detrending matches lm residuals ---
detr <- test_env$sensitivity_detrend_(x_sens, "linear")
fit_s <- stats::lm(x_sens ~ seq_along(x_sens))
resid_s <- stats::residuals(fit_s)
check("sensitivity linear detrend matches lm() residuals",
      max(abs(detr - resid_s)) < tol)

# --- 5k. Tau is consistent across all rows in a (window, detrend) group ---
for (w_val in c(30, 50, 80)) {
  for (d_val in c("none", "linear")) {
    sub <- sa[sa$window == w_val & sa$detrend == d_val, ]
    tau_unique <- unique(sub$tau)
    check(sprintf("Tau constant within (w=%d, d=%s)", w_val, d_val),
          length(tau_unique) == 1L)
  }
}

# ============================================================================
# 6. MULTIVARIATE EWS -- Key metrics verification
# ============================================================================
cat("\n=== 6. MULTIVARIATE_EWS.R (key metrics) ===\n")

set.seed(42)
tp <- test_env$generate_tipping_data(n_time = 80L, n_vars = 3L, tipping_point = 50L, saturation_point = 70L)

# --- 6a. Rolling method produces correct output ---
ews_roll <- test_env$detect_multivariate_warnings(tp, method = "rolling",
                                                   window = 40)
check("Rolling: class is multi_ews", inherits(ews_roll, "multi_ews"))
check("Rolling: method attr is 'rolling'", attr(ews_roll, "method") == "rolling")
check("Rolling: has metric column", "metric" %in% names(ews_roll))

# --- 6b. Expanding method produces correct output ---
ews_exp <- test_env$detect_multivariate_warnings(tp, method = "expanding",
                                                  window = 40)
check("Expanding: class is multi_ews", inherits(ews_exp, "multi_ews"))
check("Expanding: method attr is 'expanding'",
      attr(ews_exp, "method") == "expanding")
check("Expanding: has classification",
      !is.null(attr(ews_exp, "classification")))

# --- 6c. Verify meanSD for a known window ---
ts_data <- as.matrix(tp[, -1])
win_pct <- 40
win_pts <- round(nrow(ts_data) * win_pct / 100)
# First window: rows 1:win_pts
w_data <- ts_data[1:win_pts, ]
sds_manual <- apply(w_data, 2, stats::sd)
meanSD_manual <- mean(sds_manual)
maxSD_manual <- max(sds_manual)

# The rolling output should contain these values for the first window
roll_meanSD <- ews_roll[ews_roll$metric == "meanSD", ]
check("meanSD first value matches manual",
      abs(roll_meanSD$score[1] - meanSD_manual) < 1e-6)

roll_maxSD <- ews_roll[ews_roll$metric == "maxSD", ]
check("maxSD first value matches manual",
      abs(roll_maxSD$score[1] - maxSD_manual) < 1e-6)

# --- 6d. Verify AR1 computation (manual ar.ols) ---
ar_manual <- suppressWarnings(vapply(seq_len(ncol(w_data)), function(j) {
  col_data <- w_data[, j]
  if (stats::var(col_data) < 1e-10) return(NA_real_)
  m <- try(stats::ar.ols(col_data, order.max = 1, dmean = FALSE,
                          intercept = FALSE), silent = TRUE)
  if (inherits(m, "try-error") || length(m$ar) == 0) return(NA_real_)
  m$ar[1]
}, numeric(1)))
meanAR_manual <- mean(ar_manual, na.rm = TRUE)
maxAR_manual <- max(ar_manual, na.rm = TRUE)

roll_meanAR <- ews_roll[ews_roll$metric == "meanAR", ]
roll_maxAR <- ews_roll[ews_roll$metric == "maxAR", ]
check("meanAR first value matches manual ar.ols",
      abs(roll_meanAR$score[1] - meanAR_manual) < 1e-6)
check("maxAR first value matches manual ar.ols",
      abs(roll_maxAR$score[1] - maxAR_manual) < 1e-6)

# --- 6e. Verify eigenCOV: dominant eigenvalue of covariance matrix ---
w_scaled <- scale(w_data)
V_scaled <- stats::cov(w_scaled)
eig_manual <- eigen(V_scaled)$values[1]

roll_eigenCOV <- ews_roll[ews_roll$metric == "eigenCOV", ]
check("eigenCOV first value matches manual eigen(cov(scale(data)))",
      abs(roll_eigenCOV$score[1] - eig_manual) < 1e-6)

# --- 6f. Verify maxCOV: max off-diagonal of covariance ---
maxCOV_manual <- max(V_scaled[lower.tri(V_scaled, diag = FALSE)])

roll_maxCOV <- ews_roll[ews_roll$metric == "maxCOV", ]
check("maxCOV first value matches manual",
      abs(roll_maxCOV$score[1] - maxCOV_manual) < 1e-6)

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n", strrep("=", 50), "\n")
cat(sprintf("TOTAL: %d PASS, %d FAIL out of %d tests\n",
            PASS, FAIL, PASS + FAIL))
if (FAIL == 0L) {
  cat("ALL TESTS PASSED -- numerical equivalence verified!\n")
} else {
  cat("SOME TESTS FAILED -- review output above.\n")
}
cat(strrep("=", 50), "\n")
