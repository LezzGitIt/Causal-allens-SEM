# Key_causal_fns.R
# Core functions for the causal Allen's rule simulation project.
# Source this file from any script that needs these functions.
#
# Required packages: tidyverse, smatr, lavaan, sliR

# ── Algebraic solvers for simulation parameter grid ──────────────────────────
#
# All three functions account for the full causal structure.
# sigma_Z^2 = beta_size^2 * sigma_temp^2 + sigma_size^2  (total latent size variance)
# Cov(log_Wing, log_Mass) has two paths:
#   (1) size residual:  lambda_w * lambda_m * sigma_size^2
#   (2) temperature:    B_temp * lambda_m * beta_size * sigma_temp^2
# Ignoring these terms causes r(Wing, Mass) ≠ target_r in simulation.

solve_lambda_w <- function(
    target_r,
    target_sma,
    lambda_m   = 1,
    sigma_size = 0.10,
    sigma_mass = 0.10,
    beta_size  = 0,
    B_temp     = 0,
    sigma_temp = 0.3
){
  sigma_z_sq <- beta_size^2 * sigma_temp^2 + sigma_size^2
  var_m      <- lambda_m^2 * sigma_z_sq + sigma_mass^2
  (target_r * target_sma * var_m - B_temp * lambda_m * beta_size * sigma_temp^2) /
    (lambda_m * sigma_z_sq)
}

solve_sigma_w <- function(
    target_sma,
    lambda_w,
    lambda_m   = 1,
    sigma_size = 0.10,
    sigma_mass = 0.10,
    beta_size  = 0,
    B_temp     = 0,
    sigma_temp = 0.3
){
  sigma_z_sq <- beta_size^2 * sigma_temp^2 + sigma_size^2
  var_m      <- lambda_m^2 * sigma_z_sq + sigma_mass^2
  sigma_w_sq <- target_sma^2 * var_m -
    (B_temp + lambda_w * beta_size)^2 * sigma_temp^2 -
    lambda_w^2 * sigma_size^2
  dplyr::if_else(sigma_w_sq > 0, sqrt(sigma_w_sq), NA_real_)
}

# Returns the implied r(Wing,Mass), OLS slope, and SMA slope for a given
# (lambda_w, sigma_w) combination; useful for verifying solver output.
implied_allometry <- function(
    lambda_w,
    sigma_w,
    lambda_m   = 1,
    sigma_size = 0.10,
    sigma_mass = 0.10,
    beta_size  = 0,
    B_temp     = 0,
    sigma_temp = 0.3
){
  sigma_z_sq <- beta_size^2 * sigma_temp^2 + sigma_size^2
  var_m      <- lambda_m^2 * sigma_z_sq + sigma_mass^2
  var_w      <- (B_temp + lambda_w * beta_size)^2 * sigma_temp^2 +
                lambda_w^2 * sigma_size^2 + sigma_w^2
  cov_mw     <- lambda_m * (B_temp * beta_size * sigma_temp^2 + lambda_w * sigma_z_sq)
  tibble(
    rho      = cov_mw / sqrt(var_m * var_w),
    beta_ols = cov_mw / var_m,
    beta_sma = sqrt(var_w / var_m)
  )
}

# ── Utility helpers ──────────────────────────────────────────────────────────

format_temp <- function(df) {
  df %>% mutate(
    Temp_bin = cut(Temp_inc, breaks = 15, labels = FALSE, ordered_result = TRUE)
  ) %>%
    arrange(Temp_inc) %>%
    mutate(Temp_bin = case_when(
      Temp_bin %in% c(1:5)   ~ 5,
      Temp_bin %in% c(11:15) ~ 11,
      .default = Temp_bin
    ))
}

# ── Data generation ───────────────────────────────────────────────────────────
#
# Causal structure:
#   Temp → latent Size (Bergmann's rule: beta_size < 0)
#   Size → Mass, App1 (Wing), App2 (Tarsus)
#   Temp → App1 (direct Allen's, B_temp_w)
#   Temp → App2 (direct Allen's, B_temp_t; in parms_mat set equal to B_temp_w)
#   Mass measured with error sigma_me: M1, M2, M3 ~ N(Mass, sigma_me)
#   sigma_me = 0 → M1 = M2 = M3 = Mass (no error)

gen_causal_data <- function(
    N             = 1500,
    B_temp_w      = 0.3,
    B_temp_t      = B_temp_w / 2,
    beta_size     = -0.4,
    sigma_size    = 0.10,     # SD of latent log size; ~10% CV on raw scale
    lambda_mass   = 1.0,      # mass loading on log scale
    lambda_wing   = 0.27,     # wing loading on log scale (solve from targets)
    lambda_tarsus = 0.27,     # tarsus loading on log scale
    mu_mass       = log(50),  # mean log mass
    mu_wing       = log(180), # mean log wing (mm)
    mu_tarsus     = log(40),  # mean log tarsus (mm)
    sigma_mass    = 0.10,     # SD of log mass residual
    sigma_wing    = 0.0577,   # SD of log wing residual (solve from targets)
    sigma_tarsus  = 0.0577,   # SD of log tarsus residual
    sigma_temp    = 0.3,
    sigma_me      = 0         # measurement error on log scale (~CV for small values)
) {

  Temp     <- rnorm(N, 1, sigma_temp)
  log_Z    <- rnorm(N, beta_size * (Temp - 1), sigma_size)  # centred: E[log_Z]=0 at mean Temp
  log_Mass <- mu_mass   + lambda_mass   * log_Z + rnorm(N, 0, sigma_mass)
  log_Wing <- mu_wing   + B_temp_w * (Temp - 1) + lambda_wing   * log_Z + rnorm(N, 0, sigma_wing)
  log_Tars <- mu_tarsus + B_temp_t * (Temp - 1) + lambda_tarsus * log_Z + rnorm(N, 0, sigma_tarsus)

  M1_log <- rnorm(N, log_Mass, sigma_me)
  M2_log <- rnorm(N, log_Mass, sigma_me)
  M3_log <- rnorm(N, log_Mass, sigma_me)

  df <- tibble(
    Temp_inc   = Temp,
    Size       = log_Z,
    Mass       = exp(log_Mass),
    Append     = exp(log_Wing),
    Tarsus     = exp(log_Tars),
    Mass_log   = log_Mass,
    Append_log = log_Wing,
    Tarsus_log = log_Tars,
    M1 = exp(M1_log), M2 = exp(M2_log), M3 = exp(M3_log),
    M1_log = M1_log, M2_log = M2_log, M3_log = M3_log
  )
  # Trim log-scale outliers at ±2.5 SD per variable (~1.2% per tail removed).
  # Log-scale trimming is symmetric and preserves correlation structure;
  # raw-scale trimming is asymmetric on a lognormal and leaves the right tail intact.
  df %>%
    filter(
      abs(Mass_log   - mean(Mass_log))   < 2.5 * sd(Mass_log),
      abs(Append_log - mean(Append_log)) < 2.5 * sd(Append_log),
      abs(Tarsus_log - mean(Tarsus_log)) < 2.5 * sd(Tarsus_log)
    )
}

# Wrapper that also returns the 4×4 correlation matrix for the realism filter.
simulate_data <- function(..., N = 1500) {
  df   <- gen_causal_data(..., N = N)
  corr <- cor(df[c("Temp_inc", "Mass", "Append", "Tarsus")])
  list(data = df, corr = corr)
}

# ── Classical comparison methods ──────────────────────────────────────────────

# Ratio, Ryding (mass-as-covariate), and two SLI variants applied to one appendage.
# Coefficients are on the z-scored scale; backtransform with sd_append / sigma_temp.
#
# SLI (Standardized Length Index; Peig & Green 2009), via sliR::calc_sli() on the raw
# scale: SLI_i = App_i * (M0/Mass_i)^b. We log it to stay comparable to the other
# logged metrics: log(SLI) = App_log - b*Mass_log + b*log(M0), a size-corrected log
# length. The b*log(M0) term is an additive constant (M0 = mean(Mass) by default), which
# drops out under z-scoring, so it does not affect the fitted Temp_inc slope below.
#   - SLI (isometry): b = 0.33 (geometric similarity, length ∝ mass^(1/3))
#   - SLI (estimated): b = SMA slope of log(appendage) ~ log(mass)
apply_methods <- function(sim_df, append_log = "Append_log", append_raw = "Append") {
  sim_df2 <- sim_df %>%
    rename(App_log = all_of(append_log), App = all_of(append_raw)) %>%
    mutate(app_mass = App_log / Mass_log)

  b_sma <- coef(smatr::sma(App_log ~ Mass_log, data = sim_df2))["slope"]
  sim_df2 <- sim_df2 %>%
    sliR::calc_sli(Append = App, Mass = Mass, b_sli = 0.33,  rename_col = "sli_iso") %>%
    sliR::calc_sli(Append = App, Mass = Mass, b_sli = b_sma, rename_col = "sli_est") %>%
    mutate(sli_iso = log(sli_iso),
           sli_est = log(sli_est))

  sim_df_s <- sim_df2 %>% mutate(across(where(is.numeric), \(x) as.numeric(scale(x))))

  ryding_m  <- lm(App_log  ~ Mass_log + Temp_inc, data = sim_df_s)
  ratio_m   <- lm(app_mass ~ Temp_inc,            data = sim_df_s)
  sli_iso_m <- lm(sli_iso  ~ Temp_inc,            data = sim_df_s)
  sli_est_m <- lm(sli_est  ~ Temp_inc,            data = sim_df_s)

  tibble(
    coef_ryding  = coef(ryding_m)["Temp_inc"],
    coef_ratio   = coef(ratio_m)["Temp_inc"],
    coef_sli_iso = coef(sli_iso_m)["Temp_inc"],
    coef_sli_est = coef(sli_est_m)["Temp_inc"]
  )
}

# ── Latent SEM (lavaan) ───────────────────────────────────────────────────────

# General MIMIC model: Size =~ 1*Mass + App1 + ...; Size ~ Temp; App1 ~ Temp; ...
# Data must be pre-scaled (z-scored) by the caller before passing.
# Returns a one-row tibble: coef_sem_{label} and se_sem_{label} for each appendage.
fit_lavaan_sem <- function(df,
                           mass_name    = "Mass_log",
                           append_names = c("Append_log", "Tarsus_log"),
                           temp_name    = "Temp_inc",
                           labels       = c("wing", "tarsus"),
                           size_covs    = character(0)) {
  app_ids <- paste0("App", seq_along(append_names))
  col_map <- c(setNames(mass_name, "Mass"),
               setNames(append_names, app_ids),
               setNames(temp_name, "Temp"))

  df_s <- df %>%
    dplyr::select(all_of(c(mass_name, append_names, temp_name, size_covs))) %>%
    rename(all_of(col_map)) %>%
    mutate(across(where(is.numeric), as.numeric))

  # lavaan requires numeric predictors: encode character/factor size_covs as 0/1
  # (reference level = first alphabetical value, same convention as lm)
  for (cov in size_covs) {
    if (!is.numeric(df_s[[cov]])) {
      ref <- sort(unique(na.omit(df_s[[cov]])))[1]
      df_s[[cov]] <- as.numeric(df_s[[cov]] != ref)
    }
  }

  size_rhs <- paste(c("Temp", size_covs), collapse = " + ")
  model <- paste0(
    "Size =~ 1*Mass + ", paste(app_ids, collapse = " + "), "\n",
    "Size ~ ", size_rhs, "\n",
    paste(paste0(app_ids, " ~ Temp"), collapse = "\n")
  )

  coef_names      <- paste0("coef_sem_",      labels)
  se_names        <- paste0("se_sem_",        labels)
  lambda_names    <- paste0("lambda_sem_",    labels)
  se_lambda_names <- paste0("se_lambda_sem_", labels)
  na_ret <- as_tibble(setNames(
    as.list(rep(NA_real_, 4 * length(labels))),
    c(coef_names, se_names, lambda_names, se_lambda_names)
  ))

  fit <- tryCatch(
    lavaan::sem(model, data = df_s, estimator = "ML"),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(lavaan::sem(model, data = df_s, estimator = "ML"))
  )
  if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) return(na_ret)

  all_ests <- lavaan::parameterEstimates(fit)

  ests  <- all_ests %>% filter(op == "~",  rhs == "Temp")
  coefs <- map_dbl(app_ids, \(a) ests %>% filter(lhs == a) %>% pull(est))
  ses   <- map_dbl(app_ids, \(a) ests %>% filter(lhs == a) %>% pull(se))

  # Guard against degenerate solutions
  if (any(abs(coefs) > 1.5)) return(na_ret)

  load_ests  <- all_ests %>% filter(op == "=~", lhs == "Size", rhs != "Mass")
  lambdas    <- map_dbl(app_ids, \(a) load_ests %>% filter(rhs == a) %>% pull(est))
  se_lambdas <- map_dbl(app_ids, \(a) load_ests %>% filter(rhs == a) %>% pull(se))

  as_tibble(setNames(
    as.list(c(coefs, ses, lambdas, se_lambdas)),
    c(coef_names, se_names, lambda_names, se_lambda_names)
  ))
}

# Hierarchical MIMIC model with three repeat mass measurements.
# Mass =~ M1+M2+M3 recovers latent true mass (absorbs sigma_me into M1-M3 residuals).
# Size =~ 1*Mass+App1+App2 then cleanly estimates direct Temp effects.
# Empirically validated to recover true estimates with sigma_me = 1 g.
fit_lavaan_sem_3mass <- function(df) {
  df_s <- df %>%
    dplyr::select(M1_log, M2_log, M3_log, Append_log, Tarsus_log, Temp_inc) %>%
    rename(M1 = M1_log, M2 = M2_log, M3 = M3_log,
           App1 = Append_log, App2 = Tarsus_log, Temp = Temp_inc) %>%
    mutate(across(everything(), \(x) as.numeric(scale(x))))

  model <- '
    Mass =~ M1 + M2 + M3
    Size =~ 1*Mass + App1 + App2
    Size ~ Temp
    App1 ~ Temp
    App2 ~ Temp
  '

  fit <- tryCatch(
    lavaan::sem(model, data = df_s, estimator = "ML"),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(lavaan::sem(model, data = df_s, estimator = "ML"))
  )

  if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
    return(tibble(coef_sem3_wing = NA_real_, coef_sem3_tarsus = NA_real_))
  }

  ests <- lavaan::parameterEstimates(fit) %>% filter(op == "~", rhs == "Temp")
  w <- ests %>% filter(lhs == "App1") %>% pull(est)
  t <- ests %>% filter(lhs == "App2") %>% pull(est)
  if (abs(w) > 1.5 || abs(t) > 1.5) {
    return(tibble(coef_sem3_wing = NA_real_, coef_sem3_tarsus = NA_real_))
  }
  tibble(coef_sem3_wing = w, coef_sem3_tarsus = t)
}

# ── Measurement error helpers ─────────────────────────────────────────────────
#
# Each classical method appears twice: an ORACLE variant using the noiseless latent
# log-mass (Mass_log) and an M1 variant using a single noisy weighing (M1_log). The
# oracle isolates the bias intrinsic to the estimator from the extra cost of mass
# measurement error, exactly as ryding_oracle / ryding_m1 do for multiple regression.

# Appendage/mass ratio. `mass_col` selects the oracle (Mass_log) or noisy (M1_log) mass.
ratio_by_mass <- function(df, mass_col, suffix) {
  df_s <- df %>%
    mutate(app_mass_w = Append_log / .data[[mass_col]],
           app_mass_t = Tarsus_log / .data[[mass_col]]) %>%
    mutate(across(where(is.numeric), \(x) as.numeric(scale(x))))
  w <- lm(app_mass_w ~ Temp_inc, data = df_s)
  t <- lm(app_mass_t ~ Temp_inc, data = df_s)
  tibble(!!paste0("coef_ratio_", suffix, "_wing")   := coef(w)["Temp_inc"],
         !!paste0("coef_ratio_", suffix, "_tarsus") := coef(t)["Temp_inc"])
}

ratio_oracle <- function(df) ratio_by_mass(df, "Mass_log", "oracle")
ratio_m1     <- function(df) ratio_by_mass(df, "M1_log",   "m1")

# SLI (estimated slope), via sliR::calc_sli() on the raw scale (see apply_methods() above
# for why logging it afterwards leaves the Temp_inc slope unchanged). `mass_col` selects
# the log mass column (Mass_log/M1_log); its raw-scale counterpart (Mass/M1) already
# exists in df. Fitting the SMA to the noisy mass (M1) is what a field study would
# actually do, so the attenuated slope is part of the M1 variant.
sli_est_by_mass <- function(df, mass_col, suffix) {
  mass_raw <- sub("_log$", "", mass_col)
  mass     <- df[[mass_col]]
  b_w  <- coef(smatr::sma(df$Append_log ~ mass))["slope"]
  b_t  <- coef(smatr::sma(df$Tarsus_log ~ mass))["slope"]
  df_s <- df %>%
    sliR::calc_sli(Append = Append, Mass = .data[[mass_raw]], b_sli = b_w, rename_col = "sli_w") %>%
    sliR::calc_sli(Append = Tarsus, Mass = .data[[mass_raw]], b_sli = b_t, rename_col = "sli_t") %>%
    mutate(sli_w = log(sli_w),
           sli_t = log(sli_t)) %>%
    mutate(across(where(is.numeric), \(x) as.numeric(scale(x))))
  w <- lm(sli_w ~ Temp_inc, data = df_s)
  t <- lm(sli_t ~ Temp_inc, data = df_s)
  tibble(!!paste0("coef_sli_est_", suffix, "_wing")   := coef(w)["Temp_inc"],
         !!paste0("coef_sli_est_", suffix, "_tarsus") := coef(t)["Temp_inc"])
}

sli_est_oracle <- function(df) sli_est_by_mass(df, "Mass_log", "oracle")
sli_est_m1     <- function(df) sli_est_by_mass(df, "M1_log",   "m1")

# ── Misspecification simulation: is the SEM's error bias or imprecision? ──────
#
# Supports the Supporting Information analysis separating the SEM's two failure modes.
# Wide intervals are IMPRECISION (small N, or mass a labile size proxy) and shrink as
# 1/sqrt(N). A second latent factor shared by the appendages but not mass is
# MISSPECIFICATION, and its bias is constant in N. With only 3 indicators the
# single-factor model is just-identified (df = 0), so the two are observationally
# equivalent and cannot be told apart by fit statistics alone.

# Latent Size is driven by Temp (Bergmann); mass loads on Size only. Both appendages
# load on Size, on a second latent "Limb" factor with loading gamma, and directly on
# Temp. gamma = 0 recovers the well-specified single-factor model.
gen_sem_data <- function(N, B_temp = 0.3, beta_size = -0.6, lambda = 1, gamma = 0,
                         sigma_size = 1, sigma_limb = 1, sigma_mass = 0.5, sigma_app = 0.5) {
  Temp <- rnorm(N)
  Size <- rnorm(N, beta_size * Temp, sigma_size)
  Limb <- rnorm(N, 0, sigma_limb)
  tibble::tibble(Temp, Size,
    log_mass   = 1 * Size + rnorm(N, 0, sigma_mass),
    log_wing   = B_temp * Temp + lambda * Size + gamma * Limb + rnorm(N, 0, sigma_app),
    log_tarsus = B_temp * Temp + lambda * Size + gamma * Limb + rnorm(N, 0, sigma_app))
}

# Partial correlation cor(wing, tarsus | mass) — the single-factor strain signal.
# Both extra mass noise and a second factor inflate it, which is why they cannot be
# distinguished with three indicators.
pcor_wt_given_m <- function(d) {
  r   <- cor(d[c("log_wing", "log_tarsus", "log_mass")])
  rw  <- r[1, 3]; rt <- r[2, 3]; rwt <- r[1, 2]
  (rwt - rw * rt) / sqrt((1 - rw^2) * (1 - rt^2))
}

# Fit the single-factor MIMIC SEM alongside an ORACLE regression on the true latent
# Size. The oracle's Temp coefficient is the estimand, so bias = SEM estimate - oracle.
# No plausibility guard here: the point is to let any bias show.
fit_sf_vs_oracle <- function(d) {
  z <- d %>% dplyr::mutate(dplyr::across(c(Size, log_mass, log_wing, log_tarsus, Temp),
                                         \(x) as.numeric(scale(x))))
  oracle <- coef(lm(log_wing ~ Size + Temp, data = z))["Temp"]
  model <- 'Sz =~ 1*log_mass + log_wing + log_tarsus
            Sz ~ Temp
            log_wing ~ Temp
            log_tarsus ~ Temp'
  f <- tryCatch(suppressWarnings(lavaan::sem(model, z)), error = function(e) NULL)
  if (is.null(f) || !lavaan::lavInspect(f, "converged"))
    return(tibble::tibble(est = NA_real_, se = NA_real_, lam = NA_real_,
                          se_lam = NA_real_, oracle = oracle))
  pe <- lavaan::parameterEstimates(f)
  tibble::tibble(
    est    = pe$est[pe$lhs == "log_wing" & pe$op == "~"  & pe$rhs == "Temp"],
    se     = pe$se [pe$lhs == "log_wing" & pe$op == "~"  & pe$rhs == "Temp"],
    lam    = pe$est[pe$op == "=~" & pe$rhs == "log_wing"],
    se_lam = pe$se [pe$op == "=~" & pe$rhs == "log_wing"],
    oracle = oracle)
}

# Ryding using the true noiseless mass — oracle lower bound on collider bias.
ryding_oracle <- function(df) {
  df_s <- df %>%
    mutate(across(c(Mass_log, Append_log, Tarsus_log, Temp_inc),
                  \(x) as.numeric(scale(x))))
  w <- lm(Append_log ~ Mass_log + Temp_inc, data = df_s)
  t <- lm(Tarsus_log ~ Mass_log + Temp_inc, data = df_s)
  tibble(coef_ryding_oracle_wing   = coef(w)["Temp_inc"],
         coef_ryding_oracle_tarsus = coef(t)["Temp_inc"])
}

# Ryding using M1 only (single noisy mass measurement).
ryding_m1 <- function(df) {
  df_s <- df %>%
    mutate(across(c(M1_log, Append_log, Tarsus_log, Temp_inc),
                  \(x) as.numeric(scale(x))))
  w <- lm(Append_log ~ M1_log + Temp_inc, data = df_s)
  t <- lm(Tarsus_log ~ M1_log + Temp_inc, data = df_s)
  tibble(coef_ryding_m1_wing   = coef(w)["Temp_inc"],
         coef_ryding_m1_tarsus = coef(t)["Temp_inc"])
}

# ── Bayesian latent SEM (blavaan) ─────────────────────────────────────────────
#
# Same MIMIC model as fit_lavaan_sem() but fit with blavaan (a Stan backend; the
# default target marginalizes the latent scores, like latent_size_3ind_marginal.stan).
# Two things this buys on hard empirical data where the lavaan fit misbehaves:
#   (1) blavaan's residual-variance priors have positive support, so it CANNOT
#       return the negative-variance (Heywood) solutions lavaan silently produced;
#   (2) a prior on the appendage loadings regularizes them toward isometry,
#       preventing the loading blow-up that inflates the direct temperature effect.
#
# `loading_prior` is the tunable knob: a full Stan prior string applied to every
# appendage loading (Mass stays the fixed = 1 anchor). Tighter SD pulls the
# loadings harder toward 1 (equal, isometric loadings); a wide SD reproduces the
# unregularized lavaan behaviour. Examples: "normal(1,0.5)" (weak), "normal(1,0.25)"
# (moderate), "normal(1,0.15)" (strong / near-isometric).
#
# Returns the same columns as fit_lavaan_sem() (coef_sem_/se_sem_/lambda_sem_/
# se_lambda_sem_ per label) PLUS lci_/uci_ per label (95% highest-posterior-density
# interval of the direct temperature effect). se_ here is the HPD-implied posterior
# SD, (uci - lci) / 3.92.
fit_blavaan_sem <- function(df,
                            mass_name     = "Mass_log",
                            append_names  = c("Append_log", "Tarsus_log"),
                            temp_name     = "Temp_inc",
                            labels        = c("wing", "tarsus"),
                            size_covs     = character(0),
                            loading_prior = "normal(1,0.25)",
                            n_chains = 4, burnin = 1000, sample = 1000, seed = 1) {
  stopifnot(requireNamespace("blavaan", quietly = TRUE))
  app_ids <- paste0("App", seq_along(append_names))
  col_map <- c(setNames(mass_name, "Mass"),
               setNames(append_names, app_ids),
               setNames(temp_name, "Temp"))

  df_s <- df %>%
    dplyr::select(dplyr::all_of(c(mass_name, append_names, temp_name, size_covs))) %>%
    dplyr::rename(dplyr::all_of(col_map)) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), as.numeric))

  # lavaan/blavaan require numeric predictors: encode character/factor covs as 0/1
  for (cov in size_covs) {
    if (!is.numeric(df_s[[cov]])) {
      ref <- sort(unique(stats::na.omit(df_s[[cov]])))[1]
      df_s[[cov]] <- as.numeric(df_s[[cov]] != ref)
    }
  }

  # Attach the loading prior to each appendage loading; Mass anchors the scale.
  app_terms <- paste0("prior(\"", loading_prior, "\")*", app_ids)
  size_rhs  <- paste(c("Temp", size_covs), collapse = " + ")
  model <- paste0(
    "Size =~ 1*Mass + ", paste(app_terms, collapse = " + "), "\n",
    "Size ~ ", size_rhs, "\n",
    paste(paste0(app_ids, " ~ Temp"), collapse = "\n")
  )

  ret_names <- c(paste0("coef_sem_", labels), paste0("se_sem_", labels),
                 paste0("lambda_sem_", labels), paste0("se_lambda_sem_", labels),
                 paste0("lci_", labels), paste0("uci_", labels))
  na_ret <- tibble::as_tibble(setNames(as.list(rep(NA_real_, length(ret_names))), ret_names))

  fit <- tryCatch(
    suppressWarnings(suppressMessages(
      blavaan::bsem(model, data = df_s, n.chains = n_chains,
                    burnin = burnin, sample = sample,
                    bcontrol = list(cores = n_chains, refresh = 0), seed = seed))),
    error = function(e) NULL)
  if (is.null(fit)) return(na_ret)

  cf  <- tryCatch(lavaan::coef(fit),                error = function(e) NULL)  # posterior means
  hpd <- tryCatch(blavaan::blavInspect(fit, "hpd"), error = function(e) NULL)  # 95% HPD intervals
  if (is.null(cf) || is.null(hpd)) return(na_ret)

  getp  <- function(p)    if (p %in% names(cf))     cf[[p]]   else NA_real_
  getci <- function(p, k) if (p %in% rownames(hpd)) hpd[p, k] else NA_real_
  coefs   <- vapply(app_ids, \(a) getp(paste0(a, "~Temp")),      numeric(1))
  lcis    <- vapply(app_ids, \(a) getci(paste0(a, "~Temp"), 1),  numeric(1))
  ucis    <- vapply(app_ids, \(a) getci(paste0(a, "~Temp"), 2),  numeric(1))
  lambdas <- vapply(app_ids, \(a) getp(paste0("Size=~", a)),     numeric(1))
  lam_lo  <- vapply(app_ids, \(a) getci(paste0("Size=~", a), 1), numeric(1))
  lam_hi  <- vapply(app_ids, \(a) getci(paste0("Size=~", a), 2), numeric(1))

  tibble::as_tibble(setNames(
    as.list(c(coefs, (ucis - lcis) / 3.92, lambdas, (lam_hi - lam_lo) / 3.92, lcis, ucis)),
    ret_names))
}
