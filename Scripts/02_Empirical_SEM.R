## 02_Empirical_SEM.R ----
## Per-species latent SEM vs. ratio / SLI / multiple regression on the three empirical
## datasets, plus the identifiability diagnostics that explain where the SEM struggles.
##
## THIS SCRIPT NEEDS THE RAW DATA, which lives outside the repo (see the absolute paths
## below). It is the only step in the pipeline that does. Its outputs are committed, so
## 04_Empirical_plots.R rebuilds every empirical figure without it. See README.
##
## Run order: this script -> 03_Hierarchical_SEM.R -> 04_Empirical_plots.R
## Key output: Derived/blavaan_empirical_estimates.csv (requires blavaan_full = TRUE)
##
## Self-contained: this script loads each raw dataset directly (absolute paths),
## does only the minimal cleaning the analysis needs, then applies three
## approaches per species — the appendage/mass ratio, mass-as-covariate multiple
## regression (Ryding), and the MIMIC latent-size SEM (fit_lavaan_sem). It does
## NOT source the large SLI shape scripts (Nightjar_shape.R / Weeks_2020_ral.R /
## Atlantic_birds_shape.R, now in the SLI-allens-rule repo): those compute SLI,
## residual and plotting machinery we do not need here and reuse object names
## (Spp_metadata2, sig_age_any, ...) that collide across datasets.
##
## Diagnostic goal: the latent SEM recovers the true direct temperature effect on
## SIMULATED causal data but behaves poorly on the EMPIRICAL data. This script
## shows WHY. It (1) runs a simulation sanity check to confirm the fitting code is
## correct; (2) fits all three empirical datasets through one shared driver and
## compares the SEM against the ratio and multiple-regression approaches; and
## (3) prints the collider-relevant diagnostics — factor loading lambda and its
## SE, the allometric coupling r(Mass, appendage), and Bergmann strength
## r(Mass, Temp) — so the cause of the SEM's imprecision is visible in one place.
## It also fits each species with and without Age/Sex covariates to test whether
## covariate handling is responsible.
##
## Latent SEM (per species, all variables z-scored first):
##   Size =~ 1*Mass + Wing + {Tail | Tarsus}   # Mass anchors latent Size (loading = 1)
##   Size ~ Temp (+ Age + Sex)                 # Bergmann path + optional covariates
##   Wing ~ Temp                               # direct Allen's effect (target estimand)
##   {Tail|Tarsus} ~ Temp                      # direct Allen's effect
## With three indicators the model is just-identified (df = 0): it can only work
## when Wing/Tail/Tarsus genuinely load on the Mass-anchored size factor.

# Setup ----
library(tidyverse)
library(lavaan)
library(sliR)
source("Scripts/Key_causal_fns.R")   # provides fit_lavaan_sem() and gen_causal_data()

## Data paths (absolute; these datasets live outside this repo)
sli_data     <- "/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Academia/Methods_papers/SLI-allens-rule/Data"
ext_dir      <- "/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Academia/Datasets_external"
nightjar_rds <- file.path(sli_data, "Nightjar_temp.rds")                 # nightjar raw + cached WorldClim B.Temp
weeks_csv    <- file.path(ext_dir, "Weeks_etal_2020_Data.csv")
atlantic_csv <- file.path(ext_dir, "Ecology/Atlantic_bird_traits/ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv")

## Run parameters
min_n_obs   <- 150   # minimum complete-case observations per species
min_n_group <- 100   # minimum individuals per Age/Sex level to use it as a covariate

## Section switches (guarded so a driver script can override them before sourcing).
## The Bayesian section is off by default because blavaan fits are slow (~30-60s each).
if (!exists("run_sample_size_sim")) run_sample_size_sim <- TRUE   # nightjar N experiment (~1 min)
if (!exists("run_blavaan"))         run_blavaan         <- TRUE  # Bayesian (blavaan) section
if (!exists("blavaan_prior_grid"))  blavaan_prior_grid  <- c("normal(1,0.5)", "normal(1,0.25)", "normal(1,0.15)")  # sensitivity strip
if (!exists("blavaan_prior"))       blavaan_prior       <- "normal(1,0.25)"  # prior for the full empirical refit
if (!exists("blavaan_full"))        blavaan_full        <- FALSE  # TRUE: refit ALL species (very slow, ~30-60 min)
if (!exists("run_diagnostics"))     run_diagnostics     <- TRUE  # identifiability diagnostics (needs the blavaan CSV)

## Approach ordering, colours, and Bergmann-direction ordering shared by the
## comparison tables and every figure.
approach_levels <- c("Ratio", "Sli_iso", "Sli_est", "Ryding", "SEM")
approach_cols   <- c(Ratio = "#E69F00", Sli_iso = "#CC79A7", Sli_est = "#D55E00",
                     Ryding = "#0072B2", SEM = "#009E73")
approach_labs   <- c(Ratio = "Wing / Mass", Sli_iso = "SLI isometry",
                     Sli_est = "SLI estimated", Ryding = "Mass as covariate", SEM = "Latent SEM")
direction_levels <- c("Bergmann's", "Inverse Bergmann's", "Mixed - Wingier",
                      "Mixed - Fatter", "Stable")
study_colors    <- c("Nightjar" = "#E41A1C", "Weeks (2020)" = "#377EB8",
                     "Atlantic birds" = "#4DAF4A")


# Shared helpers ----

## Recode obvious "unknown" codings to NA so a covariate only ever contributes
## its two real levels (lavaan / lm drop NA rows listwise when the covariate is used).
clean_group <- function(x) {
  x <- as.character(x)
  if_else(x %in% c("Unknown", "Unk", "U", "unknown", "unk", "u", ""), NA_character_, x)
}

## Build the per-species z-scored list common to every dataset. `df` must already
## carry log_mass, the appendage log columns, the ratio columns, the predictor,
## and (optionally) Age/Sex. Numeric columns are z-scored within species (plain
## numeric vectors, not the 1-column matrices scale() returns); Age/Sex stay
## character. Species with fewer than min_n_obs complete rows are dropped.
make_species_list <- function(df, group_col) {
  df %>%
    add_count(.data[[group_col]], name = ".n") %>%
    filter(.n >= min_n_obs) %>%
    dplyr::select(-.n) %>%
    group_split(.data[[group_col]]) %>%
    set_names(map_chr(., \(d) as.character(unique(d[[group_col]])))) %>%
    map(\(d) mutate(d, across(where(is.numeric), \(x) as.numeric(scale(x)))))
}

## Decide which Age/Sex columns to use as covariates for a species: keep a column
## only if it has >= 2 levels each with at least min_n_group individuals.
pick_covs <- function(df, cov_cols) {
  keep <- character(0)
  for (v in cov_cols) {
    if (!v %in% names(df)) next
    lvl_ok <- df %>% filter(!is.na(.data[[v]])) %>% count(.data[[v]]) %>%
      filter(n >= min_n_group)
    if (nrow(lvl_ok) >= 2) keep <- c(keep, v)
  }
  keep
}

## Extract the predictor coefficient from a fitted lm as a one-row tidy tibble.
lm_coef <- function(mod, term, species, appendage, approach) {
  s <- summary(mod)$coefficients
  tibble(species = species, appendage = appendage, Approach = approach,
         estimate = s[term, "Estimate"], std.error = s[term, "Std. Error"],
         lambda = NA_real_, se_lambda = NA_real_)
}

## Add the Standardized Length Index columns for each appendage, BEFORE z-scoring, via
## sliR::calc_sli(). calc_sli works on the raw scale (SLI = App*(M0/Mass)^b); log_cols
## are already logged, so raw values are reconstructed with exp() rather than requiring
## the original raw columns (some dataset pipelines above select() them away before this
## is called). log(SLI) = log(App) - b*log(Mass) + b*log(M0): the additive b*log(M0) term
## is constant within a group and drops out under the z-scoring/regression that follows,
## so it does not change results. sli_iso uses the isometric b = 0.33; sli_est uses a
## per-group SMA slope of log(Append) ~ log(Mass), fit manually with smatr (as
## apply_methods()/sli_est_by_mass() in Key_causal_fns.R do) rather than calc_sli's
## control= argument: control= fits ALL groups in one joint `mass * group` SMA, which
## errors outright on a single-group df (the simulation sanity check has only one
## "species") and on any group with < 3 observations (real singleton-species rows in
## the empirical data, filtered out later by min_n_obs but still present here).
add_sli_cols <- function(df, group_col, labels, log_cols, b_iso = 0.33) {
  for (i in seq_along(labels)) {
    a <- log_cols[i]; lab <- labels[i]
    est_col <- paste0("sli_est_", lab)
    df <- df %>%
      mutate(.sli_app = exp(.data[[a]]), .sli_mass = exp(log_mass)) %>%
      sliR::calc_sli(Append = .sli_app, Mass = .sli_mass, b_sli = b_iso,
                     rename_col = paste0("sli_iso_", lab)) %>%
      group_by(.data[[group_col]]) %>%
      group_modify(\(d, ...) {
        b_sma <- tryCatch(
          coef(smatr::sma(log(.sli_app) ~ log(.sli_mass), data = d))["slope"],
          error = function(e) NA_real_)
        if (is.na(b_sma)) {
          d[[est_col]] <- NA_real_
          d
        } else {
          sliR::calc_sli(d, Append = .sli_app, Mass = .sli_mass, b_sli = b_sma, rename_col = est_col)
        }
      }) %>%
      ungroup() %>%
      mutate(across(all_of(c(paste0("sli_iso_", lab), est_col)), log)) %>%
      dplyr::select(-.sli_app, -.sli_mass)
  }
  df
}

## Ratio, SLI-isometry, SLI-estimated (direct regressions of a pre-built size-
## adjusted column on the predictor) and mass-as-covariate multiple regression
## (Ryding). All without Age/Sex covariates, matching the shape-script conventions,
## so they form a clean baseline for the SEM comparison. SLI columns follow the
## naming from add_sli_cols() (sli_iso_<label> / sli_est_<label>).
run_simple_methods <- function(sp_list, temp_name, labels, log_cols, ratio_cols) {
  imap(sp_list, \(df, sp) {
    map(seq_along(labels), \(i) {
      lab    <- labels[i]
      direct <- tibble(Approach = c("Ratio", "Sli_iso", "Sli_est"),
                       col      = c(ratio_cols[i], paste0("sli_iso_", lab), paste0("sli_est_", lab)))
      rows <- pmap(direct, \(Approach, col)
        lm_coef(lm(reformulate(temp_name, col), data = df), temp_name, sp, lab, Approach)
      ) %>% list_rbind()
      ryding_m <- lm(reformulate(c("log_mass", temp_name), log_cols[i]), data = df)
      bind_rows(rows, lm_coef(ryding_m, temp_name, sp, lab, "Ryding"))
    }) %>% list_rbind()
  }) %>% list_rbind()
}

## Classify each species' shape-shifting direction from the marginal trends of Mass
## and Wing on the predictor (temperature or year), following classify_direction()
## in the SLI-allens-rule repo. Bergmann's = mass declines with the predictor and
## wing does not increase; Inverse Bergmann's = the reverse; Mixed = both significant
## and moving the same way (both shrinking = "Fatter"; both growing = "Wingier").
classify_direction <- function(mods_tbl, p_threshold = 0.05) {
  mods_tbl %>%
    mutate(sig = p.value < p_threshold) %>%
    group_by(species) %>%
    summarise(
      n_sig    = sum(sig),
      mass_dir = estimate[dv == "mass"] < 0,   # mass declining
      wing_dir = estimate[dv == "wing"] > 0,   # wing increasing
      mass_sig = sig[dv == "mass"],
      wing_sig = sig[dv == "wing"],
      .groups  = "drop") %>%
    mutate(
      sole_decr = case_when(mass_sig & !wing_sig ~  mass_dir,
                            !mass_sig & wing_sig ~ !wing_dir,
                            .default = NA),
      Direction = case_when(
        n_sig == 0                         ~ "Stable",
        n_sig == 1 & sole_decr             ~ "Bergmann's",
        n_sig == 1 & !sole_decr            ~ "Inverse Bergmann's",
        n_sig == 2 &  mass_dir & !wing_dir ~ "Bergmann's",
        n_sig == 2 & !mass_dir &  wing_dir ~ "Inverse Bergmann's",
        n_sig == 2 &  mass_dir &  wing_dir ~ "Mixed - Wingier",
        n_sig == 2 & !mass_dir & !wing_dir ~ "Mixed - Fatter",
        TRUE ~ "Check")) %>%
    dplyr::select(species, Direction, n_sig)
}

## Per-species marginal trend of a trait on the predictor (estimate + p-value).
classify_bergmann <- function(sp_list, temp_name) {
  trend <- function(df, dv) {
    co <- summary(lm(reformulate(temp_name, dv), data = df))$coefficients
    tibble(estimate = co[temp_name, "Estimate"], p.value = co[temp_name, "Pr(>|t|)"])
  }
  imap(sp_list, \(df, sp) bind_rows(
    trend(df, "log_mass") %>% mutate(dv = "mass"),
    trend(df, "log_wing") %>% mutate(dv = "wing")) %>% mutate(species = sp)
  ) %>% list_rbind() %>% classify_direction()
}

## Fit a latent SEM across a per-species list; one tidy row per appendage.
## `fit_fn` is the fitting function (fit_lavaan_sem or fit_blavaan_sem); extra
## arguments (e.g. loading_prior for blavaan) pass through via `...`. covs_fn(df)
## supplies size_covs per species (character(0) = none), so the same runner does
## the no-covariate and covariate passes. lci/uci carry the blavaan HPD interval
## when present (NA for lavaan, where CIs are formed later from estimate +/- SE).
run_sem_list <- function(sp_list, temp_name, labels, log_cols,
                         covs_fn = function(df) character(0),
                         fit_fn = fit_lavaan_sem, ...) {
  imap(sp_list, \(df, sp) {
    covs <- covs_fn(df)
    res  <- fit_fn(df,
      mass_name = "log_mass", append_names = log_cols,
      temp_name = temp_name, labels = labels, size_covs = covs, ...)
    map(labels, \(lab) {
      grab <- function(prefix) { v <- res[[paste0(prefix, lab)]]; if (is.null(v)) NA_real_ else v }
      tibble(
        species = sp, appendage = lab, Approach = "SEM",
        estimate  = grab("coef_sem_"),   std.error = grab("se_sem_"),
        lambda    = grab("lambda_sem_"), se_lambda = grab("se_lambda_sem_"),
        lci = grab("lci_"), uci = grab("uci_"))
    }) %>% list_rbind()
  }) %>% list_rbind()
}

## Data-side diagnostics that do not depend on any model converging: sample size,
## allometric coupling of each appendage to mass, and Bergmann strength (all on
## z-scored data, so the coefficients are correlations).
diag_list <- function(sp_list, temp_name, labels, log_cols) {
  imap(sp_list, \(df, sp) {
    tibble(
      species     = sp,
      n           = nrow(df),
      r_mass_temp = cor(df$log_mass, df[[temp_name]],  use = "complete.obs"),  # Bergmann's
      !!paste0("r_mass_", labels[1]) := cor(df$log_mass, df[[log_cols[1]]], use = "complete.obs"),
      !!paste0("r_mass_", labels[2]) := cor(df$log_mass, df[[log_cols[2]]], use = "complete.obs")
    )
  }) %>% list_rbind()
}

## One driver for all three datasets. Returns the diagnostics, the SEM fit with
## and without Age/Sex covariates, and the stacked three-approach comparison.
analyse_dataset <- function(sp_list, temp_name, labels, log_cols, ratio_cols, cov_cols, study) {
  sem_nocov <- run_sem_list(sp_list, temp_name, labels, log_cols)
  sem_cov   <- run_sem_list(sp_list, temp_name, labels, log_cols,
                            covs_fn = \(df) pick_covs(df, cov_cols))
  simple    <- run_simple_methods(sp_list, temp_name, labels, log_cols, ratio_cols)
  direction <- classify_bergmann(sp_list, temp_name)

  ## Use the HPD interval when a fit supplied one (blavaan); otherwise form the CI
  ## from estimate +/- 1.96 SE (lavaan SEM, and the ratio / Ryding lm fits).
  comparison <- bind_rows(sem_nocov, simple) %>%
    mutate(LCI95 = coalesce(lci, estimate - 1.96 * std.error),
           UCI95 = coalesce(uci, estimate + 1.96 * std.error),
           Approach = factor(Approach, levels = approach_levels)) %>%
    left_join(direction, by = "species")

  list(study = study,
       diag       = diag_list(sp_list, temp_name, labels, log_cols),
       direction  = direction,
       sem_nocov  = sem_nocov,
       sem_cov    = sem_cov,
       comparison = comparison)
}

## Per-approach precision summary for the wing (the head-to-head the paper cares
## about): how often each method returns a usable estimate and how wide it is.
method_summary <- function(comparison) {
  comparison %>%
    filter(appendage == "wing") %>%
    group_by(Approach) %>%
    summarise(n_species   = n(),
              n_usable     = sum(!is.na(estimate)),
              median_est   = median(estimate,  na.rm = TRUE),
              median_se    = median(std.error, na.rm = TRUE),
              max_se       = max(std.error,    na.rm = TRUE),
              .groups = "drop")
}


# Simulation sanity check ----
## Generate data from the paper's causal DAG and push it through the SAME shared
## functions used on the empirical data (run_simple_methods + run_sem_list). If the
## SEM recovers the direct effect here — matching the oracle regression that uses
## the TRUE latent Size, while mass-as-covariate (Ryding) shows the expected
## collider bias — then the fitting code is correct and any empirical failure is a
## data / causal-assumption problem, not a bug.
set.seed(1)
sim_true_bw <- 0.05                       # true direct Temp -> Wing effect (log mm per Temp unit)
## lambda_wing/tarsus = 0.9 gives a strong shared size factor (r(Mass,Wing) ~ 0.63,
## r(Wing,Tarsus) ~ 0.61); beta_size = -0.42 gives strong Bergmann's (r(Mass,Temp) ~ -0.64).
## This is the well-behaved causal regime the SEM is designed for.
sim_raw <- gen_causal_data(N = 1500, B_temp_w = sim_true_bw, B_temp_t = sim_true_bw,
                           beta_size = -0.42, lambda_wing = 0.9, lambda_tarsus = 0.9,
                           sigma_wing = 0.10, sigma_tarsus = 0.10)

## Rename to the shared column convention and add the ratio columns, so the sim is
## just another "dataset" (a single species) fed through the common machinery.
sim_df <- sim_raw %>%
  transmute(species_ = "sim", Temp_inc, Size,
            log_mass = Mass_log, log_wing = Append_log, log_tarsus = Tarsus_log,
            wing_mass = Append_log - Mass_log, tarsus_mass = Tarsus_log - Mass_log) %>%
  add_sli_cols("species_", c("wing", "tarsus"), c("log_wing", "log_tarsus"))
sim_list <- make_species_list(sim_df, "species_")

sim_simple <- run_simple_methods(sim_list, "Temp_inc", c("wing", "tarsus"),
                                 c("log_wing", "log_tarsus"), c("wing_mass", "tarsus_mass"))
sim_sem    <- run_sem_list(sim_list, "Temp_inc", c("wing", "tarsus"),
                           c("log_wing", "log_tarsus"))
## Oracle: regress Wing on the TRUE latent Size + Temp (z-scored). Its Temp coef is
## the direct-effect estimand on the standardized scale — the target for the SEM.
sim_oracle <- coef(summary(lm(log_wing ~ Size + Temp_inc, data = sim_list[[1]])))["Temp_inc", ]

cat("\n===== SIMULATION SANITY CHECK (causal DAG, one species, N=1500) =====\n")
cat(sprintf("True B_temp_w = %.2f log mm/unit.  Bergmann r(Mass,Temp) = %+.2f,  r(Mass,Wing) = %.2f\n",
            sim_true_bw, cor(sim_raw$Mass_log, sim_raw$Temp_inc),
            cor(sim_raw$Mass_log, sim_raw$Append_log)))
cat("\nWing, standardized Temp coefficient by approach:\n")
print(bind_rows(
  sim_simple %>% filter(appendage == "wing") %>% dplyr::select(Approach, estimate, std.error, lambda, se_lambda),
  sim_sem    %>% filter(appendage == "wing") %>% dplyr::select(Approach, estimate, std.error, lambda, se_lambda),
  tibble(Approach = "Oracle(trueSize)", estimate = sim_oracle["Estimate"],
         std.error = sim_oracle["Std. Error"], lambda = NA_real_, se_lambda = NA_real_)
))
cat("=> SEM should sit near the oracle with a small SE and a well-identified loading;\n")
cat("   Ryding (mass as covariate) should be biased toward / below zero (collider bias).\n")


# Nightjar ----
## Appendages: Wing + Tail (no Tarsus in this dataset). Mass anchors latent Size.
## B.Temp is WorldClim annual mean temperature (positive coef = Allen's rule).
nj_df <- readRDS(nightjar_rds) %>%
  dplyr::select(Species, Wing = Wing.comb, Mass = Mass.comb, Tail = Tail.comb,
                B.Temp, Age, Sex) %>%
  drop_na(Wing, Mass, Tail) %>%
  mutate(log_wing = log(Wing), log_mass = log(Mass), log_tail = log(Tail),
         wing_mass = log_wing - log_mass, tail_mass = log_tail - log_mass,
         Age = clean_group(Age), Sex = clean_group(Sex)) %>%
  add_sli_cols("Species", c("wing", "tail"), c("log_wing", "log_tail"))

nj_list <- make_species_list(nj_df, "Species")
nj <- analyse_dataset(
  nj_list,
  temp_name  = "B.Temp",
  labels     = c("wing", "tail"),
  log_cols   = c("log_wing", "log_tail"),
  ratio_cols = c("wing_mass", "tail_mass"),
  cov_cols   = c("Age", "Sex"),
  study      = "Nightjar")


# Weeks (2020) ----
## Appendages: Wing + Tarsus. Predictor is Year (temporal shape-shifting;
## positive = trait increasing over time).
weeks_df <- read_csv(weeks_csv, skip = 2, show_col_types = FALSE) %>%
  rename(species_ = Taxon) %>%
  mutate(year = as.integer(format(lubridate::mdy(Date), "%Y")),
         log_mass = log(Mass), log_wing = log(Wing), log_tarsus = log(Tarsus),
         wing_mass = log_wing - log_mass, tarsus_mass = log_tarsus - log_mass,
         Age = clean_group(Age), Sex = clean_group(Sex)) %>%
  dplyr::select(species_, year, log_mass, log_wing, log_tarsus,
                wing_mass, tarsus_mass, Age, Sex) %>%
  drop_na(log_mass, log_wing, log_tarsus) %>%
  filter(year > 1978) %>%
  add_sli_cols("species_", c("wing", "tarsus"), c("log_wing", "log_tarsus"))

weeks_list <- make_species_list(weeks_df, "species_")
weeks <- analyse_dataset(
  weeks_list,
  temp_name  = "year",
  labels     = c("wing", "tarsus"),
  log_cols   = c("log_wing", "log_tarsus"),
  ratio_cols = c("wing_mass", "tarsus_mass"),
  cov_cols   = c("Age", "Sex"),
  study      = "Weeks (2020)")


# Atlantic birds ----
## Appendages: Wing + Tarsus. Predictor is B.Tavg (annual mean temperature).
atlantic_df <- read_csv(atlantic_csv, show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  rename(mass = body_mass_g, B.Tavg = annual_mean_temperature) %>%
  mutate(species_ = str_replace(binomial, " ", "_"),
         wing   = coalesce(wing_length_mm, wing_length_left_mm, wing_length_right_mm),
         tarsus = coalesce(tarsus_length_mm, tarsus_length_right_mm, tarsus_length_left_mm),
         log_mass = log(mass), log_wing = log(wing), log_tarsus = log(tarsus),
         wing_mass = log_wing - log_mass, tarsus_mass = log_tarsus - log_mass,
         Age = clean_group(age), Sex = clean_group(sex)) %>%
  dplyr::select(species_, year, B.Tavg, log_mass, log_wing, log_tarsus,
                wing_mass, tarsus_mass, Age, Sex) %>%
  drop_na(log_mass, log_wing, log_tarsus) %>%
  filter(year > 1990, !is.na(B.Tavg)) %>%
  add_sli_cols("species_", c("wing", "tarsus"), c("log_wing", "log_tarsus"))

atlantic_list <- make_species_list(atlantic_df, "species_")
atlantic <- analyse_dataset(
  atlantic_list,
  temp_name  = "B.Tavg",
  labels     = c("wing", "tarsus"),
  log_cols   = c("log_wing", "log_tarsus"),
  ratio_cols = c("wing_mass", "tarsus_mass"),
  cov_cols   = c("Age", "Sex"),
  study      = "Atlantic birds")


# Report ----
datasets <- list(nj, weeks, atlantic)

for (d in datasets) {
  cat(sprintf("\n===== %s (%d species) =====\n", d$study, nrow(d$diag)))
  cat("\n-- per-approach precision (wing) --\n")
  print(method_summary(d$comparison))
  cat("\n-- wing estimates by approach --\n")
  print(d$comparison %>% filter(appendage == "wing") %>%
          dplyr::select(species, Approach, estimate, std.error, lambda, se_lambda) %>%
          arrange(species, Approach), n = Inf)
  cat("\n-- Age/Sex effect on SEM wing coef (no-cov vs cov) --\n")
  print(d$sem_nocov %>% filter(appendage == "wing") %>%
          dplyr::select(species, est_nocov = estimate, se_nocov = std.error) %>%
          left_join(d$sem_cov %>% filter(appendage == "wing") %>%
                      dplyr::select(species, est_cov = estimate, se_cov = std.error),
                    by = "species"), n = Inf)
}


# Cross-dataset summary ----
## Link each species' SEM behaviour (converged? coefficient SE, loading SE) to its
## data-side diagnostics (allometric coupling, Bergmann strength). Poor SEM
## precision should track weak coupling (low / uncertain lambda) and/or weak
## Bergmann's — NOT anything about the code or the ratio / regression baselines.
overall <- map(datasets, \(d) {
  d$sem_nocov %>% filter(appendage == "wing") %>%
    left_join(d$diag, by = "species") %>%
    transmute(study = d$study, species,
              converged = !is.na(estimate), n,
              se_wing = std.error, lambda_wing = lambda, se_lambda_wing = se_lambda,
              r_mass_wing, bergmann = r_mass_temp)
}) %>% list_rbind()

cat("\n===== CROSS-DATASET SUMMARY (SEM wing, no covariates) =====\n")
cat(sprintf("Non-convergence / dropped: %d of %d species\n",
            sum(!overall$converged), nrow(overall)))

cat("\nMedian SEM wing-coefficient SE by dataset (converged species):\n")
print(overall %>% filter(converged) %>% group_by(study) %>%
        summarise(n_species = n(),
                  median_se_wing    = median(se_wing),
                  median_lambda     = median(abs(lambda_wing)),
                  median_se_lambda  = median(se_lambda_wing),
                  median_bergmann   = median(bergmann), .groups = "drop"))

cat("\nWorst 15 species by SEM wing-coefficient SE:\n")
print(overall %>% filter(converged) %>% arrange(desc(se_wing)) %>% head(15))

ov <- overall %>% filter(converged)
cat("\nAcross converged species, SEM wing-coefficient SE correlates with:\n")
cat(sprintf("  SE(lambda_wing)            %+.2f   (identification of the size factor)\n",
            cor(ov$se_wing, ov$se_lambda_wing, use = "complete.obs")))
cat(sprintf("  r(Mass, Wing)              %+.2f   (allometric coupling)\n",
            cor(ov$se_wing, ov$r_mass_wing,    use = "complete.obs")))
cat(sprintf("  Bergmann r(Mass, Temp)     %+.2f   (strength of the collider path)\n",
            cor(ov$se_wing, ov$bergmann,       use = "complete.obs")))


# Figure: empirical wing beta by approach, grouped by Bergmann direction ----
## One row per dataset; within a dataset, species are split into panels by their
## Bergmann-rule direction (columns), so the methods can be compared within each
## shape-shifting class. Five approaches per species: Ratio, SLI isometry, SLI
## estimated, mass-as-covariate (Ryding), and the latent SEM. Very imprecise SEM
## fits (SE >= 0.5) are dropped so the informative estimates stay legible.
library(cowplot)
ggplot2::theme_set(theme_cowplot(font_size = 11))

plot_estimates <- function(cmp, title) {
  df <- cmp %>% filter(appendage == "wing", std.error < 0.5,
                       Direction %in% direction_levels)
  ord <- df %>% group_by(species) %>% summarise(m = mean(estimate), .groups = "drop") %>%
    arrange(m) %>% pull(species)
  df %>%
    mutate(species   = factor(species, levels = ord),
           Direction = factor(Direction, levels = direction_levels)) %>%
    ggplot(aes(species, estimate, color = Approach)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(ymin = LCI95, ymax = UCI95), width = 0, alpha = 0.8,
                  position = position_dodge(width = 0.7)) +
    geom_point(size = 1.4, position = position_dodge(width = 0.7)) +
    scale_color_manual(values = approach_cols, labels = approach_labs, drop = FALSE) +
    facet_grid(cols = vars(Direction), scales = "free_x", space = "free_x") +
    labs(title = title, x = NULL, y = expression(beta[Temp] ~ "on wing (standardized)"),
         color = NULL) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5),
          legend.position = "top", panel.spacing.x = unit(4, "pt"))
}

fig_estimates <- plot_grid(
  plot_estimates(nj$comparison,       "Nightjar (predictor: temperature)"),
  plot_estimates(weeks$comparison,    "Weeks 2020 (predictor: year)"),
  plot_estimates(atlantic$comparison, "Atlantic birds (predictor: temperature)"),
  ncol = 1, rel_heights = c(1, 1, 1), labels = "auto")

ggsave("Figures/Empirical_Btemp_comparison.png", fig_estimates,
       bg = "white", width = 15, height = 15, dpi = 200)
cat("\nSaved Figures/Empirical_Btemp_comparison.png\n")

## Report the Bergmann-direction breakdown per dataset.
cat("\nSpecies per Bergmann direction:\n")
print(map(datasets, \(d) d$direction %>% count(Direction) %>% mutate(study = d$study)) %>%
        list_rbind() %>% pivot_wider(names_from = Direction, values_from = n, values_fill = 0))


# Nightjar sample-size experiment ----
## Is the nightjar SEM's imprecision just small N? Simulate multivariate-normal
## data that reproduces each nightjar species' observed correlation structure among
## (Temp, Mass, Wing, Tail), scale N up, and refit the SAME SEM. If SE(beta_wing)
## shrinks toward a usable level with high convergence, sample size is the binding
## constraint; if the loading lambda stays large and unstable, the limit is
## structural (mass is a noisy anchor / wing and tail barely share a size factor).
if (run_sample_size_sim) {
nj_cor_list <- nj_df %>%
  group_split(Species) %>%
  set_names(map_chr(., \(d) unique(d$Species))) %>%
  map(\(d) cor(dplyr::select(d, B.Temp, log_mass, log_wing, log_tail)))

## One MVN draw of size N from correlation matrix R, fit the SEM, return wing stats.
sim_nj_fit <- function(R, N) {
  X <- MASS::mvrnorm(N, mu = rep(0, 4), Sigma = R)
  colnames(X) <- c("B.Temp", "log_mass", "log_wing", "log_tail")
  res <- fit_lavaan_sem(as_tibble(X),
    mass_name = "log_mass", append_names = c("log_wing", "log_tail"),
    temp_name = "B.Temp", labels = c("wing", "tail"))
  tibble(est = res$coef_sem_wing, se = res$se_sem_wing,
         lambda = res$lambda_sem_wing)
}

set.seed(42)
N_grid <- c(430, 1000, 2500, 5000, 10000)
n_reps <- 80

nj_power <- imap(nj_cor_list, \(R, sp) {
  map(N_grid, \(N) {
    reps <- map(seq_len(n_reps), \(i) sim_nj_fit(R, N)) %>% list_rbind()
    tibble(species = sp, N = N,
           convergence = mean(!is.na(reps$est)),      # share of reps that returned a usable fit
           mean_se     = mean(reps$se,     na.rm = TRUE),  # model-based SE of the wing coef
           sd_est      = sd(reps$est,      na.rm = TRUE),  # actual spread of the point estimate
           mean_est    = mean(reps$est,    na.rm = TRUE),
           mean_lambda = mean(reps$lambda, na.rm = TRUE))
  }) %>% list_rbind()
}) %>% list_rbind()

cat("\n===== NIGHTJAR SAMPLE-SIZE EXPERIMENT =====\n")
cat(sprintf("%d reps per N; observed per-species N ~ 267-443.\n", n_reps))
print(nj_power %>% mutate(across(c(convergence, mean_se, sd_est, mean_est, mean_lambda),
                                 \(x) round(x, 3))), n = Inf)

## Two panels: (a) precision of the wing coef vs N (does more data fix it?);
## (b) the estimated loading vs N (is the large lambda a structural feature that
## persists regardless of N?).
p_se <- ggplot(nj_power, aes(N, mean_se, color = species)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
  geom_line() + geom_point(size = 2) +
  scale_x_log10() +
  annotate("text", x = min(N_grid), y = 0.11, hjust = 0, size = 3,
           color = "grey40", label = "SE = 0.1 (usable)") +
  labs(x = "Sample size per species (log scale)",
       y = "Mean SE of wing coefficient", color = NULL,
       title = "(a) Precision improves with N")

p_lambda <- ggplot(nj_power, aes(N, mean_lambda, color = species)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_line() + geom_point(size = 2) +
  scale_x_log10() +
  annotate("text", x = min(N_grid), y = 1.1, hjust = 0, size = 3,
           color = "grey40", label = "lambda = 1 (= mass anchor)") +
  labs(x = "Sample size per species (log scale)",
       y = expression("Mean wing loading" ~ lambda), color = NULL,
       title = "(b) Loading is structural, ~flat in N")

fig_power <- plot_grid(p_se + theme(legend.position = "none"),
                       p_lambda, nrow = 1, rel_widths = c(1, 1.25))
ggsave("Figures/Nightjar_sample_size.png", fig_power,
       bg = "white", width = 12, height = 5, dpi = 200)
cat("Saved Figures/Nightjar_sample_size.png\n")
}  # end run_sample_size_sim


# Bayesian latent SEM (blavaan) ----
## Fits the same MIMIC model with blavaan and a prior on the appendage loadings
## (fit_blavaan_sem, in Key_causal_fns.R). Two deliverables:
##   (1) a nightjar prior-sensitivity strip showing how the wing effect and loading
##       move as the prior tightens from weak toward near-isometric — the knob for
##       choosing a prior;
##   (2) optionally (blavaan_full = TRUE) a full empirical refit of all species with
##       the chosen `blavaan_prior`, regenerating the three-panel comparison figure.
## Off by default (blavaan fits are slow). Set run_blavaan <- TRUE to enable, and
## tune blavaan_prior_grid / blavaan_prior.
if (run_blavaan) {
  library(blavaan)
  options(mc.cores = 4)

  ## Reference ratio / Ryding wing estimates per nightjar species (for the strip).
  nj_ref <- nj$comparison %>%
    filter(appendage == "wing", Approach %in% c("Ratio", "Ryding")) %>%
    dplyr::select(species, Approach, estimate)

  ## Prior-sensitivity: refit each nightjar species at every prior width.
  ## Cached — delete Derived/nj_blavaan_prior_sens.rds (or change the grid) to refit.
  sens_cache <- "Derived/nj_blavaan_prior_sens.rds"
  if (file.exists(sens_cache)) {
    nj_prior_sens <- readRDS(sens_cache)
  } else {
    nj_prior_sens <- map(blavaan_prior_grid, \(pr) {
      run_sem_list(nj_list, "B.Temp", c("wing", "tail"), c("log_wing", "log_tail"),
                   fit_fn = fit_blavaan_sem, loading_prior = pr) %>%
        mutate(prior = pr)
    }) %>% list_rbind() %>%
      mutate(prior_sd = as.numeric(str_extract(prior, "(?<=,)[0-9.]+(?=\\))")))
    saveRDS(nj_prior_sens, sens_cache)
  }

  cat("\n===== NIGHTJAR blavaan PRIOR SENSITIVITY (wing) =====\n")
  print(nj_prior_sens %>% filter(appendage == "wing") %>%
          arrange(species, desc(prior_sd)) %>%
          dplyr::select(species, prior, estimate, lci, uci, lambda), n = Inf)

  ## Strip: wing effect (+ HPD) vs prior width, with ratio / Ryding as reference lines.
  p_sens_est <- nj_prior_sens %>% filter(appendage == "wing") %>%
    ggplot(aes(prior_sd, estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_hline(data = nj_ref, aes(yintercept = estimate, color = Approach), linetype = "dotted") +
    geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, alpha = 0.7) +
    geom_line(color = "grey40") + geom_point(size = 2) +
    scale_color_manual(values = c(Ratio = "#E69F00", Ryding = "#0072B2"), name = "reference") +
    facet_wrap(~species, scales = "free_y") +
    labs(x = "loading-prior SD  (right = looser; left = tighter / more isometric)",
         y = expression(beta[Temp] ~ "on wing (SEM, 95% HPD)"),
         title = "Nightjar: direct wing effect vs. loading-prior width")

  p_sens_lam <- nj_prior_sens %>% filter(appendage == "wing") %>%
    ggplot(aes(prior_sd, lambda, color = species)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
    geom_line() + geom_point(size = 2) +
    labs(x = "loading-prior SD  (right = looser; left = tighter)",
         y = expression("Wing loading" ~ lambda), color = NULL,
         title = "Loadings shrink toward 1 as the prior tightens")

  fig_sens <- plot_grid(p_sens_est, p_sens_lam, ncol = 1, rel_heights = c(1.2, 1), labels = "auto")
  ggsave("Figures/Nightjar_blavaan_prior_sensitivity.png", fig_sens,
         bg = "white", width = 11, height = 9, dpi = 200)
  cat("Saved Figures/Nightjar_blavaan_prior_sensitivity.png\n")

  ## Full empirical refit with the chosen prior, regenerating the comparison figure.
  if (blavaan_full) {
    cat(sprintf("\nRefitting ALL species with blavaan prior %s (slow)...\n", blavaan_prior))
    ## The expensive blavaan SEM fits are cached per dataset+prior, so re-plotting
    ## (e.g. tweaking facets or adding metrics) is cheap. Delete the RDS to refit.
    prior_tag <- gsub("[^0-9a-z]+", "", blavaan_prior)
    build_blavaan_comparison <- function(sp_list, temp_name, labels, log_cols, ratio_cols,
                                         direction, cache_tag) {
      cache_f <- sprintf("Derived/blavaan_sem_%s_%s.rds", cache_tag, prior_tag)
      if (file.exists(cache_f)) {
        sem <- readRDS(cache_f)
      } else {
        sem <- run_sem_list(sp_list, temp_name, labels, log_cols,
                            fit_fn = fit_blavaan_sem, loading_prior = blavaan_prior,
                            n_chains = 2, burnin = 500, sample = 500)   # lighter chains for the sweep
        saveRDS(sem, cache_f)
      }
      simple <- run_simple_methods(sp_list, temp_name, labels, log_cols, ratio_cols)
      bind_rows(sem, simple) %>%
        mutate(LCI95 = coalesce(lci, estimate - 1.96 * std.error),
               UCI95 = coalesce(uci, estimate + 1.96 * std.error),
               Approach = factor(Approach, levels = approach_levels)) %>%
        left_join(direction, by = "species")
    }
    nj_cmp_b  <- build_blavaan_comparison(nj_list,       "B.Temp", c("wing","tail"),   c("log_wing","log_tail"),   c("wing_mass","tail_mass"),   nj$direction,       "nightjar")
    wk_cmp_b  <- build_blavaan_comparison(weeks_list,    "year",   c("wing","tarsus"), c("log_wing","log_tarsus"), c("wing_mass","tarsus_mass"), weeks$direction,    "weeks")
    atl_cmp_b <- build_blavaan_comparison(atlantic_list, "B.Tavg", c("wing","tarsus"), c("log_wing","log_tarsus"), c("wing_mass","tarsus_mass"), atlantic$direction, "atlantic")

    fig_estimates_b <- plot_grid(
      plot_estimates(nj_cmp_b,  sprintf("Nightjar (blavaan, %s)", blavaan_prior)),
      plot_estimates(wk_cmp_b,  sprintf("Weeks 2020 (blavaan, %s)", blavaan_prior)),
      plot_estimates(atl_cmp_b, sprintf("Atlantic birds (blavaan, %s)", blavaan_prior)),
      ncol = 1, labels = "auto")
    ggsave("Figures/Empirical_Btemp_comparison_blavaan.png", fig_estimates_b,
           bg = "white", width = 15, height = 15, dpi = 200)
    cat("Saved Figures/Empirical_Btemp_comparison_blavaan.png\n")

    ## Persist the estimates and print a precision summary per dataset.
    blavaan_all <- bind_rows(
      nj_cmp_b  %>% mutate(study = "Nightjar"),
      wk_cmp_b  %>% mutate(study = "Weeks (2020)"),
      atl_cmp_b %>% mutate(study = "Atlantic birds")) %>%
      mutate(prior = blavaan_prior)
    write_csv(blavaan_all, "Derived/blavaan_empirical_estimates.csv")
    cat("Saved Derived/blavaan_empirical_estimates.csv\n")

    cat("\n===== blavaan vs. Ratio / Ryding: wing precision by dataset =====\n")
    print(blavaan_all %>% filter(appendage == "wing", !is.na(estimate)) %>%
            group_by(study, Approach) %>%
            summarise(n = n(), median_est = median(estimate),
                      median_ci_width = median(UCI95 - LCI95), .groups = "drop"))

    ## Direction-grouped figures pooling all three studies (blavaan SEM estimates).
    ## Species labels are coloured by study; within a panel species are grouped by
    ## study then ordered by mean estimate. Fig 1: Bergmann's + Inverse Bergmann's.
    ## Fig 2: the two Mixed classes side by side, panel width proportional to N species.
    library(patchwork)
    dd <- blavaan_all %>%
      filter(appendage == "wing", std.error < 0.5) %>%
      mutate(Approach = factor(Approach, levels = approach_levels),
             study    = factor(study, levels = names(study_colors)))

    direction_panel <- function(dir) {
      di   <- dd %>% filter(Direction == dir)
      meta <- di %>% group_by(species, study) %>% summarise(m = mean(estimate), .groups = "drop") %>%
        arrange(study, m)
      di %>% mutate(species = factor(species, levels = meta$species)) %>%
        ggplot(aes(species, estimate, color = Approach)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
        geom_errorbar(aes(ymin = LCI95, ymax = UCI95), width = 0, alpha = 0.8,
                      position = position_dodge(0.7)) +
        geom_point(size = 1.3, position = position_dodge(0.7)) +
        scale_color_manual(values = approach_cols, labels = approach_labs) +
        labs(title = dir, x = NULL, y = expression(beta[Temp] ~ "on wing (std.)"), color = NULL) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5,
                                         color = study_colors[as.character(meta$study)]),
              legend.position = "top")
    }
    study_cap <- "Species labels coloured by study: Nightjar (red), Weeks 2020 (blue), Atlantic birds (green)"

    fig_berg <- (direction_panel("Bergmann's") / direction_panel("Inverse Bergmann's")) +
      plot_layout(guides = "collect") + plot_annotation(tag_levels = "a", caption = study_cap) &
      theme(legend.position = "top", plot.caption = element_text(hjust = 0, size = 9))
    ggsave("Figures/Empirical_bergmann_directions.png", fig_berg,
           bg = "white", width = 13, height = 11, dpi = 200)
    cat("Saved Figures/Empirical_bergmann_directions.png\n")

    n_w <- dd %>% filter(Direction == "Mixed - Wingier") %>% distinct(species) %>% nrow()
    n_f <- dd %>% filter(Direction == "Mixed - Fatter")  %>% distinct(species) %>% nrow()
    fig_mixed <- (direction_panel("Mixed - Wingier") + direction_panel("Mixed - Fatter") +
                    plot_layout(widths = c(max(n_w, 1), max(n_f, 1)))) +
      plot_layout(guides = "collect") + plot_annotation(tag_levels = "a", caption = study_cap) &
      theme(legend.position = "top", plot.caption = element_text(hjust = 0, size = 9))
    ggsave("Figures/Empirical_mixed_directions.png", fig_mixed,
           bg = "white", width = 11, height = 5.5, dpi = 200)
    cat("Saved Figures/Empirical_mixed_directions.png\n")
  }
}  # end run_blavaan


# Identifiability diagnostics ----
## Why does the latent SEM misbehave on real data? These panels separate the TWO
## distinct failure modes and identify which drives what:
##   (1) wide credible intervals -> IMPRECISION (small N + mass a labile size
##       indicator); shrinks ~1/sqrt(N) and is fixable with more data.
##   (2) rank-order anomaly (SEM estimate above the ratio for Bergmann's species)
##       -> BIAS from misspecification: wing and the 2nd appendage share variance
##       beyond mass-defined size (a second "limb" factor) that a single mass-anchored
##       factor cannot represent. This bias does NOT shrink with N.
##
## The simulation establishing (1) vs (2) lives in 00_Simulation_data.R (its figure is
## Figures/SEM_bias_vs_variance.png, drawn by 01_Simulation_plots.R). What follows is
## the EMPIRICAL fingerprint of the same two modes.
##
## Critical caveat: with only 3 indicators (mass + 2 appendages) the single-factor
## model is just-identified (df = 0) and cannot be tested — noisy-mass and two-factor
## structures are observationally equivalent (both inflate the loading and the
## wing-2nd-appendage partial correlation). A 4th indicator separates them (Part C).

if (run_diagnostics && file.exists("Derived/blavaan_empirical_estimates.csv")) {
library(patchwork)

# A: empirical fingerprint of the two failure modes ----

## Partial correlation cor(wing, 2nd appendage | mass) per species: the single-factor
## strain signal. High values mean the appendages share variance mass cannot explain.
pcor_species <- function(sp_list, app2, study) {
  imap(sp_list, \(d, sp) {
    r   <- cor(d[c("log_wing", app2, "log_mass")], use = "complete.obs")
    rw  <- r[1, 3]; rt <- r[2, 3]; rwt <- r[1, 2]
    tibble(species = sp, study = study,
           pcor = (rwt - rw * rt) / sqrt((1 - rw^2) * (1 - rt^2)))
  }) %>% list_rbind()
}
emp_pcor <- bind_rows(
  pcor_species(nj_list,       "log_tail",   "Nightjar"),
  pcor_species(weeks_list,    "log_tarsus", "Weeks (2020)"),
  pcor_species(atlantic_list, "log_tarsus", "Atlantic birds"))

blav <- read_csv("Derived/blavaan_empirical_estimates.csv", show_col_types = FALSE) %>%
  filter(appendage == "wing")
sem_stats <- blav %>% filter(Approach == "SEM") %>%
  transmute(species, lambda, ci_width = UCI95 - LCI95, Direction)
n_tbl <- bind_rows(
  imap(nj_list,       \(d, sp) tibble(species = sp, n = nrow(d))) %>% list_rbind(),
  imap(weeks_list,    \(d, sp) tibble(species = sp, n = nrow(d))) %>% list_rbind(),
  imap(atlantic_list, \(d, sp) tibble(species = sp, n = nrow(d))) %>% list_rbind())

emp <- emp_pcor %>%
  left_join(sem_stats, by = "species") %>%
  left_join(n_tbl, by = "species") %>%
  mutate(study = factor(study, levels = names(study_colors)))

## Panel A: the 2nd-factor strain is nightjar-specific; for Weeks/Atlantic a single
## mass-anchored factor is adequate (mass mediates the wing-tarsus covariance).
pA <- emp %>%
  ggplot(aes(study, pcor, color = study)) +
  geom_boxplot(outlier.shape = NA, color = "grey70", width = 0.5) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7, size = 1.6) +
  scale_color_manual(values = study_colors, guide = "none") +
  labs(x = NULL, y = "partial cor(wing, 2nd appendage | mass)",
       title = "(a) 2nd-factor strain is nightjar-specific") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

## Panel B: the wide SEM intervals are imprecision — they shrink with per-species N.
pB <- emp %>%
  ggplot(aes(n, ci_width, color = study)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "grey30", linewidth = 0.7) +
  scale_x_log10() +
  scale_color_manual(values = study_colors, name = NULL) +
  labs(x = "sample size per species (log scale)", y = "SEM 95% interval width",
       title = "(b) Wide intervals shrink with sample size")

ggsave("Figures/SEM_empirical_fingerprint.png", (pA | pB) + plot_layout(widths = c(1, 1.3)) &
         theme(legend.position = "top"), bg = "white", width = 12, height = 5, dpi = 200)
cat("Saved Figures/SEM_empirical_fingerprint.png\n")

cat("\nMedian partial cor(wing, 2nd appendage | mass) and lambda by study:\n")
print(emp %>% group_by(study) %>%
        summarise(median_pcor = round(median(pcor), 2),
                  median_lambda = round(median(lambda, na.rm = TRUE), 2), n = n()))
cat(sprintf("\nCor(SEM interval width, log N) across species: %.2f\n",
            cor(emp$ci_width, log(emp$n), use = "complete.obs")))

## Is the SEM systematically biased vs the ratio for Bergmann's species?
berg <- blav %>% filter(Direction == "Bergmann's") %>%
  dplyr::select(species, Approach, estimate) %>%
  pivot_wider(names_from = Approach, values_from = estimate)
cat(sprintf("Bergmann's species (n=%d): SEM above ratio %.0f%%, below Ryding %.0f%%, within bracket %.0f%%\n",
            nrow(berg), 100*mean(berg$SEM > berg$Ratio, na.rm = TRUE),
            100*mean(berg$SEM < berg$Ryding, na.rm = TRUE),
            100*mean(berg$SEM <= berg$Ratio & berg$SEM >= berg$Ryding, na.rm = TRUE)))


# B: does the SEM leave the [Ryding, ratio] bracket when strain is high? ----
## Hypothesis: a species' SEM estimate falls OUTSIDE the bracket set by the ratio and
## mass-as-covariate estimates when its 2nd-factor strain is high. Tested for the
## unregularized (lavaan) SEM — where any bias is unmasked — and the regularized
## (blavaan) SEM. Result: directionally supported but weak (the prior mutes it, and
## partial-correlation strain is a confounded proxy for a 2nd factor).
raw_sem_wing <- function(sp_list, temp, app2) {
  imap(sp_list, \(d, sp) {
    z <- d %>% mutate(across(c(log_mass, log_wing, all_of(app2), all_of(temp)),
                             \(x) as.numeric(scale(x))))
    m <- sprintf("Sz =~ 1*log_mass + log_wing + %s\n Sz ~ %s\n log_wing ~ %s\n %s ~ %s",
                 app2, temp, temp, app2, temp)
    f <- tryCatch(suppressWarnings(sem(m, z)), error = function(e) NULL)
    est <- if (is.null(f) || !lavInspect(f, "converged")) NA_real_ else {
      pe <- parameterEstimates(f); pe$est[pe$lhs == "log_wing" & pe$op == "~" & pe$rhs == temp] }
    tibble(species = sp, sem_raw = est)
  }) %>% list_rbind()
}
sem_raw <- bind_rows(
  raw_sem_wing(nj_list,       "B.Temp", "log_tail"),
  raw_sem_wing(weeks_list,    "year",   "log_tarsus"),
  raw_sem_wing(atlantic_list, "B.Tavg", "log_tarsus"))

bracket <- blav %>% dplyr::select(species, Approach, estimate) %>%
  pivot_wider(names_from = Approach, values_from = estimate) %>%
  dplyr::select(species, Ratio, Ryding, SEM_blav = SEM) %>%
  left_join(sem_raw, by = "species") %>%
  left_join(emp %>% dplyr::select(species, study, pcor), by = "species") %>%
  mutate(lo = pmin(Ratio, Ryding), hi = pmax(Ratio, Ryding),
         `Unregularized (lavaan)` = sem_raw  >= lo & sem_raw  <= hi,
         `Regularized (blavaan)`  = SEM_blav >= lo & SEM_blav <= hi) %>%
  pivot_longer(c(`Unregularized (lavaan)`, `Regularized (blavaan)`),
               names_to = "sem_type", values_to = "in_bracket") %>%
  filter(!is.na(in_bracket)) %>%
  mutate(sem_type = factor(sem_type, levels = c("Unregularized (lavaan)", "Regularized (blavaan)")),
         in_bracket = factor(in_bracket, levels = c(FALSE, TRUE),
                             labels = c("outside\nbracket", "inside\nbracket")))

ggsave("Figures/SEM_bracket_test.png",
  bracket %>%
    ggplot(aes(in_bracket, pcor)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, color = "grey55") +
    geom_jitter(aes(color = study), width = 0.12, height = 0, alpha = 0.7, size = 1.5) +
    scale_color_manual(values = study_colors, name = NULL) +
    facet_wrap(~sem_type) +
    labs(x = "Is the SEM estimate inside the [mass-as-covariate, ratio] bracket?",
         y = "2nd-factor strain\n(partial cor wing–appendage | mass)",
         title = "SEM outside the bracket tends to have higher 2nd-factor strain (weakly)") +
    theme(legend.position = "top"),
  bg = "white", width = 10, height = 5, dpi = 200)
cat("\nSaved Figures/SEM_bracket_test.png\n")
for (st in levels(bracket$sem_type)) {
  s <- bracket %>% filter(sem_type == st)
  cat(sprintf("  %s: median strain outside=%.2f inside=%.2f  (Wilcoxon p=%.3f)\n", st,
      median(s$pcor[s$in_bracket == "outside\nbracket"]),
      median(s$pcor[s$in_bracket == "inside\nbracket"]),
      wilcox.test(pcor ~ in_bracket, data = s)$p.value))
}


# C: 4-indicator SEM for the Atlantic birds (mass + wing + tarsus + tail) ----
## Adding a 4th indicator over-identifies the size factor (df > 0), so the single-factor
## assumption becomes TESTABLE and the loadings are better anchored. Question: does a
## second factor show up, and does the wing estimate move?
atl4 <- read_csv(atlantic_csv, show_col_types = FALSE) %>% janitor::clean_names() %>%
  rename(mass = body_mass_g, B.Tavg = annual_mean_temperature) %>%
  mutate(species_ = str_replace(binomial, " ", "_"),
         log_wing   = log(coalesce(wing_length_mm, wing_length_left_mm, wing_length_right_mm)),
         log_tarsus = log(coalesce(tarsus_length_mm, tarsus_length_right_mm, tarsus_length_left_mm)),
         log_tail   = log(coalesce(tail_length_mm, tail_length_right_mm, tail_length_left_mm)),
         log_mass   = log(mass)) %>%
  filter(year > 1990, !is.na(B.Tavg)) %>%
  drop_na(log_mass, log_wing, log_tarsus, log_tail) %>%
  add_count(species_) %>% filter(n >= 150) %>%
  group_split(species_) %>% set_names(map_chr(., \(d) unique(d$species_))) %>%
  map(\(d) mutate(d, across(c(log_mass, log_wing, log_tarsus, log_tail, B.Tavg),
                            \(x) as.numeric(scale(x)))))

wing_row <- function(f, temp = "B.Tavg") {
  if (is.null(f) || !lavInspect(f, "converged")) return(c(est = NA, se = NA, cfi = NA, rmsea = NA))
  pe <- parameterEstimates(f) %>% filter(lhs == "log_wing", op == "~", rhs == temp)
  c(est = pe$est, se = pe$se,
    cfi = as.numeric(fitMeasures(f, "cfi")), rmsea = as.numeric(fitMeasures(f, "rmsea")))
}
atl4_res <- imap(atl4, \(d, sp) {
  m3 <- 'Sz =~ 1*log_mass + log_wing + log_tarsus
         Sz ~ B.Tavg
         log_wing ~ B.Tavg
         log_tarsus ~ B.Tavg'
  m4 <- 'Sz =~ 1*log_mass + log_wing + log_tarsus + log_tail
         Sz ~ B.Tavg
         log_wing ~ B.Tavg
         log_tarsus ~ B.Tavg
         log_tail ~ B.Tavg'
  v3 <- wing_row(tryCatch(suppressWarnings(sem(m3, d)), error = function(e) NULL))
  v4 <- wing_row(tryCatch(suppressWarnings(sem(m4, d)), error = function(e) NULL))
  tibble(species = sp, est3 = v3["est"], se3 = v3["se"],
         est4 = v4["est"], se4 = v4["se"], cfi4 = v4["cfi"], rmsea4 = v4["rmsea"])
}) %>% list_rbind()

cat(sprintf("\nAtlantic 4-indicator: %d/%d converged; single-factor fit median CFI=%.2f RMSEA=%.2f; %.0f%% good fit\n",
            sum(!is.na(atl4_res$est4)), nrow(atl4_res),
            median(atl4_res$cfi4, na.rm = TRUE), median(atl4_res$rmsea4, na.rm = TRUE),
            100 * mean(atl4_res$cfi4 > .95 & atl4_res$rmsea4 < .08, na.rm = TRUE)))

pE <- atl4_res %>% filter(!is.na(est4), !is.na(est3)) %>%
  ggplot(aes(est3, est4)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(color = "#4DAF4A", alpha = 0.8, size = 1.8) +
  coord_cartesian(xlim = c(-0.5, 1), ylim = c(-0.5, 1)) +   # 1 species blows up (est4~25); zoom to the bulk
  labs(x = "wing~Temp, 3-indicator (mass+wing+tarsus)", y = "wing~Temp, 4-indicator (+ tail)",
       title = "(a) Estimate barely moves with a 4th indicator")
## With 4 indicators the single-factor model is testable: good fit => no 2nd factor.
pF <- atl4_res %>% filter(!is.na(cfi4)) %>%
  ggplot(aes(rmsea4, cfi4)) +
  annotate("rect", xmin = -Inf, xmax = 0.08, ymin = 0.95, ymax = Inf, fill = "#4DAF4A", alpha = 0.12) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0.08, linetype = "dotted", color = "grey50") +
  geom_point(color = "#4DAF4A", alpha = 0.8, size = 1.8) +
  labs(x = "RMSEA (4-indicator single factor)", y = "CFI",
       title = "(b) Single factor fits most species (shaded = good fit)")
ggsave("Figures/Atlantic_4indicator.png", (pE | pF), bg = "white", width = 11, height = 5, dpi = 200)
cat("Saved Figures/Atlantic_4indicator.png\n")
cat(sprintf("median |est4-est3|=%.3f; median SE 3-ind=%.3f, 4-ind=%.3f (%.0f%% tighter with 4 indicators)\n",
            median(abs(atl4_res$est4 - atl4_res$est3), na.rm = TRUE),
            median(atl4_res$se3, na.rm = TRUE), median(atl4_res$se4, na.rm = TRUE),
            100 * mean(atl4_res$se4 < atl4_res$se3, na.rm = TRUE)))

}  # end run_diagnostics
