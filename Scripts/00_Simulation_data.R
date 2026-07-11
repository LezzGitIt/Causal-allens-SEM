# 00_Simulation_data.R
# Generates all simulated data, fits every method to it, and exports the results.
# Produces no figures: everything plotted lives in 01_Simulation_plots.R, which reads
# only the files exported here. Needs no external data, so it runs from a clean clone.
#
# Dependencies: Key_causal_fns.R
# Run order: this script -> 01_Simulation_plots.R -> render Scripts/Causal_allens.qmd
#
# Outputs:
#   Derived/Csv/results_tbl.csv     per-species estimates from all 5 estimators (no ME)
#   Derived/Csv/results_me_tbl.csv  per-species estimates, oracle vs M1, under mass ME
#   Derived/Csv/parms_mat.csv       full parameter grid + induced correlations + realism flag
#   Derived/Csv/allom_tbl.csv       per-species OLS and SMA allometric slopes
#   Derived/Csv/misspec_sim.csv     bias / SE vs N for the misspecification simulation
#   Derived/sim_morphology.rds      raw-scale morphology of every realistic dataset
#   Derived/sim_example.rds         one example dataset (+ its ME variant) and its fitted lines


library(tidyverse)
library(smatr)
library(lavaan)
library(MASS)
library(janitor)

source("Scripts/Key_causal_fns.R")

dir.create("Derived/Csv", recursive = TRUE, showWarnings = FALSE)

# 1. Parameter grid -------------------------------------------------------
# Solver functions (solve_lambda_w, solve_sigma_w, implied_allometry) come from Key_causal_fns.R.

### Biologically realistic allometric targets.
target_sma_vec <- c(0.3, 0.45, 0.6)

# Max is 0.65 (not 0.70): r_WM = target_r by construction, so 0.70 hits the filter ceiling and fails ~50% of empirical datasets by chance.
target_r_vec <- c(0.30, 0.50, 0.65)

# B_temp_w_vec: kept small so r(Wing,Temp) stays negative for all but the most extreme Allen's + weak-Bergmann combinations. Max before r_WT flips sign ≈ 0.091 for beta_size = -0.15, so 0.07 leaves a small safety margin.
B_temp_w_vec <- c(0.00, 0.02, 0.05, 0.07)

# beta_size: chosen so r(Mass,Temp) ≈ {-0.30, -0.50, -0.65} with sigma_temp = 0.3. Values >= |0.46| push r_MT below -0.70 (filter lower bound) and are excluded.
parms_mat <- expand_grid(
  target_sma = target_sma_vec,
  target_r   = target_r_vec,
  B_temp_w   = B_temp_w_vec,
  beta_size  = c(-0.15, -0.27, -0.42)
) %>%
  mutate(
    target_ols = target_sma * target_r,

    lambda_wing = solve_lambda_w(
      target_r = target_r, target_sma = target_sma,
      beta_size = beta_size, B_temp = B_temp_w
    ),
    sigma_wing = solve_sigma_w(
      target_sma = target_sma, lambda_w = lambda_wing,
      beta_size = beta_size, B_temp = B_temp_w
    ),

    ## Tarsus is biologically distinct from wing: tighter coupling to latent size, slightly lower allometric slope, and an equal direct Allen's response.
    B_temp_t     = B_temp_w,
    target_r_t   = pmin(target_r + 0.08, 0.62),
    target_sma_t = target_sma - 0.05,

    lambda_tarsus = solve_lambda_w(
      target_r = target_r_t, target_sma = target_sma_t,
      beta_size = beta_size, B_temp = B_temp_t
    ),
    sigma_tarsus = solve_sigma_w(
      target_sma = target_sma_t, lambda_w = lambda_tarsus,
      beta_size = beta_size, B_temp = B_temp_t
    ),

    sigma_mass = 0.10,

    Allen = factor(
      case_when(
        B_temp_w == 0.00 ~ "No support",
        B_temp_w == 0.02 ~ "Weak",
        B_temp_w == 0.05 ~ "Moderate",
        B_temp_w == 0.07 ~ "Strong"
      ),
      levels = c("Strong", "Moderate", "Weak", "No support")
    ),
    Bergmann = factor(
      case_when(
        beta_size == -0.15 ~ "Weak",
        beta_size == -0.27 ~ "Moderate",
        TRUE               ~ "Strong"
      ),
      levels = c("Weak", "Moderate", "Strong")
    )
  ) %>%
  filter(!is.na(sigma_wing), !is.na(sigma_tarsus))  # drop impossible combinations

## Verify implied allometries before simulating: rho should match target_r and beta_sma should match target_sma for every row.
allometry_check <- parms_mat %>%
  distinct(target_sma, target_r, B_temp_w, beta_size, lambda_wing, sigma_wing) %>%
  bind_cols(
    implied_allometry(
      lambda_w  = .$lambda_wing, sigma_w = .$sigma_wing,
      beta_size = .$beta_size,   B_temp  = .$B_temp_w
    )
  )
stopifnot(max(abs(allometry_check$rho - allometry_check$target_r)) < 1e-8)

## Screen the grid for empirically realistic appendage-mass and temperature-morphology correlations, using targets taken from published avian studies.
set.seed(42)
corr_tbl <- pmap(
  parms_mat %>% dplyr::select(lambda_wing, lambda_tarsus, B_temp_w, B_temp_t,
                              beta_size, sigma_wing, sigma_tarsus, sigma_mass),
  \(...) simulate_data(..., N = 1500)$corr
) %>%
  map(\(m) tibble(
    r_WM = m["Append",  "Mass"],
    r_WT = m["Append",  "Temp_inc"],
    r_MT = m["Mass",    "Temp_inc"],
    r_TM = m["Tarsus",  "Mass"],
    r_TT = m["Tarsus",  "Temp_inc"]
  )) %>%
  list_rbind()

parms_mat2 <- bind_cols(parms_mat, corr_tbl) %>%
  mutate(realistic = r_WM >= 0.25 & r_WM <= 0.70 &
                     r_WT >= -0.70 & r_WT <= 0.00 &
                     r_MT >= -0.70 & r_MT <= 0.00 &
                     r_TM >= 0.25  & r_TM <= 0.70 &
                     r_TT >= -0.70 & r_TT <= 0.00)

parms_realistic <- parms_mat2 %>% filter(realistic)
cat(sprintf("\n%d of %d parameter combinations retained as realistic.\n",
            nrow(parms_realistic), nrow(parms_mat2)))

# 2. Generate simulation datasets -----------------------------------------

set.seed(123)
df_list <- pmap(
  parms_realistic %>%
    dplyr::select(lambda_wing, lambda_tarsus, B_temp_w, B_temp_t, beta_size,
                  sigma_wing, sigma_tarsus, sigma_mass),
  \(...) gen_causal_data(..., N = 1500)
)

## Measurement-error datasets: sigma_me = 0.10 on the log scale (~10% CV at 50 g).
set.seed(789)
df_me_list <- pmap(
  parms_realistic %>%
    dplyr::select(lambda_wing, lambda_tarsus, B_temp_w, B_temp_t, beta_size,
                  sigma_wing, sigma_tarsus, sigma_mass),
  \(...) gen_causal_data(..., N = 1000, sigma_me = 0.10)
)

## Raw-scale morphology of every dataset, for the trait-distribution figure.
sim_morphology <- df_list %>%
  bind_rows(.id = "Species") %>%
  dplyr::select(Species, Mass_g = Mass, Appendage1_mm = Append, Appendage2_mm = Tarsus) %>%
  pivot_longer(-Species, names_to = "Measurement", values_to = "Size")
saveRDS(sim_morphology, "Derived/sim_morphology.rds", compress = "xz")

# 3. Example dataset ------------------------------------------------------
## One representative species (moderate Allen's, moderate Bergmann's) used for the
## data-inspection and allometric-line figures.

ex_idx   <- which(parms_realistic$Allen == "Moderate" & parms_realistic$beta_size == -0.27)
ex_idx   <- ex_idx[ceiling(length(ex_idx) / 2)]
df_ex    <- df_list[[ex_idx]]
df_ex_me <- df_me_list[[ex_idx]]

## OLS and SMA lines for the example, so the plotting script only draws them.
ex_ols_w <- lm(Append_log ~ Mass_log, data = df_ex)
ex_sma_w <- smatr::sma(Append_log ~ Mass_log, data = df_ex)
ex_ols_t <- lm(Tarsus_log ~ Mass_log, data = df_ex)
ex_sma_t <- smatr::sma(Tarsus_log ~ Mass_log, data = df_ex)

line_tbl <- tibble(
  Appendage = rep(c("App1 (wing)", "App2 (tarsus)"), each = 2),
  Method    = rep(c("OLS", "SMA"), 2),
  intercept = c(coef(ex_ols_w)[1], ex_sma_w$coef[[1]]["elevation", "coef(SMA)"],
                coef(ex_ols_t)[1], ex_sma_t$coef[[1]]["elevation", "coef(SMA)"]),
  slope     = c(coef(ex_ols_w)[2], ex_sma_w$coef[[1]]["slope", "coef(SMA)"],
                coef(ex_ols_t)[2], ex_sma_t$coef[[1]]["slope", "coef(SMA)"])
)

saveRDS(list(df_ex = df_ex, df_ex_me = df_ex_me, line_tbl = line_tbl,
             B_temp_w = parms_realistic$B_temp_w[ex_idx],
             beta_size = parms_realistic$beta_size[ex_idx]),
        "Derived/sim_example.rds", compress = "xz")

## 3b) Collider illustration (Figures/Mass_as_collider.png): NOT REGENERATED HERE.
## The archived PNG's source was lost, and its caption describes a scenario that cannot
## arise under the latent-size DAG this paper simulates: with no direct effect, the
## marginal and mass-adjusted temperature slopes both carry the sign of
## (lambda * beta_size), so they cannot have opposite signs. The caption asserts a
## positive total effect alongside a spurious negative mass-adjusted coefficient.
## Pending a decision on the intended data-generating process, the figure is committed
## as a static asset. See README "Known gaps".

# 4. Allometric scaling: OLS vs SMA ---------------------------------------
## SMA assumes error in both variables, OLS only in y. Because mass carries error, SMA
## slopes should exceed OLS slopes (OLS slope = r * SMA slope, with r < 1).

allom_raw <- map(df_list, \(df) {
  ols_w <- lm(Append_log ~ Mass_log, data = df)
  sma_w <- smatr::sma(Append_log ~ Mass_log, data = df)
  ols_t <- lm(Tarsus_log ~ Mass_log, data = df)
  sma_t <- smatr::sma(Tarsus_log ~ Mass_log, data = df)

  tibble(
    ols_slope_wing   = coef(ols_w)["Mass_log"],
    sma_slope_wing   = sma_w$coef[[1]]["slope", "coef(SMA)"],
    ols_slope_tarsus = coef(ols_t)["Mass_log"],
    sma_slope_tarsus = sma_t$coef[[1]]["slope", "coef(SMA)"],
    r_wing_mass      = cor(df$Append_log, df$Mass_log),
    r_tarsus_mass    = cor(df$Tarsus_log, df$Mass_log)
  )
}) %>% list_rbind()

allom_tbl <- bind_cols(parms_realistic, allom_raw)

# 5. Run all comparison methods -------------------------------------------

cat("Running classical methods + SEM on", length(df_list), "datasets...\n")
results_raw <- map(df_list, \(df) {
  df_s <- df %>% mutate(across(all_of(c("Mass_log", "Append_log", "Tarsus_log", "Temp_inc")),
                               \(x) as.numeric(scale(x))))
  bind_cols(
    tibble(sd_log_wing   = sd(df$Append_log),
           sd_log_tarsus = sd(df$Tarsus_log),
           sd_temp       = sd(df$Temp_inc)),
    apply_methods(df, "Append_log", "Append") %>% rename_with(\(x) paste0(x, "_wing")),
    apply_methods(df, "Tarsus_log", "Tarsus") %>% rename_with(\(x) paste0(x, "_tarsus")),
    fit_lavaan_sem(df_s)
  )
}) %>% list_rbind()

results_tbl <- bind_cols(parms_realistic, results_raw)

## Each estimator appears in an oracle (noiseless latent mass) and an M1 (one noisy
## weighing) variant, so the within-pair gap is the cost of measurement error and the
## between-pair gap is the bias intrinsic to the estimator.
cat("Running methods with measurement error on", length(df_me_list), "datasets...\n")
results_me_raw <- map(df_me_list, \(df) {
  df_m1    <- df %>% mutate(Mass_log = M1_log)
  sem_cols <- c("Mass_log", "Append_log", "Tarsus_log", "Temp_inc")
  df_s     <- df    %>% mutate(across(all_of(sem_cols), \(x) as.numeric(scale(x))))
  df_m1_s  <- df_m1 %>% mutate(across(all_of(sem_cols), \(x) as.numeric(scale(x))))

  bind_cols(
    tibble(sd_log_wing   = sd(df$Append_log),
           sd_log_tarsus = sd(df$Tarsus_log),
           sd_temp       = sd(df$Temp_inc)),
    ratio_oracle(df),
    ratio_m1(df),
    sli_est_oracle(df),
    sli_est_m1(df),
    ryding_oracle(df),
    ryding_m1(df),
    fit_lavaan_sem(df_s)    %>% rename(coef_sem_oracle_wing   = coef_sem_wing,
                                       coef_sem_oracle_tarsus = coef_sem_tarsus,
                                       se_sem_oracle_wing     = se_sem_wing,
                                       se_sem_oracle_tarsus   = se_sem_tarsus),
    fit_lavaan_sem(df_m1_s) %>% rename(coef_sem_m1_wing       = coef_sem_wing,
                                       coef_sem_m1_tarsus     = coef_sem_tarsus,
                                       se_sem_m1_wing         = se_sem_wing,
                                       se_sem_m1_tarsus       = se_sem_tarsus),
    fit_lavaan_sem_3mass(df)
  )
}) %>% list_rbind()

results_me_tbl <- bind_cols(parms_realistic, results_me_raw)

# 6. Misspecification simulation (Supporting Information) ------------------
## Crosses three data-generating scenarios with sample size. Bias is the single-factor
## SEM's wing coefficient minus an oracle regression on the true latent Size, so a
## correctly specified model has bias 0 at every N. Only the second-factor scenario
## should show bias that does not shrink as N grows.

scenarios <- tribble(
  ~scenario,             ~gamma, ~sigma_mass,
  "Well-specified",       0.0,   0.5,   # mass a good size proxy, one factor
  "Noisy mass",           0.0,   1.6,   # mass labile -> imprecision only
  "2nd (limb) factor",    1.0,   0.5)   # appendages share extra variance -> bias

## Population partial cor(wing, tarsus | mass) each scenario induces. Noisy mass and a
## second factor both inflate it, which is exactly why 3 indicators cannot separate them.
set.seed(7)
scen_pcor <- scenarios %>% rowwise() %>%
  mutate(pcor = pcor_wt_given_m(gen_sem_data(50000, gamma = gamma, sigma_mass = sigma_mass))) %>%
  ungroup()
cat("\nPopulation partial cor(wing, tarsus | mass) by scenario:\n"); print(scen_pcor)

set.seed(1)
N_grid <- c(150, 300, 600, 1200, 3000, 10000)
reps   <- 60
misspec_sim <- expand_grid(scenarios, N = N_grid) %>%
  rowwise() %>%
  mutate(fits = list(map(seq_len(reps),
           \(i) fit_sf_vs_oracle(gen_sem_data(N, gamma = gamma, sigma_mass = sigma_mass))) %>%
           list_rbind())) %>%
  mutate(Bias         = mean(fits$est - fits$oracle, na.rm = TRUE),
         `SE(B_temp)` = mean(fits$se,     na.rm = TRUE),
         `SE(lambda)` = mean(fits$se_lam, na.rm = TRUE)) %>%
  dplyr::select(-fits) %>% ungroup() %>%
  left_join(scen_pcor %>% dplyr::select(scenario, pcor), by = "scenario")

cat("\nMisspecification simulation summary:\n")
print(misspec_sim %>% mutate(across(where(is.numeric), \(x) round(x, 3))), n = Inf)

# 7. Export ---------------------------------------------------------------

write_csv(results_tbl,    "Derived/Csv/results_tbl.csv")
write_csv(results_me_tbl, "Derived/Csv/results_me_tbl.csv")
write_csv(parms_mat2,     "Derived/Csv/parms_mat.csv")
write_csv(allom_tbl,      "Derived/Csv/allom_tbl.csv")
write_csv(misspec_sim,    "Derived/Csv/misspec_sim.csv")
cat("\nSaved Derived/Csv/{results_tbl,results_me_tbl,parms_mat,allom_tbl,misspec_sim}.csv\n")
cat("Saved Derived/{sim_morphology,sim_example}.rds\n")
