## 03_Hierarchical_SEM.R ----
## Hierarchical latent-size SEM across all 95 species and 3 studies.
## Fits Scripts/latent_size_hier.stan: App1 loading pooled globally; Bergmann slope
## and direct effects pooled within study; species partial-pooled. Estimand:
## per-species d_wing and per-study mu_dw.
##
## NEEDS THE RAW DATA (sourced via 02_Empirical_SEM.R). Its output,
## Derived/hier_sem_estimates.rds, is committed so 04_Empirical_plots.R can rebuild
## the empirical figures without it.
##
## Run order: 02_Empirical_SEM.R -> this script -> 04_Empirical_plots.R

library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

## Per-species z-scored lists + Bergmann direction. The slow sections of 02 are skipped:
## we only need the per-species data lists it builds.
run_sample_size_sim <- FALSE
run_blavaan         <- FALSE
run_diagnostics     <- FALSE
source("Scripts/02_Empirical_SEM.R")

# Build sufficient statistics ----
## Each species contributes n, sum T^2, sum T*y, sum y y' with y = (mass, wing, app2).
specs <- list(
  list(lst = nj_list,       temp = "B.Temp", app2 = "log_tail",   study = 1L),
  list(lst = weeks_list,    temp = "year",   app2 = "log_tarsus", study = 2L),
  list(lst = atlantic_list, temp = "B.Tavg", app2 = "log_tarsus", study = 3L))
study_names <- c("Nightjar", "Weeks (2020)", "Atlantic birds")

ss <- imap(specs, \(sp, gi) {
  imap(sp$lst, \(d, name) {
    Y  <- as.matrix(d[, c("log_mass", "log_wing", sp$app2)])
    T  <- d[[sp$temp]]
    ok <- stats::complete.cases(Y, T); Y <- Y[ok, , drop = FALSE]; T <- T[ok]
    list(meta = tibble(species = name, study = sp$study, study_name = study_names[sp$study], n = nrow(Y)),
         Sxx = sum(T^2), Syx = as.numeric(t(Y) %*% T), Syy = t(Y) %*% Y)
  })
}) |> unlist(recursive = FALSE)

meta <- map(ss, "meta") |> list_rbind()
S <- nrow(meta)
Syy_arr <- array(0, c(S, 3, 3)); for (s in seq_len(S)) Syy_arr[s, , ] <- ss[[s]]$Syy
standata <- list(S = S, G = 3L, study = meta$study, n = meta$n,
                 Sxx = map_dbl(ss, "Sxx"),
                 Syx = do.call(rbind, map(ss, "Syx")),
                 Syy = Syy_arr)
cat(sprintf("Assembled %d species across %d studies (%s)\n", S, 3,
            paste(table(meta$study_name), collapse = "/")))

# Fit ----
## Sensible inits are essential: with random inits the non-centred loadings start at
## extreme values which, combined with the large-n species, blow up the gradient and
## every transition diverges. Starting near isometry (a~0.4-0.6, small taus) fixes it.
init_fn <- function() list(
  mu_aw = 0.4, tau_aw = 0.1, z_aw = rnorm(S, 0, 0.3),
  mu_a2 = as.array(rep(0.4, 3)),  tau_a2 = 0.1, z_a2 = rnorm(S, 0, 0.3),
  mu_dw = as.array(rep(0, 3)),    tau_dw = 0.1, z_dw = rnorm(S, 0, 0.3),
  mu_d2 = as.array(rep(0, 3)),    tau_d2 = 0.1, z_d2 = rnorm(S, 0, 0.3),
  a_mass = rep(0.6, S), kappa = rep(-0.2, S),
  mu_kappa = as.array(rep(-0.2, 3)), tau_kappa = 0.1,
  psi_mass = rep(0.6, S), psi_wing = rep(0.6, S), psi_app2 = rep(0.6, S))

fit <- stan("Scripts/latent_size_hier.stan", data = standata,
            chains = 4, iter = 2000, warmup = 1000, seed = 1, init = init_fn,
            control = list(adapt_delta = 0.97, max_treedepth = 12))
saveRDS(fit, "Derived/Rdata/hier_sem_fit.rds")

## Diagnostics
cat("\nDivergences:", sum(rstan::get_divergent_iterations(fit)), "\n")
key <- c("mu_aw", "tau_aw", paste0("mu_dw[", 1:3, "]"), "tau_dw")
print(round(summary(fit, pars = key)$summary[, c("mean", "2.5%", "97.5%", "n_eff", "Rhat")], 3))

post <- rstan::extract(fit)
qsum <- function(m) tibble(est = apply(m, 2, mean),
                           lo  = apply(m, 2, quantile, 0.025),
                           hi  = apply(m, 2, quantile, 0.975))

## d_wing = standardized direct wing effect (estimand); a_wing = standardized loading.
hier <- bind_cols(meta, qsum(post$d_wing)) |>
  mutate(lambda_wing = apply(post$a_wing, 2, mean),
         ci_width = hi - lo)
study_mu <- bind_cols(tibble(study_name = study_names), qsum(post$mu_dw))
mu_lw    <- mean(post$mu_aw)
saveRDS(list(hier = hier, study_mu = study_mu, mu_lw = mu_lw), "Derived/hier_sem_estimates.rds")
cat("\nStudy-level mean direct wing effect (mu_bw):\n"); print(study_mu)
cat("\nSaved fit + estimates.\n")
