# Fig11_collider.R
# Regenerates Figures/Mass_as_collider.png (manuscript Figure 11): the over-adjustment /
# collider illustration. Temperature vs log(appendage) stratified into mass quartiles.
#
# This figure uses the MULTIVARIATE-NORMAL simulation (sliR::sim_allometric), NOT the
# latent-size causal DAG that the rest of the paper simulates (gen_causal_data). The two
# are different models and the distinction is the point of the figure: the sign flip it
# shows — temperature has ~no marginal effect on appendage length yet a spuriously
# NEGATIVE effect once mass is held constant — is a suppression/over-adjustment pattern
# that a correlation-specified model produces but the single-factor latent-size DAG
# cannot (there both slopes share the sign of lambda * beta_size). The original figure
# was built this way; the code is ported from the SLI methods paper
# (SLI-allens-rule/Scripts/qmd/SMA_body_shape_methods.qmd, chunk `fig-condition-on-mass`).
#
# Dependency: sliR (github.com/LezzGitIt/sliR) — install with
#   remotes::install_github("LezzGitIt/sliR")
#
# Run from the project root:  Rscript Scripts/Fig11_collider.R

library(tidyverse)
library(sliR)
ggplot2::theme_set(cowplot::theme_cowplot())

## Inverse allometry (r_app_mass < 0) with no direct temperature-appendage effect
## (r_grad_app = 0) and a strong Bergmann temperature-mass relationship
## (r_grad_mass < 0). Warming shrinks mass, and because appendage scales inversely with
## mass, appendage lengthens: the effect of temperature runs entirely through mass.
set.seed(2024)
collider_df <- sim_allometric(
  n             = 1500,
  b_avg         = 0.33,
  r_app_mass    = -0.40,          # inverse allometry
  gradient      = "Temp_inc",
  gradient_dist = "normal",
  mean_gradient = 1,
  sd_gradient   = 0.18,
  r_grad_app    = 0.00,           # NO direct temperature effect on appendage
  r_grad_mass   = -0.70,          # strong Bergmann's rule
  trim_sd       = 3
)

## The contrast the figure exists to show: marginal temperature slope ~ 0 (or slightly
## positive), but conditioning on mass induces a spurious NEGATIVE slope.
marg <- coef(lm(Append_log ~ Temp_inc,            data = collider_df))["Temp_inc"]
cond <- coef(lm(Append_log ~ Mass_log + Temp_inc, data = collider_df))["Temp_inc"]
cat(sprintf("marginal Temp slope %+.3f; mass-adjusted %+.3f\n", marg, cond))
stopifnot(cond < 0, cond < marg)   # over-adjustment: conditioning drives the slope down

## Facet by mass quartile; within each stratum the Temp-appendage relationship is flat
## to negative, the visual signature of conditioning on mass.
lab_tbl <- collider_df %>%
  mutate(Mass_bin = cut_number(Mass_log, n = 4)) %>%
  distinct(Mass_bin) %>%
  mutate(lab = Mass_bin %>% as.character() %>%
           gsub("\\(|\\]|\\[", "", .) %>%
           strsplit(",") %>%
           sapply(\(x) paste("log(Mass) =", x[1], "-", x[2])))
lab_map <- setNames(lab_tbl$lab, lab_tbl$Mass_bin)

collider_df %>%
  mutate(Mass_bin = cut_number(Mass_log, n = 4)) %>%
  ggplot(aes(x = Temp_inc, y = Append_log)) +
  geom_point(aes(color = Mass_log), alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ Mass_bin, labeller = labeller(Mass_bin = lab_map)) +
  scale_color_viridis_c() +
  labs(x = "Temperature increase", y = "log(Appendage)", color = "Mass")
ggsave("Figures/Mass_as_collider.png", width = 8, height = 6, dpi = 200, bg = "white")
cat("Saved Figures/Mass_as_collider.png\n")
