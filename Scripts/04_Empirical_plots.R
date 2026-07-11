# 04_Empirical_plots.R
# Every figure built from the empirical datasets. Reads only the saved estimates
# (Derived/hier_sem_estimates.rds from 03, Derived/blavaan_empirical_estimates.csv from
# 02) and refits nothing, so it runs from a clean clone WITHOUT the external raw data.
#
# Run order: 02_Empirical_SEM.R -> 03_Hierarchical_SEM.R -> this script
#
# Manuscript figures:
#   Figures/Empirical_all_estimates.png  Fig 7  every approach per species, by Bergmann direction
#   Figures/report_study_means.png       Fig 8  study-level population mean direct App1 effect
#   Figures/report_ci_vs_n.png           Fig 9  per-species interval width vs sample size
#
# Diagnostic figures (not in the manuscript):
#   Figures/Hier_shrinkage_precision.png, Figures/Hier_population_loadings.png
#
# Note App1 = wing in all three datasets; App2 = tail (nightjars) or tarsus (others).

library(tidyverse)
library(cowplot)
library(patchwork)
ggplot2::theme_set(theme_cowplot(font_size = 11))

study_colors <- c("Nightjar" = "#E41A1C", "Weeks (2020)" = "#377EB8", "Atlantic birds" = "#4DAF4A")
study_cap    <- "Species labels coloured by study: Nightjar (red), Weeks 2020 (blue), Atlantic birds (green)"

he       <- readRDS("Derived/hier_sem_estimates.rds")
hier     <- he$hier
study_mu <- he$study_mu
mu_lw    <- he$mu_lw

blav <- read_csv("Derived/blavaan_empirical_estimates.csv", show_col_types = FALSE) %>%
  filter(appendage == "wing")
dir_tbl <- blav %>% distinct(species, Direction)

# Fig 7: all approaches per species, grouped by Bergmann direction ---------
## One figure with every estimate + CI: the reference metrics, the per-species blavaan
## SEM, and the hierarchical SEM, grouped by the direction of the confounding the SEM
## exists to remove.
approach_cols <- c("Ratio" = "#E69F00", "SLI isometry" = "#CC79A7", "SLI estimated" = "#D55E00",
                   "Mass as covariate" = "#0072B2", "SEM (blavaan)" = "#009E73",
                   "Hierarchical SEM" = "#9400D3")

other <- blav %>%
  filter(Approach %in% c("Ratio", "Sli_iso", "Sli_est", "Ryding", "SEM"), std.error < 0.5) %>%
  transmute(species, study, Direction,
            Approach = recode(Approach, Sli_iso = "SLI isometry", Sli_est = "SLI estimated",
                              Ryding = "Mass as covariate", SEM = "SEM (blavaan)"),
            estimate, LCI95, UCI95)
hierL <- hier %>%
  transmute(species, study = study_name, Approach = "Hierarchical SEM",
            estimate = est, LCI95 = lo, UCI95 = hi) %>%
  left_join(dir_tbl, by = "species")
comp <- bind_rows(other, hierL) %>%
  mutate(Approach = factor(Approach, levels = names(approach_cols)),
         study    = factor(study, levels = names(study_colors)))

## Flag species whose hierarchical SEM estimate is bracketed by the two reference metrics
## (Ratio and Mass-as-covariate) — i.e. the SEM lands between them rather than outside both.
bracket <- comp %>%
  filter(Approach %in% c("Ratio", "Mass as covariate", "Hierarchical SEM")) %>%
  dplyr::select(species, study, Direction, Approach, estimate) %>%
  pivot_wider(names_from = Approach, values_from = estimate) %>%
  mutate(between = `Hierarchical SEM` >= pmin(Ratio, `Mass as covariate`) &
                   `Hierarchical SEM` <= pmax(Ratio, `Mass as covariate`))
## Asterisk y-position: just above the tallest CI drawn for each species.
ast_pos <- comp %>% group_by(species) %>%
  summarise(y = max(UCI95, na.rm = TRUE), .groups = "drop")
bracket <- bracket %>% left_join(ast_pos, by = "species")

prop_tbl <- bracket %>% filter(!is.na(between)) %>%
  group_by(Direction) %>%
  summarise(n_between = sum(between), n_total = n(), prop_between = mean(between), .groups = "drop")

panel <- function(dir) {
  d    <- comp %>% filter(Direction == dir)
  meta <- d %>% group_by(species, study) %>% summarise(m = mean(estimate), .groups = "drop") %>%
    arrange(study, m)
  ast  <- bracket %>% filter(Direction == dir, between) %>%
    mutate(species = factor(species, levels = meta$species))
  d %>% mutate(species = factor(species, levels = meta$species)) %>%
    ggplot(aes(species, estimate, color = Approach)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(ymin = LCI95, ymax = UCI95), width = 0, alpha = 0.8,
                  position = position_dodge(0.75)) +
    geom_point(size = 1.1, position = position_dodge(0.75)) +
    geom_text(data = ast, aes(x = species, y = y + 0.03), label = "*",
              inherit.aes = FALSE, size = 5, vjust = 0) +
    scale_color_manual(values = approach_cols) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 2.5))) +
    labs(title = dir, x = NULL, y = expression(beta[Temp] ~ "on App1 (std.)"), color = NULL) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5,
                                     color = study_colors[as.character(meta$study)]),
          legend.position = "top")
}
fig_all <- (panel("Bergmann's") / panel("Inverse Bergmann's")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", caption = study_cap) &
  theme(legend.position = "top", plot.caption = element_text(hjust = 0, size = 9))
ggsave("Figures/Empirical_all_estimates.png", fig_all, bg = "white", width = 13, height = 8, dpi = 200)
cat("Saved Figures/Empirical_all_estimates.png\n")
cat("Proportion of species with hierarchical SEM bracketed by Ratio & Mass-as-covariate:\n")
print(prop_tbl)

# Fig 8: study-level "average Allen's rule" --------------------------------
## Population mean direct App1 effect per study. The Weeks (2020) effect is per year;
## the other two are per degree C, so magnitudes are not comparable across studies.

p_means <- study_mu %>% mutate(study_name = factor(study_name, levels = names(study_colors))) %>%
  ggplot(aes(est, study_name, color = study_name)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_color_manual(values = study_colors, guide = "none") +
  labs(x = "population mean direct App1 effect  (95% CI, standardized)", y = NULL)
ggsave("Figures/report_study_means.png", p_means, bg = "white", width = 7, height = 2.6, dpi = 200)

# Fig 9: SEM interval width vs per-species sample size ---------------------
## Imprecision is a sample-size problem. Intervals are the per-species blavaan HPD
## widths, NOT the hierarchical SEM (which borrows strength across species and so partly
## decouples width from a species' own N). Per-study slopes plus the pooled slope show
## the relationship holds WITHIN each dataset, not just across them.

n_tbl <- hier %>% dplyr::select(species, study = study_name, n)
ci_ps <- blav %>% filter(Approach == "SEM") %>% transmute(species, ci_width = UCI95 - LCI95)

d <- n_tbl %>% left_join(ci_ps, by = "species") %>%
  mutate(study = factor(study, levels = names(study_colors)))
r <- cor(d$ci_width, log(d$n), use = "complete.obs")

p_ci <- d %>% ggplot(aes(n, ci_width, color = study)) +
  geom_point(alpha = 0.8, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linewidth = 0.9) +
  scale_x_log10() + scale_color_manual(values = study_colors, name = NULL) +
  annotate("text", x = min(d$n), y = min(d$ci_width), hjust = 0, vjust = 0, size = 3.4,
           label = sprintf("pooled r = %.2f with log N", r)) +
  labs(x = "sample size per species (log scale)", y = "per-species SEM 95% interval width",
       caption = "Black line = pooled fit across all species; coloured lines = within-study fits") +
  theme(legend.position = "top", plot.caption = element_text(hjust = 0, size = 8))
ggsave("Figures/report_ci_vs_n.png", p_ci, bg = "white", width = 7, height = 4.0, dpi = 200)

cat("Saved Figures/report_study_means.png and Figures/report_ci_vs_n.png\n")
cat(sprintf("interval width vs log N: pooled r = %.2f\n", r))
cat("within-study slopes (ci_width ~ log10 N):\n")
print(d %>% group_by(study) %>%
        summarise(n_spp = n(), slope = coef(lm(ci_width ~ log10(n)))[2],
                  r = cor(ci_width, log(n)), .groups = "drop"))

# Diagnostics (not in the manuscript) --------------------------------------

## Shrinkage + precision: per-species (blavaan) vs hierarchical.
sem_ps <- blav %>% filter(Approach == "SEM") %>%
  transmute(species, ps_est = estimate, ps_ci = UCI95 - LCI95)
shr <- hier %>%
  transmute(species, study = study_name, hier_est = est, hier_ci = ci_width, n) %>%
  left_join(sem_ps, by = "species") %>%
  mutate(study = factor(study, levels = names(study_colors)))

## Horizontal lines = each study's population mean: hierarchical estimates (y) are pulled
## toward their study's line, while per-species estimates (x) spread.
study_lines <- study_mu %>% mutate(study = factor(study_name, levels = names(study_colors)))
pA <- shr %>% ggplot(aes(ps_est, hier_est, color = study)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(data = study_lines, aes(yintercept = est, color = study),
             linetype = "dotted", linewidth = 0.8, show.legend = FALSE) +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_color_manual(values = study_colors, name = NULL) +
  coord_cartesian(xlim = c(-0.6, 0.8), ylim = c(-0.6, 0.8)) +
  labs(x = "per-species SEM estimate", y = "hierarchical SEM estimate",
       title = "(a) Shrinkage toward study means (dotted lines)")
pB <- shr %>% ggplot(aes(ps_ci, hier_ci, color = study)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_color_manual(values = study_colors, name = NULL) +
  labs(x = "per-species SEM 95% interval width", y = "hierarchical 95% interval width",
       title = "(b) ...and tightens intervals (points below the line)")
ggsave("Figures/Hier_shrinkage_precision.png",
       (pA | pB) + plot_layout(guides = "collect") & theme(legend.position = "top"),
       bg = "white", width = 12, height = 5.5, dpi = 200)
cat("Saved Figures/Hier_shrinkage_precision.png\n")
cat(sprintf("median CI width: per-species %.3f -> hierarchical %.3f (%.0f%% narrower)\n",
            median(shr$ps_ci, na.rm = TRUE), median(shr$hier_ci, na.rm = TRUE),
            100 * (1 - median(shr$hier_ci, na.rm = TRUE) / median(shr$ps_ci, na.rm = TRUE))))

## Study-level means + global pooling of the App1 loading.
pC <- study_mu %>% mutate(study_name = factor(study_name, levels = names(study_colors))) %>%
  ggplot(aes(est, study_name, color = study_name)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15) +
  geom_point(size = 3) +
  scale_color_manual(values = study_colors, guide = "none") +
  labs(x = expression("population mean direct App1 effect" ~ mu[d1] ~ "(95% CI)"),
       y = NULL, title = "(a) Study-level 'average Allen's rule'")
pD <- hier %>% mutate(study = factor(study_name, levels = names(study_colors))) %>%
  ggplot(aes(n, lambda_wing, color = study)) +
  geom_hline(yintercept = mu_lw, linetype = "dashed", color = "grey40") +
  annotate("text", x = min(hier$n), y = mu_lw, vjust = -0.6, hjust = 0, size = 3,
           color = "grey40", label = "global pooled mean") +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_x_log10() + scale_color_manual(values = study_colors, name = NULL) +
  labs(x = "sample size per species (log scale)",
       y = "standardized App1 loading (cor with latent size)",
       title = "(b) Loadings pooled globally (weakly-coupled species shrink to the mean)")
ggsave("Figures/Hier_population_loadings.png",
       (pC | pD) + plot_layout(widths = c(1, 1.2)) & theme(legend.position = "top"),
       bg = "white", width = 12, height = 5, dpi = 200)
cat("Saved Figures/Hier_population_loadings.png\n")
