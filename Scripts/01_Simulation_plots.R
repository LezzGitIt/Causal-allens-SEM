# 01_Simulation_plots.R
# Every figure built from simulated data. Reads only the files exported by
# 00_Simulation_data.R, so it needs no external raw data and runs from a clean clone.
#
# Run order: 00_Simulation_data.R -> this script -> render Scripts/Causal_allens.qmd
#
# Manuscript figures:
#   Figures/dag_comparison.png          Fig 1   three causal assumptions as DAGs
#   Figures/Morphology_distributions.png Fig 2  raw-scale trait distributions
#   Figures/fig_corr_diagnostic.png     Fig 3   induced correlations + realism filter
#   Figures/fig_methods_comparison.png  Fig 4   five estimators vs the true direct effect
#   Figures/fig_lambda_influence.png    Fig 5   effect of the appendage-size loading
#   Figures/fig_covariate_comparison.png Fig 6  oracle vs M1 under mass measurement error
#   Figures/SEM_bias_vs_variance.png    Fig 10  misspecification bias is invariant to N
#
# Figure 11 (Figures/Mass_as_collider.png) is a static asset; see 00_Simulation_data.R.
#
# Diagnostic figures (not in the manuscript):
#   Figures/fig_inspect_variables.png, fig_inspect_mass_me.png,
#   Figures/fig_allom_lines.png, fig_allom_compare.png

library(tidyverse)
library(cowplot)
library(ggpubr)
library(dagitty)
library(ggdag)

ggplot2::theme_set(theme_cowplot())

# Read data ---------------------------------------------------------------

allen_levels    <- c("Strong", "Moderate", "Weak", "No support")
bergmann_levels <- c("Weak", "Moderate", "Strong")

read_results <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(Allen    = factor(Allen,    levels = allen_levels),
           Bergmann = factor(Bergmann, levels = bergmann_levels))
}

results_tbl    <- read_results("Derived/Csv/results_tbl.csv")
results_me_tbl <- read_results("Derived/Csv/results_me_tbl.csv")
parms_mat2     <- read_csv("Derived/Csv/parms_mat.csv",  show_col_types = FALSE)
allom_tbl      <- read_results("Derived/Csv/allom_tbl.csv")
misspec_sim    <- read_csv("Derived/Csv/misspec_sim.csv", show_col_types = FALSE)

sim_morphology <- readRDS("Derived/sim_morphology.rds")
sim_example    <- readRDS("Derived/sim_example.rds")

# Shared labels -----------------------------------------------------------

allen_labs  <- c("No support" = "No support\n(β_A1 = β_A2 = 0.00)",
                 "Weak"       = "Weak\n(β_A1 = β_A2 = 0.02)",
                 "Moderate"   = "Moderate\n(β_A1 = β_A2 = 0.05)",
                 "Strong"     = "Strong\n(β_A1 = β_A2 = 0.07)")
allen_to_bw <- c("Strong" = 0.07, "Moderate" = 0.05, "Weak" = 0.02, "No support" = 0.00)
berg_colors <- c("Weak" = "#92C5DE", "Moderate" = "#4393C3", "Strong" = "#D6604D")

## Reference line: the true direct effect, one value per Allen's level (B_temp_t = B_temp_w).
ref_tbl <- expand_grid(
  Allen     = factor(names(allen_to_bw), levels = allen_levels),
  Appendage = factor(c("Appendage 1", "Appendage 2"), levels = c("Appendage 1", "Appendage 2"))
) %>%
  mutate(true_b = allen_to_bw[as.character(Allen)])

## Standardized coefficients are back-transformed to log mm per degree by multiplying
## by sd(log appendage) / sd(temperature), matching the units of B_temp_w / B_temp_t.
pivot_estimates <- function(tbl, labels) {
  tbl %>%
    pivot_longer(cols = starts_with("coef_"), names_to = "col_name", values_to = "estimate") %>%
    mutate(
      Appendage  = if_else(str_ends(col_name, "_wing"), "Appendage 1", "Appendage 2"),
      method_key = col_name %>% str_remove("^coef_") %>% str_remove("_(wing|tarsus)$"),
      Method     = labels[method_key],
      Method     = factor(Method, levels = unname(labels)),
      Appendage  = factor(Appendage, levels = c("Appendage 1", "Appendage 2")),
      sd_append  = if_else(Appendage == "Appendage 1", sd_log_wing, sd_log_tarsus),
      estimate   = estimate * sd_append / sd_temp
    ) %>%
    filter(!is.na(Method))
}

method_labels <- c("ratio"   = "Appendage / mass",
                   "ryding"  = "Mass as covariate",
                   "sli_iso" = "SLI (isometry)",
                   "sli_est" = "SLI (estimated)",
                   "sem"     = "Latent SEM")

## Every classical method carries an oracle (noiseless latent mass) and an M1 (single
## noisy weighing) variant, so the gap within a pair is the cost of measurement error
## and the gap between pairs is the bias intrinsic to the estimator. The three latent-SEM
## variants come last, so the proposed method reads as the rightmost block of each facet.
method_labels_me <- c(
  "ryding_oracle"  = "Covariate (true mass)",
  "ryding_m1"      = "Covariate (M1)",
  "ratio_oracle"   = "Ratio (true mass)",
  "ratio_m1"       = "Ratio (M1)",
  "sli_est_oracle" = "SLI (true mass)",
  "sli_est_m1"     = "SLI (M1)",
  "sem_oracle"     = "SEM (true mass)",
  "sem_m1"         = "SEM (M1)",
  "sem3"           = "SEM (3x mass)"
)

## Method family drives the boxplot fill; colours match the empirical figures
## (04_Empirical_plots.R) so the same method reads the same across the paper.
me_family <- c(
  "Covariate (true mass)" = "Mass as covariate", "Covariate (M1)" = "Mass as covariate",
  "Ratio (true mass)"     = "Appendage / mass",  "Ratio (M1)"     = "Appendage / mass",
  "SLI (true mass)"       = "SLI (estimated)",   "SLI (M1)"       = "SLI (estimated)",
  "SEM (true mass)"       = "Latent SEM",        "SEM (M1)"       = "Latent SEM",
  "SEM (3x mass)"         = "Latent SEM"
)
me_family_cols <- c("Mass as covariate" = "#FC8D62", "Appendage / mass" = "#E69F00",
                    "SLI (estimated)"   = "#D55E00", "Latent SEM"       = "#8DA0CB")

results_long    <- pivot_estimates(results_tbl,    method_labels)
results_me_long <- pivot_estimates(results_me_tbl, method_labels_me)

# Fig 1: dag_comparison ---------------------------------------------------
# App1 = first appendage (e.g. wing), App2 = second appendage (e.g. tarsus).
# Short names prevent label cutoff in ggdag nodes.

node_fill  <- c("Exposure" = "#FDB462",   # orange    – temperature
                "Outcome"  = "#80B1D3",   # blue      – appendages
                "Latent"   = "#CCCCCC",   # grey      – latent size
                "Observed" = "#F5F5F5")   # off-white – mass / other observed
node_col   <- c("Exposure" = "#9E4A00", "Outcome" = "#1A4F7A",
                "Latent"   = "#555555", "Observed" = "#555555")

make_dag_gg <- function(dag, title) {
  lat <- latents(dag); exp <- exposures(dag); out <- outcomes(dag)
  td  <- tidy_dagitty(dag)
  td$data <- td$data %>%
    mutate(node_class = case_when(
      name %in% lat ~ "Latent",
      name %in% exp ~ "Exposure",
      name %in% out ~ "Outcome",
      TRUE          ~ "Observed"
    ))

  ggplot(td, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(edge_colour = "grey30", edge_width = 0.85,
                   arrow_directed = grid::arrow(length = grid::unit(9, "pt"), type = "closed")) +
    geom_dag_node(aes(colour = node_class, fill = node_class), size = 19) +
    geom_dag_text(colour = "black", size = 3.6, fontface = "bold") +
    scale_fill_manual(values = node_fill,  aesthetics = "fill") +
    scale_colour_manual(values = node_col, aesthetics = "colour") +
    theme_dag(base_size = 11) +
    guides(fill = "none", colour = "none") +
    ggtitle(title) +
    theme(plot.title       = element_text(size = 11, hjust = 0.5, face = "plain",
                                          margin = margin(b = 4)),
          plot.background  = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA))
}

## Model A: mass placed at centre (not top) so Mass->App1 and Temp->App2 don't cross.
model_a <- dagitty('
dag {
  bb = "-1,-1,1,1"
  Temp [exposure, pos = "-0.8,  0.0"]
  Mass [          pos =  "0.0,  0.0"]
  App1 [outcome,  pos =  "0.8,  0.4"]
  App2 [outcome,  pos =  "0.8, -0.4"]
  Temp -> Mass
  Temp -> App1
  Temp -> App2
  Mass -> App1
  Mass -> App2
}')

## Model B: both appendages point INTO Mass, making it a collider.
model_b <- dagitty('
dag {
  bb = "-1,-1,1,1"
  Temp [exposure, pos = "-0.8,  0.0"]
  App1 [outcome,  pos =  "0.8,  0.4"]
  App2 [outcome,  pos =  "0.8, -0.4"]
  Mass [          pos =  "0.1,  0.0"]
  Temp -> App1
  Temp -> App2
  Temp -> Mass
  App1 -> Mass
  App2 -> Mass
}')

## Model C: Size is raised above the Temp-App1 horizontal line so the Size->App1 arrow
## separates from Temp->App1, making both direct effects visible.
model_c <- dagitty('
dag {
  bb = "-1,-1,1,1"
  Temp [exposure, pos = "-0.8,  0.0"]
  Size [latent,   pos = "-0.05, 0.30"]
  Mass [          pos =  "0.8,  0.55"]
  App1 [outcome,  pos =  "0.8,  0.0"]
  App2 [outcome,  pos =  "0.8, -0.55"]
  Temp -> Size
  Temp -> App1
  Temp -> App2
  Size -> Mass
  Size -> App1
  Size -> App2
}')

plot_grid(make_dag_gg(model_a, "(a) Traditional approach"),
          make_dag_gg(model_b, "(b) Mass as collider"),
          make_dag_gg(model_c, "(c) Latent body size"),
          nrow = 1, rel_widths = c(1, 1, 1))
ggsave("Figures/dag_comparison.png", width = 13, height = 4.5, dpi = 150, bg = "white")

# Fig 2: Morphology_distributions -----------------------------------------

sim_morphology %>%
  ggplot(aes(x = Size, color = Species)) +
  geom_density() +
  facet_wrap(~Measurement, scales = "free") +
  labs(x = NULL) +
  guides(color = "none")
ggsave("Figures/Morphology_distributions.png", width = 8, height = 4, dpi = 200, bg = "white")

# Fig 3: fig_corr_diagnostic ----------------------------------------------

corr_colors <- c("TRUE" = "#4DAF4A", "FALSE" = "#E41A1C")

make_corr_panel <- function(df, x_var, y_var, y_lab, hlines, x_lab) {
  df %>%
    mutate(x = .data[[x_var]]) %>%
    ggplot(aes(x = x, y = .data[[y_var]], color = realistic, shape = Allen)) +
    geom_hline(yintercept = hlines, linetype = "dashed", color = "grey60") +
    geom_jitter(width = 0.15, size = 2, alpha = 0.85) +
    facet_wrap(~Bergmann, nrow = 1) +
    scale_color_manual(values = corr_colors) +
    labs(x = x_lab, y = y_lab, color = "Realistic") +
    theme(legend.position = "none", axis.title.x = element_text(size = 9))
}

p_WM <- make_corr_panel(parms_mat2, "lambda_wing",   "r_WM", "r(App1, Mass)", c(0.25, 0.70),  expression(lambda[App1]))
p_WT <- make_corr_panel(parms_mat2, "lambda_wing",   "r_WT", "r(App1, Temp)", c(-0.70, 0.00), expression(lambda[App1]))
p_TM <- make_corr_panel(parms_mat2, "lambda_tarsus", "r_TM", "r(App2, Mass)", c(0.25, 0.70),  expression(lambda[App2]))
p_TT <- make_corr_panel(parms_mat2, "lambda_tarsus", "r_TT", "r(App2, Temp)", c(-0.70, 0.00), expression(lambda[App2]))
p_MT <- make_corr_panel(parms_mat2, "lambda_wing",   "r_MT", "r(Mass, Temp)", c(-0.70, 0.00), expression(lambda[A1]))

top_panels <- ggarrange(p_WM, p_WT, p_TM, p_TT, labels = c("a", "b", "c", "d"),
                        nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")
bot_panel  <- ggarrange(p_MT, labels = "e", ncol = 1, legend = "none")
ggarrange(top_panels, bot_panel, nrow = 2, heights = c(2, 1))
ggsave("Figures/fig_corr_diagnostic.png", width = 11, height = 10, dpi = 300, bg = "white")

# Fig 4: fig_methods_comparison -------------------------------------------

results_long %>%
  ggplot(aes(x = Method, y = estimate)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_hline(data = ref_tbl, aes(yintercept = true_b),
             color = "#33A02C", linetype = "dashed", linewidth = 0.9) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Bergmann), width = 0.15, alpha = 0.55, size = 1.5) +
  facet_grid(Appendage ~ Allen, labeller = labeller(Allen = allen_labs), scales = "free_y") +
  scale_color_manual(values = berg_colors) +
  labs(x = NULL, y = expression(beta[T] ~ "(log mm/°C)"), color = "Bergmann's rule") +
  theme(axis.text.x = element_text(angle = 55, vjust = 0.6, size = 8),
        strip.text  = element_text(size = 8.5), legend.position = "top")
ggsave("Figures/fig_methods_comparison.png", width = 9, height = 7, dpi = 300, bg = "white")

# Fig 5: fig_lambda_influence ---------------------------------------------

results_long %>%
  filter(method_key %in% c("sem", "ryding")) %>%
  mutate(
    lambda = if_else(Appendage == "Appendage 1", lambda_wing, lambda_tarsus),
    Method = factor(method_labels[method_key], levels = c("Mass as covariate", "Latent SEM"))
  ) %>%
  ggplot(aes(x = lambda, y = estimate, color = Method)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
  geom_hline(data = ref_tbl, aes(yintercept = true_b),
             color = "#33A02C", linetype = "dashed", linewidth = 0.9) +
  geom_point(alpha = 0.35, size = 1.5) +
  geom_smooth(method = "loess", span = 0.9, se = TRUE, linewidth = 1, alpha = 0.15) +
  facet_grid(Appendage ~ Allen, labeller = labeller(Allen = allen_labs), scales = "free_y") +
  scale_color_manual(values = c("Mass as covariate" = "#FC8D62", "Latent SEM" = "#8DA0CB")) +
  labs(x = expression("Appendage loading " * (lambda)),
       y = expression(beta[T] ~ "(log mm/°C)"), color = "Method") +
  theme(strip.text = element_text(size = 8.5), legend.position = "top")
ggsave("Figures/fig_lambda_influence.png", width = 12, height = 6, dpi = 300, bg = "white")

# Fig 6: fig_covariate_comparison -----------------------------------------

results_me_long %>%
  mutate(approach = factor(me_family[as.character(Method)], levels = names(me_family_cols))) %>%
  ggplot(aes(x = Method, y = estimate, fill = approach)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_hline(data = ref_tbl, aes(yintercept = true_b),
             color = "#33A02C", linetype = "dashed", linewidth = 0.9) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_jitter(aes(color = Bergmann), width = 0.15, alpha = 0.5, size = 1.2) +
  facet_grid(Appendage ~ Allen, labeller = labeller(Allen = allen_labs), scales = "free_y") +
  scale_fill_manual(values = me_family_cols) +
  scale_color_manual(values = berg_colors) +
  labs(x = NULL, y = expression(beta[T] ~ "(log mm/°C)"),
       fill = "Method", color = "Bergmann's rule") +
  guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1)) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 7.5),
        strip.text      = element_text(size = 8.5),
        legend.position = "top", legend.box = "vertical", legend.spacing.y = unit(0, "pt"))
ggsave("Figures/fig_covariate_comparison.png", width = 13, height = 7.5, dpi = 300, bg = "white")

# Fig 10: SEM_bias_vs_variance (Supporting Information) --------------------
## Bias is the single-factor SEM's wing coefficient minus an oracle regression on the
## true latent Size. It is ~0 at every N when the model is correctly specified (whether
## or not mass is noisy), but constant and large when a second appendage factor exists.
## The standard errors shrink as 1/sqrt(N) in all three scenarios: more data buys
## precision, never correctness.

scen_cols <- c("Well-specified" = "#009E73", "Noisy mass" = "#E69F00",
               "2nd (limb) factor" = "#D55E00")

misspec_sim %>%
  pivot_longer(c(Bias, `SE(B_temp)`, `SE(lambda)`), names_to = "quantity", values_to = "value") %>%
  mutate(quantity = factor(quantity, levels = c("Bias", "SE(B_temp)", "SE(lambda)")),
         scenario = factor(scenario, levels = names(scen_cols))) %>%
  ggplot(aes(N, value, color = scenario)) +
  geom_hline(data = ~filter(.x, quantity == "Bias"), aes(yintercept = 0),
             linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_x_log10() +
  scale_color_manual(values = scen_cols, name = NULL) +
  facet_wrap(~quantity, scales = "free_y") +
  labs(x = "Sample size per species (log scale)", y = NULL,
       title = "Sample size fixes imprecision (SE), not misspecification bias",
       subtitle = "Bias in the direct wing effect stays constant with N only when a 2nd appendage factor is present") +
  theme(legend.position = "top")
ggsave("Figures/SEM_bias_vs_variance.png", width = 13, height = 5, dpi = 200, bg = "white")

# Diagnostics (not in the manuscript) --------------------------------------

df_ex    <- sim_example$df_ex
df_ex_me <- sim_example$df_ex_me

## Marginal distributions of the example dataset.
df_ex %>%
  dplyr::select(Temp_inc, Mass, Append, Tarsus) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = c("Temp_inc", "Mass", "Append", "Tarsus"),
                           labels = c("Temperature", "Mass (g)", "App1 (mm)", "App2 (mm)"))) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white", linewidth = 0.2) +
  facet_wrap(~variable, scales = "free") +
  labs(title = sprintf("Example dataset: B_temp_w = %.2f, beta_size = %.2f",
                       sim_example$B_temp_w, sim_example$beta_size),
       x = NULL, y = "Count")
ggsave("Figures/fig_inspect_variables.png", width = 8, height = 5, dpi = 200, bg = "white")

## How much measurement error widens the repeat mass measurements.
bind_rows(
  df_ex    %>% dplyr::select(M1, M2, M3) %>% mutate(sigma_me = "No error"),
  df_ex_me %>% dplyr::select(M1, M2, M3) %>% mutate(sigma_me = "sigma_me = 0.10 (log)")
) %>%
  pivot_longer(c(M1, M2, M3), names_to = "rep", values_to = "mass") %>%
  ggplot(aes(x = mass, fill = sigma_me)) +
  geom_histogram(bins = 40, position = "identity", alpha = 0.6, color = "white", linewidth = 0.2) +
  facet_wrap(~rep) +
  scale_fill_manual(values = c("No error" = "#4393C3", "sigma_me = 0.10 (log)" = "#D6604D")) +
  labs(title = "Mass measurements: effect of measurement error", x = "Mass (g)", y = "Count",
       fill = expression(sigma[me]))
ggsave("Figures/fig_inspect_mass_me.png", width = 8, height = 4, dpi = 200, bg = "white")

## OLS vs SMA lines on the example dataset. SMA exceeds OLS because mass carries error.
df_ex %>%
  pivot_longer(c(Append_log, Tarsus_log), names_to = "appendage", values_to = "app_log") %>%
  mutate(Appendage = if_else(appendage == "Append_log", "App1 (wing)", "App2 (tarsus)"),
         Appendage = factor(Appendage, levels = c("App1 (wing)", "App2 (tarsus)"))) %>%
  ggplot(aes(x = Mass_log, y = app_log)) +
  geom_point(alpha = 0.15, size = 0.9, color = "grey40") +
  geom_abline(data = sim_example$line_tbl,
              aes(slope = slope, intercept = intercept, color = Method, linetype = Method),
              linewidth = 1.1) +
  scale_color_manual(values = c("OLS" = "#D6604D", "SMA" = "#4393C3")) +
  scale_linetype_manual(values = c("OLS" = "dashed", "SMA" = "solid")) +
  facet_wrap(~Appendage, scales = "free") +
  labs(x = "log(Mass)", y = "log(Appendage)",
       title = sprintf("B_temp_w = %.2f, beta_size = %.2f",
                       sim_example$B_temp_w, sim_example$beta_size)) +
  theme(legend.position = "top")
ggsave("Figures/fig_allom_lines.png", width = 8, height = 4.5, dpi = 250, bg = "white")

## SMA vs OLS slope across the grid, as a function of the appendage loading.
allom_long <- bind_rows(
  allom_tbl %>% transmute(Appendage = "App1 (wing)",   lambda = lambda_wing,
                          Allen, Bergmann, OLS_slope = ols_slope_wing,   SMA_slope = sma_slope_wing),
  allom_tbl %>% transmute(Appendage = "App2 (tarsus)", lambda = lambda_tarsus,
                          Allen, Bergmann, OLS_slope = ols_slope_tarsus, SMA_slope = sma_slope_tarsus)
) %>%
  mutate(Appendage = factor(Appendage, levels = c("App1 (wing)", "App2 (tarsus)")))

allom_long %>%
  pivot_longer(c(OLS_slope, SMA_slope), names_to = "method", values_to = "slope") %>%
  mutate(method = factor(if_else(method == "OLS_slope", "OLS", "SMA"), levels = c("OLS", "SMA"))) %>%
  ggplot(aes(x = lambda, y = slope, fill = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8, position = position_dodge(0.7)) +
  scale_fill_manual(values = c("OLS" = "#D6604D", "SMA" = "#4393C3")) +
  facet_wrap(~Appendage, scales = "free") +
  labs(x = expression("Appendage loading " * (lambda)), y = "Slope (log-log)",
       fill = "Method", title = "Allometric slope by loading and method") +
  theme(legend.position = "top")
ggsave("Figures/fig_allom_compare.png", width = 9, height = 4.5, dpi = 250, bg = "white")

cat("Saved all simulation figures to Figures/\n")
