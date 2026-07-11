# Reevaluating Allen's rule: a causal modeling perspective

**Aaron Skinner**

## Overview

This repository contains the analysis code for a methods paper on estimating the direct
effect of temperature on appendage length — Allen's rule, the ecogeographic pattern of
relatively longer appendages in warmer climates. Testing the rule requires controlling
for body size, but its usual proxy, body mass, is a noisy indicator of a latent
body-size factor rather than a cause of the appendages. Conditioning on observed mass
therefore induces collider bias.

I use causal directed acyclic graphs and a factorial simulation to demonstrate this bias,
compare four analytical approaches, and propose a latent structural equation model (SEM)
that estimates body size from the shared variance among mass and multiple appendages. The
approaches are then applied to three empirical bird datasets, with a hierarchical Bayesian
SEM that pools information across 95 species and three studies.

## Repository structure

```
Scripts/
  qmd/
    Causal_allens.qmd          # Manuscript (Quarto)
  Key_causal_fns.R             # Shared functions (sourced by the numbered scripts)
  00_Simulation_data.R         # Simulate under the causal DAG; fit all methods
  01_Simulation_plots.R        # Simulation figures
  02_Empirical_SEM.R           # Per-species lavaan / blavaan SEM on the 3 datasets
  03_Hierarchical_SEM.R        # Hierarchical latent-size SEM (latent_size_hier.stan)
  04_Empirical_plots.R         # Empirical figures
  Fig11_collider.R             # Collider / over-adjustment illustration
  latent_size_hier.stan        # Stan model for the hierarchical SEM

Suppfiles/                     # Bibliography, journal metadata, title-page partial
_extensions/                   # Quarto elsevier journal-format extension (needed to render)
```

Manuscript `.qmd` files live in `Scripts/qmd/`, separate from the analysis `.R` scripts,
so that Quarto's per-render byproducts don't clutter `Scripts/`. All scripts use paths
relative to the project root (`_quarto.yml` sets `execute-dir: project`).

## Reproducing the analysis

Run the simulation scripts (no external data required), then, after obtaining the raw
data (see Data availability), the empirical scripts:

``` r
source("Scripts/00_Simulation_data.R")
source("Scripts/01_Simulation_plots.R")
source("Scripts/Fig11_collider.R")
source("Scripts/02_Empirical_SEM.R")
source("Scripts/03_Hierarchical_SEM.R")
source("Scripts/04_Empirical_plots.R")
```

Then render the manuscript from the project root:

``` bash
quarto render Scripts/qmd/Causal_allens.qmd
```

## Approaches compared

| Approach | Size adjustment | Unbiased in simulation? |
|---|---|---|
| Appendage/mass ratio | `log(App) − 1·log(Mass)` | No — strong positive bias |
| SLI (isometry) | `log(App) − (1/3)·log(Mass)` | No — moderate bias |
| SLI (estimated exponent) | `log(App) − b_SMA·log(Mass)` | Nearly, with true mass |
| Mass as covariate | `App ~ Temp + Mass` | No — collider bias, negative |
| Latent SEM | latent Size from shared variance | Yes |

The SLI is one approach fitted with two exponents, so the four approaches yield five
estimators. Under measurement error each appears in an *oracle* (noiseless mass) and *M1*
(single noisy weighing) variant, separating estimator bias from the cost of measurement
error.

## Data availability

- **Nightjar data**: Skinner et al. (2025) *Journal of Biogeography*
  <https://doi.org/10.1111/jbi.15176>. Data at
  <https://datadryad.org/dataset/doi:10.5061/dryad.pnvx0k6xw>.
- **Weeks et al. (2020)**: Shared morphological consequences of global warming in North
  American migratory birds. *Ecology Letters* 23, 316–325
  <https://doi.org/10.1111/ele.13434>. Data at
  <https://datadryad.org/dataset/doi:10.5061/dryad.8pk0p2nhw>.
- **Atlantic bird traits**: ATLANTIC BIRD TRAITS dataset. Rodrigues et al. (2019)
  *Ecology* 100, e02647 <https://doi.org/10.1002/ecy.2647>.

## Dependencies

R packages: `tidyverse`, `lavaan`, `blavaan`, `rstan`, `smatr`, `cowplot`, `patchwork`,
`ggpubr`, `dagitty`, `ggdag`, `MASS`, `janitor`, and `sliR` (Figure 11 only;
`remotes::install_github("LezzGitIt/sliR")`).

## Citation

[To be added upon publication]
