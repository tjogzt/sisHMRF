# sisHMRF

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.21287416.svg?v=1)](https://doi.org/10.5281/zenodo.21287416)

**Spatial Immune State Inference via Hidden Markov Random Fields**

Continuous-latent-variable hidden Markov random field for inferring immune
state landscapes from spatial transcriptomics data.

## Installation

```r
remotes::install_github("tjogzt/sisHMRF")
```

## Quick Start

```r
library(sisHMRF)

# Build immune feature matrix from expression data
feats <- sis_features(expr_mat)
#       returns list(F_mat, feature_names, feature_categories, n_features)

# Fit the spatial immune state model
fit <- sis_fit(F_mat = feats$F_mat, coords = coords, P = 3, k = 6)

# Visualize results
plot(fit, coords, type = "panel")       # latent dimensions
plot(fit, coords, type = "archetypes")  # immune archetypes (IS1-IS4)
plot(fit, coords, type = "gradient")    # immune transition zones

# Extract immune scores
scores <- sis_scores(fit)  # IAS, ISS, IBI
```

## Features

- **GMRF-based continuous latent states** — models immune archetypes as a smooth spatial continuum rather than discrete clusters
- **26 immune features** — 17 pathway signatures + 9 cell-type scores covering the full cancer-immunity cycle
- **Spatial scoring** — Immune Activation Score (IAS), Immune Suppression Score (ISS), Immune Balance Index (IBI)
- **Model selection** — BIC-based selection of latent dimension count (`sis_select_P`)
- **Publication-ready visualizations** — ggplot2-based plots with 7pt base theme

## Citation

```r
citation("sisHMRF")
```

## License

MIT © Tao Zhu
