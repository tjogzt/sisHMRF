# sisHMRF

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.21287416.svg?v=1)](https://doi.org/10.5281/zenodo.21287416)

**Spatial Immune Scoring with Hidden Markov Random Field**

A VI-EM framework for deconvolving spatially-organized immune states in spatial transcriptomics data.

## Installation

```r
remotes::install_github("tjogzt/sisHMRF")
```

## Quick Start

```r
library(sisHMRF)

# Load example BRCA Visium data (500 spots, 28 immune features)
data("sisHMRF_data")

# Fit the model
fit <- sis_fit(
  Y = sisHMRF_data$F_mat,
  adjacency = build_spatial_graph(sisHMRF_data$coords, method = "knn", k = 6),
  P = 5, n_init = 3
)

# Visualize results
sis_plot_scores(fit, coords = sisHMRF_data$coords)
sis_plot_panel(fit, coords = sisHMRF_data$coords)
```

## Features

- **VI-EM core engine** — Variational inference with expectation-maximization for spatial HMRF with immune-informed priors
- **28 immune features** — Across functional categories (checkpoint, cytolytic, antigen presentation, etc.)
- **Spatial scoring** — Immune Archetype Score (IAS), Immune State Score (ISS), Immune Boundary Index (IBI) with empirical Bayesian confidence intervals
- **Model selection** — BIC-based selection of the number of spatial domains (sis_select_P)
- **Publication-ready visualizations** — 8 ggplot2-based plot types

## Citation

```r
citation("sisHMRF")
```

## License

MIT © Tao Zhu
