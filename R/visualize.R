#' Visualization for sis-HMRF results
#'
#' @description
#' Publication-quality ggplot2 visualizations for continuous spatial immune
#' states inferred by sis-HMRF. All plotting functions return ggplot objects
#' suitable for further customization or export via `ggplot2::ggsave()`.
#'
#' @name sis_plots
#' @importFrom ggplot2 aes
NULL

# Silence R CMD check about ggplot2 aes() variables
utils::globalVariables(c(
  "x", "y", "z", "se", "gradient", "archetype",
  "IAS", "ISS", "IBI",
  "Dimension", "Feature", "Loading",
  "iteration", "elbo"
))


# ---- Color palettes ----

#' @keywords internal
#' @noRd
.sis_pal_div <- function(...) {
  ggplot2::scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7",
                                  high = "#B2182B", midpoint = 0, ...)
}

#' @keywords internal
#' @noRd
.sis_fill_div <- function(...) {
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7",
                                 high = "#B2182B", midpoint = 0, ...)
}

#' @keywords internal
#' @noRd
.sis_pal_seq <- function(...) {
  ggplot2::scale_color_viridis_c(option = "plasma", ...)
}

#' @keywords internal
#' @noRd
.sis_theme_map <- function(base_size = 7) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = base_size, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = base_size - 1),
      legend.title     = ggplot2::element_text(size = base_size - 1),
      legend.text      = ggplot2::element_text(size = base_size - 2),
      legend.key.size  = ggplot2::unit(0.3, "cm"),
      panel.grid       = ggplot2::element_blank(),
      axis.text        = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank()
    )
}


# ---- Plot functions ----

#' Plot a single latent immune dimension
#'
#' @param fit A fitted `sisHMRF` object.
#' @param coords N x 2 matrix of spatial coordinates.
#' @param dim_idx Integer. Which latent dimension to plot (1-indexed).
#' @param title Optional plot title.
#' @param point_size Point size for scatter (default: 1.5).
#'
#' @return A `ggplot` object.
#'
#' @export
sis_plot_dimension <- function(fit, coords, dim_idx = 1, title = NULL,
                                point_size = 1.5) {
  if (is.null(title)) {
    title <- paste("Immune State Dimension", dim_idx)
  }

  df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    z = fit$Z_mean[, dim_idx]
  )

  ggplot2::ggplot(df, aes(x = x, y = y, color = z)) +
    ggplot2::geom_point(size = point_size) +
    .sis_pal_div() +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    .sis_theme_map()
}


#' Plot all latent immune dimensions as a panel
#'
#' @param fit A fitted `sisHMRF` object.
#' @param coords N x 2 matrix of spatial coordinates.
#' @param dim_labels Optional character vector of dimension labels.
#' @param ncol Number of columns in the panel (default: 2).
#'
#' @return A `ggplot` object (combined via patchwork) or list of ggplots.
#'
#' @export
sis_plot_panel <- function(fit, coords, dim_labels = NULL, ncol = 2) {
  P <- ncol(fit$Z_mean)
  if (is.null(dim_labels)) dim_labels <- paste0("Dim ", seq_len(P))

  plots <- lapply(seq_len(P), function(p) {
    sis_plot_dimension(fit, coords, dim_idx = p, title = dim_labels[p])
  })

  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(plots, ncol = ncol)
  } else {
    plots
  }
}


#' Plot immune archetypes (post-hoc clustering)
#'
#' @param fit A fitted `sisHMRF` object.
#' @param coords N x 2 matrix of spatial coordinates.
#' @param K Number of archetypes (default: 4).
#'
#' @return A `ggplot` object.
#'
#' @export
sis_plot_archetypes <- function(fit, coords, K = 4) {
  clusters <- cluster_immune_states(fit, K = K)

  df <- data.frame(
    x         = coords[, 1],
    y         = coords[, 2],
    archetype = clusters
  )

  ggplot2::ggplot(df, aes(x = x, y = y, color = archetype)) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_color_brewer(palette = "Set2", name = "Archetype") +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = "Immune State Archetypes", x = NULL, y = NULL) +
    .sis_theme_map()
}


#' Plot immune gradient map
#'
#' @param fit A fitted `sisHMRF` object.
#' @param coords N x 2 matrix of spatial coordinates.
#'
#' @return A `ggplot` object.
#'
#' @export
sis_plot_gradient <- function(fit, coords) {
  grad <- compute_immune_gradient(fit)

  df <- data.frame(
    x        = coords[, 1],
    y        = coords[, 2],
    gradient = grad
  )

  ggplot2::ggplot(df, aes(x = x, y = y, color = gradient)) +
    ggplot2::geom_point(size = 1.5) +
    .sis_pal_seq() +
    ggplot2::labs(
      title    = "Immune State Gradient Map",
      subtitle = "High values indicate immune transition zones",
      x = NULL, y = NULL
    ) +
    .sis_theme_map()
}


#' Plot composite immune scores (IAS / ISS / IBI)
#'
#' @param fit A fitted `sisHMRF` object.
#' @param coords N x 2 matrix of spatial coordinates.
#'
#' @return A `ggplot` object (3-panel via patchwork) or list.
#'
#' @export
sis_plot_scores <- function(fit, coords) {
  scores <- sis_scores(fit)
  df <- data.frame(
    x   = coords[, 1],
    y   = coords[, 2],
    IAS = scores[, "IAS"],
    ISS = scores[, "ISS"],
    IBI = scores[, "IBI"]
  )

  p1 <- ggplot2::ggplot(df, aes(x = x, y = y, color = IAS)) +
    ggplot2::geom_point(size = 1.2) +
    .sis_pal_div() +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = "IAS (Immune Activation)", x = NULL, y = NULL) +
    .sis_theme_map()

  p2 <- ggplot2::ggplot(df, aes(x = x, y = y, color = ISS)) +
    ggplot2::geom_point(size = 1.2) +
    .sis_pal_div() +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = "ISS (Immune Suppression)", x = NULL, y = NULL) +
    .sis_theme_map()

  p3 <- ggplot2::ggplot(df, aes(x = x, y = y, color = IBI)) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_color_gradientn(
      colors = c("#2166AC", "#F7F7F7", "#4DAF4A"),
      limits = c(0, 1), name = "IBI"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = "IBI (Immune Balance Index)", x = NULL, y = NULL) +
    .sis_theme_map()

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p1 + p2 + p3 + patchwork::plot_layout(ncol = 3)
  } else {
    list(IAS = p1, ISS = p2, IBI = p3)
  }
}


#' Plot ELBO convergence
#'
#' @param fit A fitted `sisHMRF` object.
#'
#' @return A `ggplot` object.
#'
#' @export
sis_plot_elbo <- function(fit) {
  df <- data.frame(
    iteration = seq_along(fit$elbo),
    elbo      = fit$elbo
  )

  ggplot2::ggplot(df, aes(x = iteration, y = elbo)) +
    ggplot2::geom_line(linewidth = 0.5, color = "#2166AC") +
    ggplot2::labs(title = "ELBO Convergence", x = "Iteration", y = "ELBO") +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(size = 7, face = "bold"),
      axis.title  = ggplot2::element_text(size = 6),
      axis.text   = ggplot2::element_text(size = 5)
    )
}


#' Plot feature loading heatmap
#'
#' @param fit A fitted `sisHMRF` object.
#' @param feature_names Optional character vector of feature names.
#'
#' @return A `ggplot` object.
#'
#' @export
sis_plot_loadings <- function(fit, feature_names = NULL) {
  W <- fit$W
  D <- nrow(W)
  P <- ncol(W)

  if (!is.null(feature_names)) {
    rn <- feature_names[seq_len(D)]
  } else {
    rn <- paste0("F", seq_len(D))
  }
  rownames(W) <- rn
  colnames(W) <- paste0("Dim", seq_len(P))

  W_long <- reshape2::melt(W, varnames = c("Feature", "Dimension"),
                           value.name = "Loading")

  ggplot2::ggplot(W_long, aes(x = Dimension, y = Feature, fill = Loading)) +
    ggplot2::geom_tile() +
    .sis_fill_div() +
    ggplot2::labs(title = "Feature Loadings (W matrix)") +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 7, face = "bold"),
      axis.text.x  = ggplot2::element_text(size = 6),
      axis.text.y  = ggplot2::element_text(size = 5),
      legend.title = ggplot2::element_text(size = 6),
      legend.text  = ggplot2::element_text(size = 5),
      legend.key.size = ggplot2::unit(0.3, "cm"),
      panel.grid   = ggplot2::element_blank()
    )
}


#' Plot posterior uncertainty (standard error)
#'
#' @param fit A fitted `sisHMRF` object.
#' @param coords N x 2 matrix of spatial coordinates.
#' @param dim_idx Which latent dimension (default: 1).
#'
#' @return A `ggplot` object.
#'
#' @export
sis_plot_uncertainty <- function(fit, coords, dim_idx = 1) {
  df <- data.frame(
    x  = coords[, 1],
    y  = coords[, 2],
    se = sqrt(fit$Z_var[, dim_idx])
  )

  ggplot2::ggplot(df, aes(x = x, y = y, color = se)) +
    ggplot2::geom_point(size = 1.5) +
    .sis_pal_seq() +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = paste("Posterior SE: Dimension", dim_idx),
      x = NULL, y = NULL
    ) +
    .sis_theme_map()
}


# ---- S3 plot method ----

#' Plot a fitted sisHMRF model
#'
#' @param x A fitted `sisHMRF` object.
#' @param coords Required. N x 2 matrix of spatial coordinates.
#' @param type Plot type: `"panel"` (all immune dimensions), `"gradient"`,
#'   `"elbo"`, `"scores"`, `"archetypes"`, `"loadings"`, `"uncertainty"`.
#' @param ... Additional arguments passed to the specific plot function.
#'
#' @return A `ggplot` object.
#'
#' @export
plot.sisHMRF <- function(x, coords = NULL, type = "panel", ...) {
  if (is.null(coords)) stop("coords must be provided for spatial plots")

  switch(
    type,
    panel       = sis_plot_panel(x, coords, ...),
    gradient    = sis_plot_gradient(x, coords),
    elbo        = sis_plot_elbo(x),
    scores      = sis_plot_scores(x, coords),
    archetypes  = sis_plot_archetypes(x, coords, ...),
    loadings    = sis_plot_loadings(x, ...),
    uncertainty = sis_plot_uncertainty(x, coords, ...),
    stop("Unknown plot type: ", type,
         ". Use panel, gradient, elbo, scores, archetypes, loadings, or uncertainty.")
  )
}
