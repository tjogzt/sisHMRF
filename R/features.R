#' Build immune feature matrix for sis-HMRF
#'
#' @description
#' Constructs the immune feature matrix \eqn{F \in \mathbb{R}^{N \times D}} from
#' spatial transcriptomics expression data. Features span seven modules of the
#' cancer-immunity cycle: innate immunity, interferon response, adaptive immunity,
#' immunosuppression, antigen presentation, chemokine signaling, and cell-type
#' composition. Features are computed via simplified centroid scoring (mean
#' expression of each gene set minus global background, z-scored across spots).
#'
#' @param expr_mat Numeric matrix (N x G). Log-normalized expression matrix.
#'   Rows = spots, columns = genes (gene symbols as column names).
#' @param include_pathways Logical. Include pathway-level features (default: TRUE).
#' @param include_cell_types Logical. Include cell-type scores (default: TRUE).
#' @param scale_features Logical. Scale features to unit variance (default: TRUE).
#' @param verbose Logical. Print progress (default: TRUE).
#'
#' @return A list with components:
#'   \item{F_mat}{N x D feature matrix (z-scored across features).}
#'   \item{feature_names}{Character vector of D feature names.}
#'   \item{feature_categories}{Character vector of categories ("Pathway" or "CellType").}
#'   \item{n_features}{Total number of features.}
#'
#' @examples
#' \dontrun{
#' feats <- sis_features(count_matrix)
#' head(feats$F_mat[, 1:5])
#' }
#'
#' @export
sis_features <- function(
    expr_mat,
    include_pathways   = TRUE,
    include_cell_types = TRUE,
    scale_features     = TRUE,
    verbose            = TRUE
) {
  if (!is.matrix(expr_mat)) {
    expr_mat <- as.matrix(expr_mat)
  }

  gene_sets <- get_immune_gene_sets()
  blocks    <- list()
  meta_cat  <- c()
  meta_name <- c()

  # ---- Pathway activity ----
  if (include_pathways) {
    pathway_names <- c(
      "Inflammasome", "NFkB_Signaling", "TNFa_Signaling", "TLR_Signaling",
      "IFNalpha_Response", "IFNgamma_Response",
      "T_Cell_Activation", "CD8_Cytotoxicity", "B_Cell_Plasma",
      "Treg_Signature", "Immune_Checkpoint",
      "TGFb_Signaling", "IL10_AntiInflammatory",
      "MHC_I_Presentation", "MHC_II_Presentation",
      "ProInflammatory_Chemokines", "T_Cell_Recruiting_Chemokines"
    )
    if (verbose) message("Computing pathway activity scores (", length(pathway_names), " pathways)...")
    pw_sets <- gene_sets[pathway_names]
    pw_scores <- .compute_pathway_scores(expr_mat, pw_sets)
    blocks[[length(blocks) + 1]] <- pw_scores
    meta_cat  <- c(meta_cat,  rep("Pathway", ncol(pw_scores)))
    meta_name <- c(meta_name, colnames(pw_scores))
  }

  # ---- Cell-type scores ----
  if (include_cell_types) {
    ct_names <- c(
      "CD8_T_Cell", "CD4_T_Cell", "Macrophage_M1", "Macrophage_M2",
      "NK_Cell", "Neutrophil", "Dendritic_Cell", "Endothelial",
      "Fibroblast"
    )
    if (verbose) message("Computing cell-type scores (", length(ct_names), " cell types)...")
    ct_sets <- gene_sets[ct_names]
    ct_scores <- .compute_ct_scores(expr_mat, ct_sets)
    blocks[[length(blocks) + 1]] <- ct_scores
    meta_cat  <- c(meta_cat,  rep("CellType", ncol(ct_scores)))
    meta_name <- c(meta_name, colnames(ct_scores))
  }

  if (length(blocks) == 0) stop("No features selected. Set include_pathways or include_cell_types to TRUE.")

  # ---- Combine ----
  F_mat <- do.call(cbind, blocks)

  # Remove zero-variance features
  fvars <- apply(F_mat, 2, var, na.rm = TRUE)
  keep  <- fvars > 1e-10
  F_mat     <- F_mat[, keep, drop = FALSE]
  meta_cat  <- meta_cat[keep]
  meta_name <- meta_name[keep]

  # Scale
  if (scale_features) {
    F_mat <- scale(F_mat, center = TRUE, scale = TRUE)
  }

  if (verbose) {
    message(sprintf("Built immune feature matrix: %d spots x %d features", nrow(F_mat), ncol(F_mat)))
  }

  structure(
    list(
      F_mat              = F_mat,
      feature_names      = meta_name,
      feature_categories = meta_cat,
      n_features         = ncol(F_mat)
    ),
    class = "sisFeatures"
  )
}


# ---- Internal: pathway scoring ----
#' @keywords internal
#' @noRd
.compute_pathway_scores <- function(expr_mat, gene_sets) {
  N <- nrow(expr_mat)
  S <- length(gene_sets)
  scores <- matrix(NA, nrow = N, ncol = S)
  colnames(scores) <- names(gene_sets)

  global_mean <- colMeans(expr_mat)

  for (s in seq_len(S)) {
    genes <- intersect(gene_sets[[s]], colnames(expr_mat))
    if (length(genes) < 3) {
      scores[, s] <- 0
      next
    }
    set_mean <- rowMeans(expr_mat[, genes, drop = FALSE])
    bg_mean  <- mean(global_mean[genes])
    scores[, s] <- set_mean - bg_mean
  }

  # z-score across pathways per spot
  scores <- t(scale(t(scores)))
  scores[is.na(scores)] <- 0
  scores
}


# ---- Internal: cell-type scoring ----
#' @keywords internal
#' @noRd
.compute_ct_scores <- function(expr_mat, ct_sets) {
  N <- nrow(expr_mat)
  C <- length(ct_sets)
  scores <- matrix(NA, nrow = N, ncol = C)
  colnames(scores) <- names(ct_sets)

  for (c in seq_len(C)) {
    genes <- intersect(ct_sets[[c]], colnames(expr_mat))
    if (length(genes) == 0) {
      scores[, c] <- 0
      next
    }
    scores[, c] <- rowMeans(expr_mat[, genes, drop = FALSE])
  }

  # Normalize per spot to approximate compositional constraints
  scores <- scores / (rowSums(abs(scores)) + 1e-8)
  scores
}


#' @export
print.sisFeatures <- function(x, ...) {
  cat("sisHMRF Immune Feature Matrix\n")
  cat("==============================\n")
  cat(sprintf("Spots:    %d\n", nrow(x$F_mat)))
  cat(sprintf("Features: %d\n", x$n_features))
  cat("\nFeature breakdown:\n")
  tab <- table(x$feature_categories)
  for (n in names(tab)) {
    cat(sprintf("  %s: %d\n", n, tab[n]))
  }
  invisible(x)
}
