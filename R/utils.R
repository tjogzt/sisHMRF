#' Spatial graph construction and utility functions
#'
#' @name sis_utils
NULL


#' Build spatial neighborhood graph
#'
#' Constructs a symmetric k-nearest-neighbor or radius-based adjacency
#' matrix from spatial coordinates.
#'
#' @param coords N x 2 matrix of spatial coordinates.
#' @param method Graph construction method: `"knn"` (default) or `"radius"`.
#' @param k Number of nearest neighbors (for `method = "knn"`).
#' @param radius Distance threshold (for `method = "radius"`).
#'
#' @return Sparse adjacency matrix (dgCMatrix, symmetric, binary).
#'
#' @importFrom FNN get.knn
#' @importFrom Matrix sparseMatrix
#'
#' @export
build_spatial_graph <- function(coords, method = "knn", k = 6, radius = NULL) {
  N <- nrow(coords)

  if (method == "knn") {
    knn_res <- FNN::get.knn(coords, k = k)
    i_idx <- rep(seq_len(N), each = k)
    j_idx <- as.vector(t(knn_res$nn.index))

    A <- Matrix::sparseMatrix(
      i = i_idx, j = j_idx, x = 1, dims = c(N, N)
    )
    A <- (A + Matrix::t(A)) > 0
    A <- as(as(A, "dMatrix"), "dgCMatrix")

  } else if (method == "radius") {
    if (is.null(radius)) stop("radius must be specified for radius method")
    dist_mat <- as.matrix(dist(coords))
    A <- dist_mat <= radius
    diag(A) <- 0L
    A <- as(A, "dgCMatrix")

  } else {
    stop("Unknown method: ", method)
  }

  diag(A) <- 0L
  A
}


#' Compute node degrees from adjacency matrix
#'
#' @param A Sparse adjacency matrix.
#' @return Integer vector of node degrees.
#'
#' @keywords internal
#' @export
compute_degrees <- function(A) {
  as.integer(Matrix::rowSums(A))
}


#' Get canonical immune gene signatures
#'
#' Returns a named list of gene sets covering major immune pathways and
#' cell types for feature construction. Includes 26 gene sets spanning seven
#' modules of the cancer-immunity cycle. All gene symbols are HGNC (human).
#'
#' @return Named list of gene sets (character vectors of gene symbols).
#'
#' @examples
#' gs <- get_immune_gene_sets()
#' names(gs)
#'
#' @export
get_immune_gene_sets <- function() {

  list(
    # --- Innate immunity ---
    "Inflammasome"  = c("NLRP3", "PYCARD", "CASP1", "IL1B", "IL18", "AIM2"),
    "NFkB_Signaling" = c("NFKB1", "NFKB2", "RELA", "RELB", "IKBKG",
                         "CHUK", "IKBKB", "NFKBIA"),
    "TNFa_Signaling" = c("TNF", "TNFRSF1A", "TRADD", "TRAF2", "RIPK1",
                          "TAB2", "TAB3", "MAP3K7"),
    "TLR_Signaling"  = c("TLR2", "TLR3", "TLR4", "TLR7", "TLR8", "TLR9",
                         "MYD88", "TIRAP", "TRIF", "TRAF6", "IRAK1", "IRAK4"),

    # --- Interferon response ---
    "IFNalpha_Response" = c("IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM2",
                            "IFITM3", "ISG15", "MX1", "MX2", "OAS1", "OAS2",
                            "OAS3", "OASL", "RSAD2", "STAT1", "STAT2", "IRF9"),
    "IFNgamma_Response" = c("STAT1", "IRF1", "IRF8", "GBP1", "GBP2", "GBP4",
                            "GBP5", "CXCL9", "CXCL10", "CXCL11", "IDO1",
                            "CIITA", "HLA-A", "HLA-B", "HLA-C", "B2M",
                            "TAP1", "TAP2", "PSMB8", "PSMB9"),

    # --- Adaptive immunity ---
    "T_Cell_Activation" = c("CD3D", "CD3E", "CD3G", "CD247", "ZAP70", "LCK",
                            "LAT", "ITK", "CD2", "CD28", "TRAC", "TRBC1"),
    "CD8_Cytotoxicity"  = c("CD8A", "CD8B", "GZMA", "GZMB", "GZMH", "GZMK",
                            "PRF1", "FASLG", "GNLY", "NKG7"),
    "B_Cell_Plasma"     = c("CD19", "CD79A", "CD79B", "MS4A1", "BLK", "PAX5",
                            "SDC1", "MZB1", "JCHAIN", "IGHG1", "IGHA1", "IGKC"),

    # --- Immunosuppression ---
    "Treg_Signature"     = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",
                             "TNFRSF18", "IL10", "EBi3", "LAG3", "ENTPD1"),
    "Immune_Checkpoint"  = c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3",
                             "HAVCR2", "TIGIT", "VSIR", "BTLA", "IDO1",
                             "CD276", "VTCN1"),
    "TGFb_Signaling"     = c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2",
                             "SMAD2", "SMAD3", "SMAD4", "SMAD7"),
    "IL10_AntiInflammatory" = c("IL10", "IL10RA", "IL10RB", "TGFB1",
                                "IL4", "IL13", "CSF2", "VEGFA"),

    # --- Antigen presentation ---
    "MHC_I_Presentation"  = c("HLA-A", "HLA-B", "HLA-C", "B2M",
                              "TAP1", "TAP2", "TAPBP", "CALR", "PDIA3"),
    "MHC_II_Presentation" = c("HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1",
                              "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA",
                              "CTSL", "CTSS"),

    # --- Chemokines ---
    "ProInflammatory_Chemokines"   = c("CCL2", "CCL3", "CCL4", "CCL5", "CCL7",
                                       "CCL8", "CCL11", "CCL13", "CXCL1", "CXCL2",
                                       "CXCL3", "CXCL5", "CXCL6", "CXCL8"),
    "T_Cell_Recruiting_Chemokines" = c("CXCL9", "CXCL10", "CXCL11", "CXCL13",
                                       "CCL19", "CCL21"),

    # --- Cell-type markers ---
    "CD8_T_Cell"      = c("CD8A", "CD8B", "GZMK", "CCL5", "NKG7", "CD3D"),
    "CD4_T_Cell"      = c("CD4", "IL7R", "CD40LG", "CCR7", "LEF1"),
    "Macrophage_M1"   = c("CD68", "IL1B", "TNF", "IL6", "CCL3", "IL12B",
                          "NOS2", "SOCS3"),
    "Macrophage_M2"   = c("CD163", "MSR1", "MRC1", "IL10", "TGFB1",
                          "CCL18", "ARG1", "CHI3L1"),
    "NK_Cell"         = c("NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1", "NCAM1"),
    "Neutrophil"      = c("FCGR3B", "CSF3R", "CXCR1", "CXCR2", "ELANE", "MPO"),
    "Dendritic_Cell"  = c("CLEC9A", "XCR1", "CLEC10A", "FCER1A", "CLEC4C",
                          "IRF7", "IRF8", "BATF3"),
    "Endothelial"     = c("PECAM1", "CDH5", "VWF", "ENG", "CLDN5", "ESAM"),
    "Fibroblast"      = c("COL1A1", "COL1A2", "DCN", "LUM", "FAP", "ACTA2",
                          "PDGFRA", "PDGFRB", "THY1")
  )
}


#' Predict latent states for new spots (out-of-sample)
#'
#' Given a fitted model and new feature data with optional spatial context,
#' predicts posterior latent states. For spots with no spatial neighbors,
#' falls back to the factor model with no spatial regularization.
#'
#' @param object A fitted `sisHMRF` object.
#' @param newdata N_new x D feature matrix for new spots.
#' @param new_coords Optional N_new x 2 coordinate matrix.
#' @param ... Not used.
#'
#' @return N_new x P matrix of predicted latent states.
#'
#' @export
predict.sisHMRF <- function(object, newdata, new_coords = NULL, ...) {
  if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
  N_new <- nrow(newdata)
  P     <- object$P

  Z_pred <- matrix(0, N_new, P)

  for (i in seq_len(N_new)) {
    # Factor-model prediction: Z = (W^T W)^{-1} W^T (F_i - b)
    WT_W_inv <- solve(crossprod(object$W))
    data_term <- crossprod(object$W, newdata[i, ] - object$b)
    Z_pred[i, ] <- as.vector(WT_W_inv %*% data_term)
  }

  # If coords provided, apply one smoothing pass using trained GMRF
  if (!is.null(new_coords) && is.matrix(new_coords)) {
    # Build graph among new spots
    A_new <- build_spatial_graph(new_coords, method = "knn", k = 6)
    degrees_new <- as.integer(Matrix::rowSums(A_new))

    rho  <- object$rho
    tau2 <- object$tau2

    Z_smooth <- Z_pred
    for (i in seq_len(N_new)) {
      di <- degrees_new[i]
      if (di > 0) {
        nb <- which(A_new[i, ] > 0)
        spatial_term <- (rho / tau2) * colSums(Z_smooth[nb, , drop = FALSE])
        # One-step smoothing (not full E-step)
        lambda <- di / (di + tau2)
        Z_smooth[i, ] <- (1 - lambda) * Z_pred[i, ] + lambda * spatial_term
      }
    }
    Z_pred <- Z_smooth
  }

  Z_pred
}
