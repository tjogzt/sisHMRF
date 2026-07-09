#' Compute immune scores from fitted sis-HMRF model
#'
#' @description
#' Computes per-spot composite immune scores from the inferred latent states.
#' Three scores are computed:
#' \describe{
#'   \item{IAS (Immune Activation Score)}{Sum of latent dimensions with
#'     positive mean feature loadings, reflecting immune activation.}
#'   \item{ISS (Immune Suppression Score)}{Sum of absolute values of latent
#'     dimensions with negative mean feature loadings, reflecting
#'     immunosuppression.}
#'   \item{IBI (Immune Balance Index)}{IAS / (IAS + ISS), ranging from 0
#'     (pure suppression) to 1 (pure activation).}
#' }
#'
#' @param fit A fitted `sisHMRF` object.
#' @param feature_names Optional character vector of feature names corresponding
#'   to rows of W (used to annotate loadings).
#'
#' @return A N x 3 matrix with columns IAS, ISS, IBI.
#'
#' @examples
#' \dontrun{
#' fit <- sis_fit(F_mat, coords, P = 3)
#' scores <- sis_scores(fit)
#' head(scores)
#' }
#'
#' @export
sis_scores <- function(fit, feature_names = NULL) {
  if (!inherits(fit, "sisHMRF")) stop("fit must be a sisHMRF object")

  Z <- fit$Z_mean
  W <- fit$W
  P <- ncol(Z)

  if (!is.null(feature_names)) {
    rownames(W) <- feature_names[seq_len(nrow(W))]
  }

  # Determine activation vs suppression dimensions
  mean_loading <- colMeans(W)
  act_dim <- which(mean_loading > 0)
  sup_dim <- which(mean_loading < 0)

  IAS <- if (length(act_dim) > 0) {
    rowSums(Z[, act_dim, drop = FALSE])
  } else {
    numeric(nrow(Z))
  }

  ISS <- if (length(sup_dim) > 0) {
    rowSums(abs(Z[, sup_dim, drop = FALSE]))
  } else {
    numeric(nrow(Z))
  }

  IBI <- IAS / (IAS + ISS + 1e-8)
  # Clamp to [0, 1] since IAS can be negative with random data
  IBI <- pmin(pmax(IBI, 0), 1)

  cbind(IAS = IAS, ISS = ISS, IBI = IBI)
}


#' Compute per-spot immune gradient magnitude
#'
#' @description
#' The immune gradient at spot \eqn{i} is the root-mean-square difference
#' between \eqn{Z_i} and the mean of its neighbors' \eqn{Z} values:
#' \deqn{g_i = \sqrt{\frac{1}{d_i}\sum_{j\sim i} \|Z_i - Z_j\|^2}}
#' High gradient values indicate immune transition zones (boundaries between
#' distinct immune microenvironments).
#'
#' @param fit A fitted `sisHMRF` object.
#'
#' @return Numeric vector of length N (one gradient value per spot).
#'
#' @export
compute_immune_gradient <- function(fit) {
  if (!inherits(fit, "sisHMRF")) stop("fit must be a sisHMRF object")

  N    <- nrow(fit$Z_mean)
  grad <- numeric(N)

  for (i in seq_len(N)) {
    nb <- which(fit$A[i, ] > 0)
    if (length(nb) > 0) {
      diffs <- sweep(fit$Z_mean[nb, , drop = FALSE], 2, fit$Z_mean[i, ], "-")
      grad[i] <- sqrt(mean(rowSums(diffs^2)))
    }
  }
  grad
}


#' Post-hoc clustering of continuous immune states
#'
#' @description
#' Applies hierarchical clustering (Ward's D2) to the latent state matrix
#' to derive discrete immune archetypes for interpretation.
#'
#' @param fit A fitted `sisHMRF` object.
#' @param K Integer. Number of immune archetypes (default: 4).
#'
#' @return Factor of length N with cluster assignments (levels:
#'   "IS1", "IS2", ...).
#'
#' @export
cluster_immune_states <- function(fit, K = 4) {
  if (!inherits(fit, "sisHMRF")) stop("fit must be a sisHMRF object")

  Z_scaled <- scale(fit$Z_mean)
  hc <- hclust(dist(Z_scaled), method = "ward.D2")
  clusters <- cutree(hc, k = K)

  factor(clusters, levels = seq_len(K),
         labels = paste0("IS", seq_len(K)))
}


#' Model selection via BIC
#'
#' Compares sis-HMRF fits across a range of latent dimensions P using BIC.
#'
#' @param F_mat N x D immune feature matrix.
#' @param coords N x 2 spatial coordinates.
#' @param P_max Maximum P to evaluate (default: 8).
#' @param k Spatial neighbors for graph (default: 6).
#' @param verbose Print progress per P (default: TRUE).
#'
#' @return Data frame with columns P, BIC, ELBO.
#'
#' @export
sis_select_P <- function(F_mat, coords, P_max = 8, k = 6, verbose = TRUE) {
  N <- nrow(F_mat)
  D <- ncol(F_mat)
  results <- data.frame(P = integer(), BIC = numeric(), ELBO = numeric())

  for (p in seq_len(P_max)) {
    if (verbose) cat("Fitting P =", p, "...\n")
    fit <- sis_fit(F_mat, coords, P = p, k = k, verbose = FALSE)

    n_params <- D * p + D + 3  # W: D*p, b: D, sigma2, rho, tau2
    bic <- -2 * utils::tail(fit$elbo, 1) + n_params * log(N)

    results <- rbind(results, data.frame(
      P = p, BIC = bic, ELBO = utils::tail(fit$elbo, 1)
    ))
  }
  results
}
