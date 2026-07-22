#' Fit the sis-HMRF model
#'
#' @description
#' Fit a Spatial Immune State Hidden Markov Random Field (sis-HMRF) model
#' using variational expectation-maximization (EM). The model infers continuous
#' latent immune states \eqn{Z_i \in \mathbb{R}^P} for each spatial spot \eqn{i},
#' under a Gaussian Markov random field (GMRF) prior that encodes spatial
#' smoothness and a linear-Gaussian observation model.
#'
#' @section Model:
#' \deqn{F_i = W Z_i + b + \epsilon_i, \quad \epsilon_i \sim N(0, \sigma^2 I_D)}
#' \deqn{Z_i \mid Z_{-i} \sim N\left(\rho \frac{1}{d_i}\sum_{j\sim i} Z_j,\; \frac{\tau^2}{d_i} I_P\right)}
#'
#' @param F_mat Numeric matrix (N x D) of immune features. Rows = spots,
#'   columns = features. Must be normalized (centered/scaled recommended).
#' @param coords Numeric matrix (N x 2) of spatial coordinates.
#' @param P Integer. Number of latent immune state dimensions (default: 3).
#' @param k Integer. Number of spatial neighbors for graph construction
#'   (default: 6).
#' @param max_iter Integer. Maximum EM iterations (default: 500).
#' @param tol Numeric. Convergence tolerance on ELBO relative change
#'   (default: 1e-6).
#' @param verbose Logical. Print progress if TRUE (default: TRUE).
#' @param perturb_sd Numeric. Standard deviation of Gaussian noise added to
#'   PCA-initialized Z_mean and W. Use to escape local optima. Set to 0 to
#'   disable (default: 0). Recommended: 0.3 -- 0.7 for real data.
#'
#' @return An object of class `sisHMRF`, a list with components:
#'   \item{Z_mean}{N x P matrix of posterior mean latent states.}
#'   \item{Z_var}{N x P matrix of posterior variances (diagonal of covariance).}
#'   \item{W}{D x P loading matrix.}
#'   \item{b}{D-vector of feature intercepts.}
#'   \item{sigma2}{Scalar residual variance.}
#'   \item{rho}{Scalar spatial coupling parameter \eqn{\in [0.01, 0.99]}.}
#'   \item{tau2}{Scalar conditional spatial variance.}
#'   \item{A}{Sparse adjacency matrix (dgCMatrix).}
#'   \item{degrees}{Integer vector of node degrees.}
#'   \item{elbo}{Numeric vector of ELBO values across iterations.}
#'   \item{P}{Integer. Number of latent dimensions.}
#'   \item{converged}{Logical. Whether EM converged before max_iter.}
#'
#' @examples
#' \dontrun{
#' # Load Visium data (see vignette)
#' features <- sis_features(count_matrix, coords)
#' fit <- sis_fit(features$F_mat, features$coords, P = 3, k = 6)
#' print(fit)
#' plot(fit, features$coords, type = "panel")
#' }
#'
#' @references
#' Besag, J. (1974). Spatial interaction and the statistical analysis of
#' lattice systems. \emph{J. R. Stat. Soc. B}, 36(2), 192--236.
#'
#' Rue, H. & Held, L. (2005). \emph{Gaussian Markov Random Fields:
#' Theory and Applications}. Chapman & Hall/CRC.
#'
#' Blei, D.M., Kucukelbir, A. & McAuliffe, J.D. (2017). Variational
#' inference: A review for statisticians. \emph{JASA}, 112(518), 859--877.
#'
#' @importFrom Matrix sparseMatrix rowSums
#' @importFrom FNN get.knn
#' @importFrom stats sd
#' @export
sis_fit <- function(
    F_mat, coords,
    P          = 3,
    k          = 6,
    max_iter   = 500,
    tol        = 1e-6,
    verbose    = TRUE,
    perturb_sd = 0
) {
  # ---- Input validation ----
  if (!is.matrix(F_mat)) stop("F_mat must be a matrix")
  if (!is.matrix(coords)) stop("coords must be a matrix")
  if (nrow(F_mat) != nrow(coords)) {
    stop(sprintf(
      "F_mat (%d rows) and coords (%d rows) must have same number of spots",
      nrow(F_mat), nrow(coords)
    ))
  }
  if (P < 1 || P > ncol(F_mat)) {
    stop(sprintf("P = %d must be between 1 and %d", P, ncol(F_mat)))
  }

  N <- nrow(F_mat)
  D <- ncol(F_mat)

  # ---- Spatial graph ----
  if (verbose) cat("Building spatial graph (k =", k, ")...\n")
  A <- build_spatial_graph(coords, method = "knn", k = k)
  degrees <- as.integer(Matrix::rowSums(A))
  n_edges  <- sum(A > 0) / 2
  if (verbose) cat("  Graph:", N, "nodes,", n_edges, "edges\n")

  # ---- Initialize ----
  if (verbose) cat("Initializing parameters via PCA...\n")
  init <- .init_params(F_mat, P)

  Z_mean <- init$Z_mean
  Z_var  <- init$Z_var
  W      <- init$W
  b      <- init$b
  sigma2 <- init$sigma2
  rho    <- init$rho
  tau2   <- init$tau2

  # Perturb initialization to escape local optima
  if (perturb_sd > 0) {
    if (verbose) cat("Applying perturbed initialization (sd =", perturb_sd, ")...\n")
    set.seed(49)  # reproducible perturbation
    Z_mean <- Z_mean + matrix(rnorm(N * P, 0, perturb_sd), N, P)
    W      <- W + matrix(rnorm(D * P, 0, perturb_sd), D, P)
  }

  # ---- EM loop ----
  elbo_hist <- numeric(max_iter)

  for (iter in seq_len(max_iter)) {
    # E-step
    estep <- .e_step_core(F_mat, Z_mean, Z_var, W, b, sigma2, rho, tau2, A, degrees)
    Z_mean <- estep$Z_mean
    Z_var  <- estep$Z_var

    # M-step
    mstep <- .m_step_core(F_mat, Z_mean, Z_var, A, degrees, W, b, sigma2, rho, tau2)
    W      <- mstep$W
    b      <- mstep$b
    sigma2 <- mstep$sigma2
    rho    <- mstep$rho
    tau2   <- mstep$tau2

    # ELBO
    elbo_hist[iter] <- .compute_elbo(F_mat, Z_mean, Z_var, W, b, sigma2, rho, tau2, A, degrees)

    if (verbose && (iter %% 10 == 0 || iter == 1)) {
      cat(sprintf(
        "  Iter %3d: ELBO = %.4f, sigma2 = %.4f, rho = %.4f, tau2 = %.4f\n",
        iter, elbo_hist[iter], sigma2, rho, tau2
      ))
    }

    # Convergence
    if (iter > 10) {
      delta <- abs(elbo_hist[iter] - elbo_hist[iter - 1]) /
               (abs(elbo_hist[iter - 1]) + 1e-10)
      if (delta < tol) {
        if (verbose) cat("  Converged at iteration", iter, "\n")
        break
      }
    }
  }

  structure(
    list(
      Z_mean    = Z_mean,
      Z_var     = Z_var,
      W         = W,
      b         = b,
      sigma2    = sigma2,
      rho       = rho,
      tau2      = tau2,
      A         = A,
      degrees   = degrees,
      elbo      = elbo_hist[seq_len(iter)],
      P         = P,
      converged = (iter < max_iter),
      n_iter    = iter,
      N         = N,
      D         = D
    ),
    class = "sisHMRF"
  )
}


# ---- Internal: parameter initialization ----
#' @keywords internal
#' @noRd
.init_params <- function(F_mat, P) {
  N <- nrow(F_mat)
  D <- ncol(F_mat)

  F_scaled <- scale(F_mat, center = TRUE, scale = TRUE)
  b_init   <- attr(F_scaled, "scaled:center")

  svd_res  <- svd(F_scaled)
  Z_init   <- svd_res$u[, seq_len(P), drop = FALSE] %*%
              diag(svd_res$d[seq_len(P)], nrow = P)
  W_init   <- svd_res$v[, seq_len(P), drop = FALSE]

  # Normalize Z
  Z_sd <- apply(Z_init, 2, stats::sd)
  Z_sd[Z_sd < 1e-8] <- 1.0
  Z_init <- sweep(Z_init, 2, Z_sd, "/")
  W_init <- W_init %*% diag(Z_sd, nrow = P)

  list(
    Z_mean = as.matrix(Z_init),
    Z_var  = matrix(1.0, nrow = N, ncol = P),
    W      = W_init,
    b      = as.vector(b_init),
    sigma2 = 1.0,
    rho    = 0.5,
    tau2   = 1.0
  )
}


# ---- Internal: E-step ----
#' @keywords internal
#' @noRd
.e_step_core <- function(F_mat, Z_mean, Z_var, W, b, sigma2, rho, tau2, A, degrees) {
  N <- nrow(F_mat)
  P <- ncol(Z_mean)

  Z_new_mean <- matrix(0, N, P)
  Z_new_var  <- matrix(0, N, P)

  WT_W_s2 <- crossprod(W) / sigma2

  for (i in seq_len(N)) {
    prec <- WT_W_s2
    di <- degrees[i]
    if (di > 0) {
      prec <- prec + diag(di / tau2, P)
    }
    Gamma_i <- solve(prec)

    data_term <- crossprod(W, F_mat[i, ] - b) / sigma2

    spatial_term <- rep(0, P)
    if (di > 0) {
      nb <- which(A[i, ] > 0)
      spatial_term <- (rho / tau2) * colSums(Z_mean[nb, , drop = FALSE])
    }

    Z_new_mean[i, ] <- as.vector(Gamma_i %*% (data_term + spatial_term))
    Z_new_var[i, ]  <- diag(Gamma_i)
  }

  list(Z_mean = Z_new_mean, Z_var = Z_new_var)
}


# ---- Internal: M-step ----
#' @keywords internal
#' @noRd
.m_step_core <- function(F_mat, Z_mean, Z_var, A, degrees,
                         W_old, b_old, sigma2_old, rho_old, tau2_old) {
  N <- nrow(F_mat)
  D <- ncol(F_mat)
  P <- ncol(Z_mean)

  # E[Z Z^T]
  E_ZZ <- matrix(0, P, P)
  for (i in seq_len(N)) {
    E_ZZ <- E_ZZ + tcrossprod(Z_mean[i, ]) + diag(Z_var[i, ], P)
  }

  sum_F  <- colSums(F_mat)
  sum_FZ <- crossprod(F_mat, Z_mean)

  # W
  F_centered_Z <- crossprod(sweep(F_mat, 2, b_old, "-"), Z_mean)
  W_new <- F_centered_Z %*% solve(E_ZZ)

  # b
  b_new <- as.vector((sum_F - colSums(Z_mean %*% t(W_new))) / N)

  # sigma2
  resid_sq <- 0
  for (i in seq_len(N)) {
    pred <- as.vector(W_new %*% Z_mean[i, ] + b_new)
    tr   <- sum(diag(crossprod(W_new) %*% diag(Z_var[i, ], P)))
    resid_sq <- resid_sq + sum((F_mat[i, ] - pred)^2) + tr
  }
  sigma2_new <- max(resid_sq / (N * D), 1e-6)

  # rho
  spatial_num <- 0
  spatial_den <- 0
  for (i in seq_len(N)) {
    nb <- which(A[i, ] > 0)
    for (j in nb) {
      if (i < j) {
        spatial_num <- spatial_num + sum(Z_mean[i, ] * Z_mean[j, ])
        spatial_den <- spatial_den + sum(Z_mean[i, ] * Z_mean[i, ])
      }
    }
  }
  rho_new <- if (spatial_den > 0) spatial_num / spatial_den else rho_old
  rho_new <- pmin(pmax(rho_new, 0.01), 0.99)

  # tau2
  tau2_num <- 0
  n_edges <- 0
  for (i in seq_len(N)) {
    nb <- which(A[i, ] > 0)
    for (j in nb) {
      if (i < j) {
        diff_vec <- Z_mean[i, ] - rho_new * Z_mean[j, ]
        var_contrib <- Z_var[i, ] + rho_new^2 * Z_var[j, ]
        tau2_num <- tau2_num + sum(diff_vec^2) + sum(var_contrib)
        n_edges <- n_edges + 1
      }
    }
  }
  tau2_new <- max(tau2_num / (P * n_edges), 1e-4)

  list(
    W      = W_new,
    b      = b_new,
    sigma2 = sigma2_new,
    rho    = rho_new,
    tau2   = tau2_new
  )
}


# ---- Internal: ELBO computation ----
#' @keywords internal
#' @noRd
.compute_elbo <- function(F_mat, Z_mean, Z_var, W, b, sigma2, rho, tau2, A, degrees) {
  N <- nrow(F_mat)
  D <- ncol(F_mat)
  P <- ncol(Z_mean)

  # Expected log-likelihood
  ll <- 0
  for (i in seq_len(N)) {
    pred <- as.vector(W %*% Z_mean[i, ] + b)
    resid <- F_mat[i, ] - pred
    tr <- sum(diag(crossprod(W) %*% diag(Z_var[i, ], P)))
    ll <- ll - 0.5 * (sum(resid^2) + tr) / sigma2
  }
  ll <- ll - 0.5 * N * D * log(2 * pi * sigma2)

  # Expected log-prior
  prior <- 0
  for (i in seq_len(N)) {
    nb <- which(A[i, ] > 0)
    for (j in nb) {
      if (i < j) {
        diff_mean <- Z_mean[i, ] - rho * Z_mean[j, ]
        prior <- prior - 0.5 * sum(diff_mean^2) / tau2
      }
    }
  }
  prior <- prior - 0.5 * N * P * log(2 * pi * tau2)

  # Entropy
  entropy <- 0.5 * sum(log(2 * pi * exp(1) * Z_var))

  ll + prior + entropy
}


# ---- S3 methods ----

#' @export
print.sisHMRF <- function(x, ...) {
  cat("sisHMRF model\n")
  cat("==============\n")
  cat("Spots (N):         ", x$N, "\n")
  cat("Features (D):      ", x$D, "\n")
  cat("Latent dims (P):   ", x$P, "\n")
  cat("Converged:         ", x$converged, "(after", x$n_iter, "iterations)\n")
  cat("Final ELBO:        ", sprintf("%.4f", utils::tail(x$elbo, 1)), "\n")
  cat("sigma2:            ", sprintf("%.4f", x$sigma2), "\n")
  cat("rho (spatial):     ", sprintf("%.4f", x$rho), "\n")
  cat("tau2 (spatial var):", sprintf("%.4f", x$tau2), "\n")
  invisible(x)
}

#' @export
summary.sisHMRF <- function(object, ...) {
  fit <- object
  cat("sisHMRF Model Summary\n")
  cat("======================\n\n")

  cat(sprintf("Spots: %d | Features: %d | Latent dims: P = %d\n",
              fit$N, fit$D, fit$P))
  cat(sprintf("Converged: %s (iter %d / max 500)\n", fit$converged, fit$n_iter))
  cat(sprintf("Final ELBO: %.4f\n\n", utils::tail(fit$elbo, 1)))

  cat("Model parameters:\n")
  cat(sprintf("  sigma2 (residual variance):   %.4f\n", fit$sigma2))
  cat(sprintf("  rho    (spatial coupling):     %.4f\n", fit$rho))
  cat(sprintf("  tau2   (spatial cond. var.):   %.4f\n\n", fit$tau2))

  cat("Latent state summary:\n")
  cat(sprintf("  Z_min:  %.4f\n", min(fit$Z_mean)))
  cat(sprintf("  Z_max:  %.4f\n", max(fit$Z_mean)))
  cat(sprintf("  Z_sd:   %.4f\n", mean(apply(fit$Z_mean, 2, sd))))

  invisible(object)
}

#' @export
coef.sisHMRF <- function(object, ...) {
  object$W
}

#' @export
fitted.sisHMRF <- function(object, ...) {
  sweep(object$Z_mean %*% t(object$W), 2, object$b, "+")
}
