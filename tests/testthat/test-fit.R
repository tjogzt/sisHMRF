test_that("sis_fit runs on simulated data", {
  set.seed(49)
  N <- 100
  D <- 20
  P <- 2

  # Simulate data with spatial structure
  coords <- cbind(x = runif(N), y = runif(N))
  Z_true <- cbind(
    sin(coords[, 1] * pi),
    cos(coords[, 2] * pi)
  )
  W_true <- matrix(rnorm(D * P, sd = 0.5), D, P)
  b_true <- rnorm(D, sd = 0.2)

  F_mat <- t(W_true %*% t(Z_true) + b_true) +
            matrix(rnorm(N * D, sd = 0.3), N, D)
  colnames(F_mat) <- paste0("F", seq_len(D))

  fit <- sis_fit(F_mat, coords, P = P, k = 5, verbose = FALSE)

  expect_s3_class(fit, "sisHMRF")
  expect_equal(fit$P, P)
  expect_equal(nrow(fit$Z_mean), N)
  expect_equal(ncol(fit$Z_mean), P)
  expect_equal(nrow(fit$W), D)
  expect_equal(ncol(fit$W), P)
  expect_true(fit$converged)
  expect_gt(fit$sigma2, 0)
  expect_true(fit$rho >= 0.01 && fit$rho <= 0.99)
  expect_gt(fit$tau2, 0)
  expect_length(fit$elbo, fit$n_iter)
})

test_that("sis_fit handles edge cases", {
  N <- 50
  D <- 10
  coords <- cbind(seq_len(N), rep(1, N))
  F_mat <- matrix(rnorm(N * D), N, D)

  # P = 1
  fit1 <- sis_fit(F_mat, coords, P = 1, k = 4, verbose = FALSE)
  expect_s3_class(fit1, "sisHMRF")
  expect_equal(ncol(fit1$Z_mean), 1)

  # P must be <= D
  expect_error(
    sis_fit(F_mat, coords, P = D + 1, verbose = FALSE),
    "must be between"
  )

  # Dimension mismatch
  expect_error(
    sis_fit(F_mat[1:10, ], coords, verbose = FALSE),
    "same number of spots"
  )
})

test_that("S3 methods work", {
  set.seed(49)
  N <- 60
  D <- 15
  coords <- cbind(runif(N), runif(N))
  F_mat <- matrix(rnorm(N * D), N, D)

  fit <- sis_fit(F_mat, coords, P = 2, k = 5, verbose = FALSE)

  expect_output(print(fit), "sisHMRF model")
  expect_output(summary(fit), "Model Summary")

  W <- coef(fit)
  expect_equal(dim(W), c(D, 2))

  fv <- fitted(fit)
  expect_equal(dim(fv), c(N, D))
})

test_that("predict works for new data", {
  set.seed(49)
  N <- 80
  D <- 12
  P <- 2

  # Generate data with actual factor structure so prediction is meaningful
  coords <- cbind(runif(N), runif(N))
  Z_true <- cbind(sin(coords[, 1] * pi), cos(coords[, 2] * pi))
  W_true <- matrix(rnorm(D * P, sd = 0.5), D, P)
  b_true <- rnorm(D, sd = 0.2)
  F_mat <- t(W_true %*% t(Z_true) + b_true) +
            matrix(rnorm(N * D, sd = 0.3), N, D)

  fit <- sis_fit(F_mat, coords, P = P, k = 5, verbose = FALSE)

  # Predict on same data (factor-model only, no spatial smoothing)
  pred <- predict(fit, F_mat, coords)
  expect_equal(dim(pred), c(N, P))
  # Predictions should be finite and well-behaved
  expect_true(all(is.finite(pred)))

  # Predict on new data without coords
  new_F <- matrix(rnorm(20 * D), 20, D)
  pred_new <- predict(fit, new_F)
  expect_equal(dim(pred_new), c(20, P))
})
