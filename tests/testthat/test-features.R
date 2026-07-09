test_that("get_immune_gene_sets returns valid structure", {
  gs <- get_immune_gene_sets()

  expect_type(gs, "list")
  expect_gt(length(gs), 20)
  expect_named(gs)

  # Every set should be a character vector of gene symbols
  for (s in names(gs)) {
    expect_type(gs[[s]], "character")
    expect_gt(length(gs[[s]]), 1)
  }
})

test_that("sis_features runs on simulated data", {
  set.seed(49)
  N <- 100
  G <- 500

  # Include some immune genes in the data
  immune_genes <- c("CD8A", "GZMB", "PRF1", "CD3D", "CD274", "FOXP3",
                    "CTLA4", "CD4", "NKG7", "IL10", "TGFB1", "CXCL9",
                    "CXCL10", "HLA-A", "HLA-B", "PDCD1")
  all_genes <- c(immune_genes, paste0("Random", seq_len(G - length(immune_genes))))

  expr_mat <- matrix(rnorm(N * G), N, G,
                     dimnames = list(NULL, all_genes))

  features <- sis_features(expr_mat, verbose = FALSE)

  expect_s3_class(features, "sisFeatures")
  expect_true(is.matrix(features$F_mat))
  expect_equal(nrow(features$F_mat), N)
  expect_gt(features$n_features, 10)
  expect_output(print(features), "Immune Feature Matrix")
})

test_that("sis_scores computes valid scores", {
  set.seed(49)
  N <- 60
  D <- 15
  coords <- cbind(runif(N), runif(N))
  F_mat <- matrix(rnorm(N * D), N, D)

  fit <- sis_fit(F_mat, coords, P = 3, k = 5, verbose = FALSE)
  scores <- sis_scores(fit)

  expect_equal(dim(scores), c(N, 3))
  expect_true("IAS" %in% colnames(scores))
  expect_true("ISS" %in% colnames(scores))
  expect_true("IBI" %in% colnames(scores))

  # IBI should be clamped to [0, 1]
  expect_true(all(scores[, "IBI"] >= 0 & scores[, "IBI"] <= 1))
})

test_that("compute_immune_gradient is non-negative", {
  set.seed(49)
  N <- 60
  D <- 15
  coords <- cbind(runif(N), runif(N))
  F_mat <- matrix(rnorm(N * D), N, D)

  fit <- sis_fit(F_mat, coords, P = 2, k = 5, verbose = FALSE)
  grad <- compute_immune_gradient(fit)

  expect_length(grad, N)
  expect_true(all(grad >= 0))
})

test_that("cluster_immune_states returns valid clusters", {
  set.seed(49)
  N <- 80
  D <- 15
  coords <- cbind(runif(N), runif(N))
  F_mat <- matrix(rnorm(N * D), N, D)

  fit <- sis_fit(F_mat, coords, P = 3, k = 5, verbose = FALSE)

  for (K in c(2, 4, 6)) {
    cl <- cluster_immune_states(fit, K = K)
    expect_s3_class(cl, "factor")
    expect_length(cl, N)
    expect_equal(nlevels(cl), K)
  }
})

test_that("build_spatial_graph produces valid adjacency", {
  N <- 50
  coords <- cbind(runif(N), runif(N))

  A <- build_spatial_graph(coords, method = "knn", k = 5)

  expect_s4_class(A, "dgCMatrix")
  expect_equal(dim(A), c(N, N))
  expect_true(Matrix::isSymmetric(A))
  expect_equal(sum(Matrix::diag(A)), 0)  # no self-loops
})
