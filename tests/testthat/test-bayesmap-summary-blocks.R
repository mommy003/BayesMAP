test_that("one-block summary sampler matches dense summary sampler", {
  set.seed(123)

  N <- 80L
  m <- 20L
  G <- 6L
  C <- 3L

  X <- matrix(
    rnorm(N * m),
    nrow = N,
    ncol = m
  )

  X <- scale(X)

  beta_true <- numeric(m)
  beta_true[c(2, 7, 14)] <- c(0.20, -0.15, 0.18)

  y <- as.numeric(
    X %*% beta_true +
      rnorm(N)
  )

  y <- as.numeric(scale(y))

  bhat <- as.numeric(
    crossprod(X, y) / N
  )

  LD <- crossprod(X) / N
  LD <- (LD + t(LD)) / 2

  phenotype_variance <- as.numeric(
    crossprod(y) / N
  )

  L <- matrix(
    rbinom(m * G, 1L, 0.15),
    nrow = m,
    ncol = G
  )

  A <- matrix(
    rbinom(m * C, 1L, 0.20),
    nrow = m,
    ncol = C
  )

  B <- matrix(
    rbinom(G * C, 1L, 0.25),
    nrow = G,
    ncol = C
  )

  set.seed(999)

  fit_dense <- bayesmap_summary_R(
    bhat = bhat,
    LD = LD,
    N = N,
    L = L,
    A = A,
    B = B,
    phenotype_variance = phenotype_variance,
    niter = 60,
    burnin = 20,
    thin = 2,
    store_beta = FALSE,
    store_delta = FALSE,
    verbose = FALSE
  )

  set.seed(999)

  fit_block <- bayesmap_summary_blocks_R(
    bhat = bhat,
    LD_blocks = list(LD),
    block_indices = list(seq_len(m)),
    N = N,
    L = L,
    A = A,
    B = B,
    phenotype_variance = phenotype_variance,
    niter = 60,
    burnin = 20,
    thin = 2,
    store_beta = FALSE,
    store_delta = FALSE,
    verbose = FALSE
  )

  expect_s3_class(
    fit_block,
    "BayesMAPSummaryBlocks"
  )

  expect_equal(
    fit_block$n_blocks,
    1L
  )

  expect_equal(
    fit_block$block_sizes,
    m
  )

  expect_equal(
    fit_block$par,
    fit_dense$par,
    tolerance = 1e-12
  )

  expect_equal(
    fit_block$pip_snp,
    fit_dense$pip_snp,
    tolerance = 1e-12
  )

  expect_equal(
    fit_block$pip_gene,
    fit_dense$pip_gene,
    tolerance = 1e-12
  )

  expect_equal(
    fit_block$pip_cell,
    fit_dense$pip_cell,
    tolerance = 1e-12
  )
})


test_that("two-block summary sampler matches block-diagonal dense LD", {
  set.seed(456)

  N <- 100L
  m1 <- 10L
  m2 <- 12L
  m <- m1 + m2
  G <- 6L
  C <- 3L

  X1 <- matrix(
    rnorm(N * m1),
    nrow = N,
    ncol = m1
  )

  X2 <- matrix(
    rnorm(N * m2),
    nrow = N,
    ncol = m2
  )

  X1 <- scale(X1)
  X2 <- scale(X2)

  y <- as.numeric(scale(rnorm(N)))

  bhat1 <- as.numeric(
    crossprod(X1, y) / N
  )

  bhat2 <- as.numeric(
    crossprod(X2, y) / N
  )

  bhat <- c(bhat1, bhat2)

  LD1 <- crossprod(X1) / N
  LD1 <- (LD1 + t(LD1)) / 2

  LD2 <- crossprod(X2) / N
  LD2 <- (LD2 + t(LD2)) / 2

  index1 <- seq_len(m1)
  index2 <- m1 + seq_len(m2)

  LD_full <- matrix(
    0,
    nrow = m,
    ncol = m
  )

  LD_full[index1, index1] <- LD1
  LD_full[index2, index2] <- LD2

  phenotype_variance <- as.numeric(
    crossprod(y) / N
  )

  L <- matrix(
    rbinom(m * G, 1L, 0.15),
    nrow = m,
    ncol = G
  )

  A <- matrix(
    rbinom(m * C, 1L, 0.20),
    nrow = m,
    ncol = C
  )

  B <- matrix(
    rbinom(G * C, 1L, 0.25),
    nrow = G,
    ncol = C
  )

  set.seed(888)

  fit_dense <- bayesmap_summary_R(
    bhat = bhat,
    LD = LD_full,
    N = N,
    L = L,
    A = A,
    B = B,
    phenotype_variance = phenotype_variance,
    niter = 60,
    burnin = 20,
    thin = 2,
    store_beta = FALSE,
    store_delta = FALSE,
    verbose = FALSE
  )

  set.seed(888)

  fit_blocks <- bayesmap_summary_blocks_R(
    bhat = bhat,
    LD_blocks = list(LD1, LD2),
    block_indices = list(index1, index2),
    N = N,
    L = L,
    A = A,
    B = B,
    phenotype_variance = phenotype_variance,
    niter = 60,
    burnin = 20,
    thin = 2,
    store_beta = FALSE,
    store_delta = FALSE,
    verbose = FALSE
  )

  expect_equal(
    fit_blocks$n_blocks,
    2L
  )

  expect_equal(
    fit_blocks$block_sizes,
    c(m1, m2)
  )

  expect_equal(
    fit_blocks$par,
    fit_dense$par,
    tolerance = 1e-12
  )

  expect_equal(
    fit_blocks$pip_snp,
    fit_dense$pip_snp,
    tolerance = 1e-12
  )

  expect_equal(
    fit_blocks$pip_gene,
    fit_dense$pip_gene,
    tolerance = 1e-12
  )

  expect_equal(
    fit_blocks$pip_cell,
    fit_dense$pip_cell,
    tolerance = 1e-12
  )

  expect_true(
    all(is.finite(fit_blocks$par))
  )

  expect_true(
    all(fit_blocks$pip_snp >= 0 &
          fit_blocks$pip_snp <= 1)
  )
})
