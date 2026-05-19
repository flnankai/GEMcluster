test_that("gem_fit returns a valid fitted object", {
  skip_if_not_installed("mvtnorm")
  skip_if_not_installed("huge")

  set.seed(2)
  Mu <- dense_block_means(K = 3, p = 6, signal = 0.8)
  Sigma <- make_ar_cov(p = 6, rho = 0.4)
  dat <- simulate_elliptical_mixture(
    n = 36,
    Mu = Mu,
    Sigma = Sigma,
    family = "gaussian",
    seed = 2
  )

  fit <- gem_fit(
    dat$X,
    K = 3,
    max_outer = 1,
    init_n_boot = 1,
    init_nstart = 2,
    init_max_iter = 5,
    outer_nstart = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "gem_fit")
  expect_equal(length(fit$cluster), nrow(dat$X))
  expect_equal(dim(fit$Tau), c(nrow(dat$X), 3L))
  expect_equal(length(predict(fit)), nrow(dat$X))
  expect_equal(dim(gem_predict(fit, type = "posterior")), c(nrow(dat$X), 3L))
})

test_that("gem_select_k_gap returns K-selection summary", {
  skip_if_not_installed("mvtnorm")
  skip_if_not_installed("huge")

  set.seed(3)
  Mu <- dense_block_means(K = 3, p = 5, signal = 0.9)
  Sigma <- make_ar_cov(p = 5, rho = 0.3)
  dat <- simulate_elliptical_mixture(
    n = 30,
    Mu = Mu,
    Sigma = Sigma,
    family = "gaussian",
    seed = 3
  )

  gap <- gem_select_k_gap(
    dat$X,
    K_grid = 2:3,
    B = 1,
    control = list(
      max_outer = 1,
      init_n_boot = 1,
      init_nstart = 2,
      init_max_iter = 5,
      outer_nstart = 1,
      verbose = FALSE
    ),
    verbose = FALSE
  )

  expect_s3_class(gap, "gem_gap")
  expect_true(gap$best_K %in% 2:3)
  expect_equal(nrow(gap$summary), 2L)
})
