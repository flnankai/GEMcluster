test_that("simulation helpers return consistent dimensions", {
  skip_if_not_installed("mvtnorm")

  set.seed(1)
  Mu <- dense_block_means(K = 3, p = 6, signal = 0.5)
  Sigma <- make_ar_cov(p = 6, rho = 0.5)
  dat <- simulate_elliptical_mixture(
    n = 30,
    Mu = Mu,
    Sigma = Sigma,
    family = "gaussian",
    seed = 1
  )

  expect_equal(dim(dat$X), c(30L, 6L))
  expect_equal(length(dat$label), 30L)
  expect_equal(sort(unique(dat$label)), 1:3)
  expect_equal(dim(Sigma), c(6L, 6L))
})

test_that("cluster accuracy is invariant to label permutation", {
  truth <- c(1, 1, 2, 2, 3, 3)
  pred <- c(2, 2, 3, 3, 1, 1)
  expect_equal(cluster_accuracy(pred, truth)$accuracy, 1)
})
