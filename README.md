# GEMcluster

`GEMcluster` implements generalized expectation-maximization clustering for
high-dimensional semiparametric elliptical mixtures. The method is designed for
settings where component distributions may be heavy tailed, the radial
generator is not specified parametrically, and the common precision-shape
matrix must be regularized in high dimension.

The package follows the method in:

Feng, L. and Zhuang, D. (2026). *Semiparametric Elliptical Mixture Clustering
for High-Dimensional Data*. arXiv:2605.08995.
https://arxiv.org/abs/2605.08995

## Main functions

- `gem_fit()`: fit GEM clustering for a fixed number of clusters.
- `gem_predict()`: predict labels or posterior probabilities from a fitted
  model.
- `gem_select_k_gap()`: select the number of clusters with the GAP-LSE rule.
- `simulate_elliptical_mixture()`: generate Gaussian, t, Laplace, or slash
  elliptical mixture samples for examples and simulation checks.

## Installation

From a local source package:

```r
install.packages("GEMcluster_0.1.0.tar.gz", repos = NULL, type = "source")
```

From a local source directory:

```r
install.packages("GEMcluster", repos = NULL, type = "source")
```

The package uses Rcpp and RcppArmadillo. On Windows, Rtools must be available.

## Quick start

```r
library(GEMcluster)

set.seed(1)
K <- 3
p <- 20
Mu <- dense_block_means(K = K, p = p, signal = 0.8)
Sigma <- make_ar_cov(p = p, rho = 0.5)
dat <- simulate_elliptical_mixture(
  n = 90,
  Mu = Mu,
  Sigma = Sigma,
  family = "t",
  t_df = 5,
  seed = 1
)

fit <- gem_fit(
  dat$X,
  K = K,
  max_outer = 5,
  init_n_boot = 2,
  init_nstart = 4,
  outer_nstart = 1
)

table(fit$cluster)
cluster_accuracy(fit$cluster, dat$label)$accuracy
```

Select `K` with the GAP-LSE rule:

```r
gap <- gem_select_k_gap(
  dat$X,
  K_grid = 2:5,
  B = 20,
  control = list(
    max_outer = 5,
    init_n_boot = 2,
    init_nstart = 4,
    outer_nstart = 1
  )
)

gap$best_K
gap$summary
```

For large simulation studies, increase `max_outer`, `init_n_boot`,
`init_nstart`, `outer_nstart`, and `B` according to the desired accuracy and
computational budget.
