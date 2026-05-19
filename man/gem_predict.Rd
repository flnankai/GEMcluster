\name{gem_predict}
\alias{gem_predict}
\alias{predict.gem_fit}
\title{Predict from a GEM fit}
\usage{
gem_predict(object, newdata = NULL, type = c("class", "posterior"), ...)
\method{predict}{gem_fit}(object, newdata = NULL, type = c("class", "posterior"), ...)
}
\arguments{
\item{object}{A fitted object returned by \code{gem_fit}.}
\item{newdata}{Optional new data matrix. If NULL, fitted observations are used.}
\item{type}{Either \code{"class"} or \code{"posterior"}.}
\item{...}{Ignored.}
}
\description{
Returns predicted cluster labels or posterior membership probabilities.
}
\value{
If \code{type = "class"}, an integer vector of labels. If
\code{type = "posterior"}, an n by K matrix of posterior probabilities.
}
\examples{
if (requireNamespace("mvtnorm", quietly = TRUE) &&
    requireNamespace("huge", quietly = TRUE)) {
  set.seed(3)
  Mu <- dense_block_means(K = 3, p = 5, signal = 0.8)
  Sigma <- make_ar_cov(p = 5, rho = 0.4)
  dat <- simulate_elliptical_mixture(
    n = 30, Mu = Mu, Sigma = Sigma, family = "gaussian", seed = 3
  )
  fit <- gem_fit(
    dat$X, K = 3,
    max_outer = 1,
    init_n_boot = 1,
    init_nstart = 2,
    init_max_iter = 5,
    outer_nstart = 1
  )
  head(gem_predict(fit))
  dim(gem_predict(fit, type = "posterior"))
}
}
