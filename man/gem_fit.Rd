\name{gem_fit}
\alias{gem_fit}
\title{Fit GEM clustering}
\usage{
gem_fit(X, K, max_outer = 25L, tol = 2e-3, init_n_boot = 5L,
  init_nstart = 6L, init_max_iter = 50L, outer_nstart = 3L,
  eta_mu = 0.7, eta_omega = 0.7, lambda_u = NULL,
  lambda_omega = NULL, lambda_u_scale = 0.55,
  lambda_omega_scale = 0.28, m_gr_upper = 8L,
  bandwidth = NULL, verbose = FALSE)
}
\arguments{
\item{X}{n by p data matrix.}
\item{K}{Number of clusters.}
\item{max_outer}{Maximum number of GEM outer iterations.}
\item{tol}{Convergence tolerance.}
\item{init_n_boot}{Number of permutation bootstraps for sparse K-median tau selection.}
\item{init_nstart}{Number of sparse K-median starts.}
\item{init_max_iter}{Maximum sparse K-median iterations.}
\item{outer_nstart}{Number of outer GEM starts; the largest final objective is retained.}
\item{eta_mu}{Mean update damping.}
\item{eta_omega}{Precision update damping.}
\item{lambda_u}{Optional POET threshold.}
\item{lambda_omega}{Optional graphical-lasso threshold base.}
\item{lambda_u_scale}{Default POET threshold scale.}
\item{lambda_omega_scale}{Default graphical-lasso threshold scale.}
\item{m_gr_upper}{Maximum factor count in the GR factor-number estimator.}
\item{bandwidth}{Optional kernel bandwidth for the nonparametric generator.}
\item{verbose}{If TRUE, print progress.}
}
\description{
Fits the semiparametric elliptical mixture GEM clustering method using sparse
K-median initialization and Rcpp kernels.
}
\details{
The algorithm estimates cluster-specific centers, mixing proportions, a common
sparse precision-shape matrix, and an unknown common radial generator. The
precision block combines weighted Tyler-type shape estimation, POET
regularization, and graphical lasso. The fitted object is selected across
\code{outer_nstart} starts by the final GEM objective.
}
\value{
An object of class \code{"gem_fit"}, represented as a list with entries
\code{cluster}, \code{Tau}, \code{pi}, \code{Mu}, \code{Sigma}, \code{Omega},
\code{generator}, \code{iterations}, \code{objective}, and initialization
diagnostics.
}
\examples{
if (requireNamespace("mvtnorm", quietly = TRUE) &&
    requireNamespace("huge", quietly = TRUE)) {
  set.seed(1)
  Mu <- dense_block_means(K = 3, p = 6, signal = 0.8)
  Sigma <- make_ar_cov(p = 6, rho = 0.5)
  dat <- simulate_elliptical_mixture(
    n = 36, Mu = Mu, Sigma = Sigma, family = "gaussian", seed = 1
  )

  fit <- gem_fit(
    dat$X, K = 3,
    max_outer = 1,
    init_n_boot = 1,
    init_nstart = 2,
    init_max_iter = 5,
    outer_nstart = 1
  )
  table(fit$cluster)
}
}
