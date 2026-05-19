\name{simulate_elliptical_mixture}
\alias{simulate_elliptical_mixture}
\alias{make_ar_cov}
\alias{dense_block_means}
\title{Simulation helpers}
\usage{
simulate_elliptical_mixture(n, Mu, Sigma, family = .elliptical_family_levels,
  t_df = 5, slash_df = 4, pi = NULL, seed = NULL)
make_ar_cov(p, rho = 0.5)
dense_block_means(K = 3, p = 200, signal = 0.5)
}
\arguments{
\item{n}{Total sample size.}
\item{Mu}{K by p mean matrix.}
\item{Sigma}{Common p by p scatter/covariance matrix.}
\item{family}{One of gaussian, t, laplace, slash.}
\item{t_df}{Degrees of freedom for t.}
\item{slash_df}{Tail parameter for slash; must be greater than 2.}
\item{pi}{Optional class proportions.}
\item{seed}{Optional random seed.}
\item{p}{Dimension.}
\item{rho}{AR correlation.}
\item{K}{Number of classes.}
\item{signal}{Signal size multiplying the block mean pattern.}
}
\description{
Utility functions for generating the elliptical mixture simulation settings
used with GEMcluster.
}
\details{
\code{simulate_elliptical_mixture()} uses the total sample size \code{n}; when
\code{pi = NULL}, class sizes are as balanced as possible. The t and slash
families are scaled so that the target scatter is \code{Sigma} when the
degrees-of-freedom parameter is greater than 2.
}
\value{
\code{simulate_elliptical_mixture()} returns a list with the data matrix
\code{X}, class labels \code{label}, and family name. \code{make_ar_cov()}
returns a covariance matrix. \code{dense_block_means()} returns a K by p mean
matrix for the dense block-constant design.
}
\examples{
if (requireNamespace("mvtnorm", quietly = TRUE)) {
  Mu <- dense_block_means(K = 3, p = 8, signal = 0.5)
  Sigma <- make_ar_cov(p = 8, rho = 0.5)
  dat <- simulate_elliptical_mixture(
    n = 30, Mu = Mu, Sigma = Sigma, family = "t", t_df = 5, seed = 1
  )
  dim(dat$X)
  table(dat$label)
}
}
