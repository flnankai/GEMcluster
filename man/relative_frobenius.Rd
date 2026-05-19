\name{relative_frobenius}
\alias{relative_frobenius}
\title{Relative Frobenius error}
\usage{
relative_frobenius(A, B)
}
\arguments{
\item{A}{Estimated matrix.}
\item{B}{Target matrix.}
}
\description{
Computes \eqn{\|A-B\|_F / \|B\|_F}.
}
\value{
A nonnegative numeric scalar.
}
\examples{
A <- diag(2)
B <- matrix(c(1, 0.1, 0.1, 1), 2, 2)
relative_frobenius(A, B)
}
