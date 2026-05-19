\name{cluster_accuracy}
\alias{cluster_accuracy}
\title{Clustering accuracy}
\usage{
cluster_accuracy(pred, truth)
}
\arguments{
\item{pred}{Predicted labels.}
\item{truth}{True labels.}
}
\description{
Computes clustering accuracy after the best label permutation.
}
\value{
A list with \code{accuracy} and the relabeled vector \code{mapped}.
}
\examples{
truth <- c(1, 1, 2, 2, 3, 3)
pred <- c(2, 2, 3, 3, 1, 1)
cluster_accuracy(pred, truth)$accuracy
}
