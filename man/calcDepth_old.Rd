% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcDepth_old.R
\name{calcDepth_old}
\alias{calcDepth_old}
\title{Calculating the K-sign depth (deprecated)}
\usage{
calcDepth_old(res, K, transform = FALSE)
}
\arguments{
\item{res}{[\code{numeric}]\cr
numeric vector of residuals}

\item{K}{[\code{integer(1)}]\cr
the parameter K for the sign depth}

\item{transform}{[\code{logical(1)}]\cr
Should the depth be transformed by the formula \eqn{N * (res - 0.5^(K-1))},
so that the distribution of the K-sign depth converges for large number
of data points N? Default is \code{FALSE}.}
}
\value{
[\code{numeric(1)}] the calculated K-sign depth of \code{res}.
}
\description{
\code{calcDepth_old} calculates the K-sign depth of a given
vector of residuals in time complexity O(N^K)
}
\examples{
calcDepth_old(rnorm(10), 3)
calcDepth_old(runif(100, -1, 1), 4, transform = TRUE)

}
\references{
Horn M. Sign depth for parameter tests in multiple regression.
TU Dortmund University; 2021.
}
