
#' @title The Distribution of the Sign Depth
#' @name depthDistribution
#'
#' @description
#' Density, distribution function and quantile function for the distribution of the Sign Depth.
#'
#' @param x,q [\code{numeric}]\cr vector of quantiles.
#' @param p [\code{numeric}]\cr vector of probabilities.
#' @param n [\code{integer(1)}]\cr number of observations.
#' @param k [\code{integer(1)}]\cr parameter of the Sign Depth.
#' Currently only \eqn{k = 2, 3, 4, 5} is supported. Default is \eqn{k = 3}.
#' @param transform [\code{logical(1)}]\cr
#' Shall the values be transformed so that the distribution converges for large n?
#'
#' @return \code{ddepth} gives the density, \code{pdepth} the distribution function
#' and \code{qdepth} the quantile function.
#'
#' @details
#' For \eqn{n = k, \ldots, 25} the depth distributions are calculated exactly.
#' For \eqn{n = 26, \ldots, 100} the distribution is simulated.
#' For \eqn{n > 100} the values of \eqn{n = 100} are token
#' (the convergence of the distribution is sufficiently at that point).
#'
#' If y is a non-transformed value of one of these functions, \eqn{n * (y - 0.5^(k-1))} is the transformed value.
#' For further understanding of the transformed statistic, see references.
#'
#' @references
#' Kustosz C., Leucht A. and Mueller Ch. H. (2016). Tests based on Simplicial depth for AR(1) models with explosion.
#' Journal of Time Series Analysis. In press.
#'
#' @examples
#' qdepth(0.05, 20)
#' ddepth(0.25, 100)
#' pdepth(0.25, 50)
#'
#' @name depthDistribution
#' @rdname depthDistribution
NULL


#' @export
#' @rdname depthDistribution
qdepth <- function(p, n, k = 3, transform = FALSE) {
  assert_numeric(p, lower = 0, upper = 1, any.missing = FALSE, min.len = 1)
  assert_integerish(k, lower = 2, upper = 5, any.missing = FALSE, len = 1)
  assert_integerish(n, lower = k, len = 1, any.missing = FALSE)
  assert_logical(transform, any.missing = FALSE, len = 1)

  n <- min(n, 100)
  p <- round(p, 4)
  assign("dat", get(paste0("quants", k)))

  res <- dat[as.character(p), as.character(n)]
  if(transform) {
    res <- n * (res - (1/2)^(k-1))
  }
  return(res)
}


#' @export
#' @rdname depthDistribution
pdepth <- function(q, n, k = 3, transform = FALSE) {
  assert_numeric(q, any.missing = FALSE, min.len = 1)
  assert_integerish(k, lower = 2, upper = 5, any.missing = FALSE, len = 1)
  assert_integerish(n, lower = k, len = 1, any.missing = FALSE)
  assert_logical(transform, any.missing = FALSE, len = 1)

  n <- min(n, 100)
  q <- round(q, 4)
  assign("dat", get(paste0("quants", k)))

  dat <- dat[[as.character(n)]]
  if(transform) {
    dat <- n * (dat - (1/2)^(k-1))
  }
  pf <- ecdf(dat)
  return(pf(q))
}


#' @export
#' @rdname depthDistribution
ddepth <- function(x, n, k = 3, transform = FALSE) {
  assert_numeric(x, any.missing = FALSE, min.len = 1)
  assert_integerish(k, lower = 2, upper = 5, any.missing = FALSE, len = 1)
  assert_integerish(n, lower = k, len = 1, any.missing = FALSE)
  assert_logical(transform, any.missing = FALSE, len = 1)

  n <- min(n, 100)
  assign("dat", get(paste0("quants", k)))

  dat <- dat[[as.character(n)]]
  if(transform) {
    dat <- n * (dat - (1/2)^(k-1))
  }
  res <- density(dat, n = 2^13)
  inds <- sapply(x, function(y) which.min(abs(res$x - y)))
  return(ifelse(inds %in% c(1, length(res$x)), 0, res$y[inds]))
}