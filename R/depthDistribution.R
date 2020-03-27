
#' @title The Distribution of the Sign Depth
#' @name depthDistribution
#'
#' @description
#' Density, distribution function and quantile function for the distribution of the K Sign Depth.
#'
#' @param x,q [\code{numeric}]\cr vector of quantiles.
#' @param p [\code{numeric}]\cr vector of probabilities.
#' @param N [\code{integer(1)}]\cr number of observations.
#' @param K [\code{integer(1)}]\cr parameter of the Sign Depth.
#' Currently only \eqn{K = 2, 3, 4, 5} is supported. Default is \eqn{K = 3}.
#' @param transform [\code{logical(1)}]\cr
#' Shall the values be transformed so that the distribution converges for large N?
#'
#' @return \code{ddepth} gives the density, \code{pdepth} the distribution function
#' and \code{qdepth} the quantile function.
#'
#' @details
#' For \eqn{N = K, \ldots, 25} the depth distributions are calculated exactly.
#' For \eqn{N = 26, \ldots, 100} the distribution is simulated.
#' For \eqn{N > 100} the values of \eqn{n = 100} are token
#' (the convergence of the distribution is sufficiently at that point).
#'
#' If y is a non-transformed value of one of these functions, \eqn{N * (y - 0.5^(K-1))} is the transformed value.
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
#' @rdname depthDistribution
NULL


#' @export
#' @rdname depthDistribution
qdepth <- function(p, N, K = 3, transform = FALSE) {
  assert_numeric(p, lower = 0, upper = 1, any.missing = FALSE, min.len = 1)
  assert_integerish(K, lower = 2, upper = 5, any.missing = FALSE, len = 1)
  assert_integerish(N, lower = K, len = 1, any.missing = FALSE)
  assert_logical(transform, any.missing = FALSE, len = 1)

  N <- min(N, 100)
  p <- round(p, 4)
  assign("dat", get(paste0("quants", K)))

  res <- dat[as.character(p), as.character(N)]
  if(transform) {
    res <- N * (res - (1/2)^(K-1))
  }
  return(res)
}


#' @export
#' @rdname depthDistribution
pdepth <- function(q, N, K = 3, transform = FALSE) {
  assert_numeric(q, any.missing = FALSE, min.len = 1)
  assert_integerish(K, lower = 2, upper = 5, any.missing = FALSE, len = 1)
  assert_integerish(N, lower = K, len = 1, any.missing = FALSE)
  assert_logical(transform, any.missing = FALSE, len = 1)

  N <- min(N, 100)
  q <- round(q, 4)
  assign("dat", get(paste0("quants", K)))

  dat <- dat[[as.character(N)]]
  if(transform) {
    dat <- N * (dat - (1/2)^(K-1))
  }
  pf <- ecdf(dat)
  return(pf(q))
}


#' @export
#' @rdname depthDistribution
ddepth <- function(x, N, K = 3, transform = FALSE) {
  assert_numeric(x, any.missing = FALSE, min.len = 1)
  assert_integerish(K, lower = 2, upper = 5, any.missing = FALSE, len = 1)
  assert_integerish(N, lower = K, len = 1, any.missing = FALSE)
  assert_logical(transform, any.missing = FALSE, len = 1)

  N <- min(N, 100)
  assign("dat", get(paste0("quants", K)))

  dat <- dat[[as.character(N)]]
  if(transform) {
    dat <- N * (dat - (1/2)^(K-1))
  }
  res <- density(dat, n = 2^13)
  inds <- sapply(x, function(y) which.min(abs(res$x - y)))
  return(ifelse(inds %in% c(1, length(res$x)), 0, res$y[inds]))
}
