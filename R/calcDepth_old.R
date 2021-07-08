#' @title Calculating the K-sign depth (deprecated)
#'
#' @description
#' \code{calcDepth_old} calculates the K-sign depth of a given
#' vector of residuals in time complexity O(N^K)
#'
#' @param res [\code{numeric}]\cr
#'   numeric vector of residuals
#' @param K [\code{integer(1)}]\cr
#'   the parameter K for the sign depth
#' @param transform [\code{logical(1)}]\cr
#'   Should the depth be transformed by the formula \eqn{N * (res - 0.5^(K-1))},
#'   so that the distribution of the K-sign depth converges for large number
#'   of data points N? Default is \code{FALSE}.
#' @return [\code{numeric(1)}] the calculated K-sign depth of \code{res}.
#'
#'
#' @references
#' Horn M. Sign depth for parameter tests in multiple regression.
#' TU Dortmund University; 2021.
#'
#' @examples
#' calcDepth_old(rnorm(10), 3)
#' calcDepth_old(runif(100, -1, 1), 4, transform = TRUE)
#'
#' @export
calcDepth_old <- function(res, K, transform = FALSE) {
  assert_numeric(res, min.len = K, any.missing = FALSE)
  assert_integerish(K, lower = 2, len = 1, any.missing = FALSE)
  assert_logical(transform, any.missing = FALSE, len = 1)

  res2 <- res[res != 0]
  N <- length(res)
  M <- length(res2)


  warning("It is recommended to use the linear implementation in calcDepth()!")
  if (N == M) {
    erg <- calcDepth_def(res, K)
  } else {
    erg <- (choose(M, K) * calcDepth_def(res2, K) + choose(N - M, K)) / choose(N, K)
  }
  if (transform) erg <- length(res) * (erg - (1/2)^(K - 1))

  return(erg)
}
