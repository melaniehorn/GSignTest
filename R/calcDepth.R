#' @title Calculating the K-sign depth
#'
#' @description
#' \code{calcDepth} calculates the (approximative) K-sign depth of a given
#' vector of residuals.
#'
#' @param res [\code{numeric}]\cr
#'   numeric vector of residuals
#' @param K [\code{integer(1)}]\cr
#'   the parameter K for the sign depth
#' @param transform [\code{logical(1)}]\cr
#'   Should the depth be transformed by the formula \eqn{N * (res - 0.5^(K-1))},
#'   so that the distribution of the K-sign depth converges for large number
#'   of data points N? Default is \code{FALSE}.
#' @param linear [\code{logical(1)}]\cr
#'   Should the depth be calculated in linear time complexity? Default is \code{TRUE}.
#' @return [\code{numeric(1)}] the calculated K-sign depth of \code{res}.
#'
#' @note The implementation of the linear calculation is based on a work of
#'   Dennis Malcherczyk.
#'
#' @references
#' Horn M. Sign depth for parameter tests in multiple regression.
#' TU Dortmund University; 2021.
#'
#' Malcherczy D. K-sign depth: Asymptotic distribution, efficient computation
#' and applications. TU Dortmund University; 2021.
#'
#' @examples
#' calcDepth(rnorm(10), 3)
#'
#' # Difference of exact and approximative implementation:
#' res <- rnorm(100)
#' calcDepth(res, 4)
#' calcDepth(res, 4, transform = TRUE)
#' @export
calcDepth <- function(res, K, transform = FALSE, linear = TRUE) {
  assert_numeric(res, min.len = K, any.missing = FALSE)
  assert_integerish(K, lower = 2, len = 1, any.missing = FALSE)
  assert_logical(transform, any.missing = FALSE, len = 1)
  assert_logical(linear, any.missing = FALSE, len = 1)

  res2 <- res[res != 0]
  N <- length(res)
  M <- length(res2)


  if (!linear) {
    warning("It is recommended to use the linear implementation!")
    if (N == M) {
      erg <- calcDepth_def(res, K)
    } else {
      erg <- (choose(M, K) * calcDepth_def(res2, K) + choose(N - M, K)) / choose(N, K)
    }
    if (transform) erg <- length(res) * (erg - (1/2)^(K - 1))
  } else {
    signRes <- sign(res2)
    if (transform) {
      return((RcppCalcKDepthBlock(signRes, K = K) / choose(M, K) - (1/2)^(K-1)) * M)
    } else {
      return(RcppCalcKDepthBlock(signRes, K = K) / choose(M, K))
    }
    if (M != N) {
      erg <- (choose(M, K) * erg + choose(N - M, K)) / choose(N, K)
    }
  }
  return(erg)
}
