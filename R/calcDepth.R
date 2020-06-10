#' @title K Sign Depth
#'
#' @description
#' \code{calcDepth} calculates the (approximative) K Sign Depth of a given
#' vector of residuals.
#'
#' @param res [\code{numeric}]\cr
#'   numeric vector of residuals
#' @param K [\code{integer(1)}]\cr
#'   the K of the depth-formula
#' @param transform [\code{logical(1)}]\cr
#'   Should the depth be transformed by the formula \eqn{N * (res - 0.5^(K-1))},
#'   so that the distribution of the K Sign Depth converges for large number
#'   of observations N? Default is \code{FALSE}.
#' @param linear [\code{logical(1)}]\cr
#'   Should the depth be calculated in linear runtime? At the moment, for exact
#'   results this method is available only for K = 2, 3, 4, 5. Default is \code{TRUE}.
#' @param exact [\code{logical(1)}]\cr
#'   Should the depth be calculated exactly? If \code{FALSE}, terms of scale
#'   o(1) will ne ignored. If \code{linear} is \code{FALSE} always the exact
#'   result is computed, independent of the value of \code{exact}. Default is
#'   \code{TRUE}.
#' @return [\code{numeric(1)}] the calculated K Sign Depth of \code{res}.
#'
#' @note The implementation of the linear calculation is based on a work of
#'   Kevin Leckey, Dennis Malcherczyk and Christine Mueller.
#'
#' @examples
#' calcDepth(rnorm(10), 3)
#'
#' # Difference of exact and approximative implementation:
#' res <- rnorm(100)
#' calcDepth(res, 4, exact = TRUE)
#' calcDepth(res, 4, exact = FALSE)
#' @export
calcDepth <- function(res, K, transform = FALSE, linear = TRUE, exact = TRUE) {
  assert_numeric(res, min.len = K, any.missing = FALSE)
  assert_integerish(K, lower = 2, len = 1, any.missing = FALSE)
  assert_logical(transform, any.missing = FALSE, len = 1)
  assert_logical(linear, any.missing = FALSE, len = 1)
  assert_logical(exact, any.missing = FALSE, len = 1)

  res2 <- res[res != 0]
  N <- length(res)
  M <- length(res2)


  if (!linear) {
    if (K %in% 3:5)
      warning("It is recommended to use the linear implementation!")
    if (N == M) {
      erg <- calcDepth_def(res, K)
    } else {
      erg <- (choose(M, K) * calcDepth_def(res2, K) + choose(N - M, K)) / choose(N, K)
    }
    if (transform) erg <- length(res) * (erg - (1/2)^(K - 1))
  } else {
    erg <- asymp_K_depth(res2, K)
    if (exact) {
      if (K == 4) {
        S <- cumsum(sign(res2))
        erg <- erg + M * linearprod(res2, 4) / (choose(M, 4) * 8) +
          M^3 / ((M - 1) * (M - 2) * (M - 3)) * 3 * (1 / (4 * M^2) - 1 / (2 * M^3)) *
          (S[M]^2 - M)
      } else if (K == 5) {
        S <- cumsum(sign(res2))
        z1 <- sum(S^2)
        erg <- erg + M / (16 * choose(M, 5)) * (M * linearprod(res2, 4) -
            2 * prod_one_factor(res2, 4, 4) + 2 * prod_one_factor(res2, 4, 3) -
            2 * prod_one_factor(res2, 4, 2) + 2 * prod_one_factor(res2, 4, 1)) +
          M^4 / (32 * choose(M, 5)) * ((1 / (2 * M) * (S[M]^2 - M) -
              4 / (3 * M^2) * (S[M]^2 - M) - 1/M * S[M]^2 + 1/M^2 * (z1 + sum((S[M] - S)^2)) +
              8 / (3 * M^2) * S[M]^2 - 8 / (3 * M^3) * (z1 + sum((S[M] - S)^2))))
      } else if (K >= 6) {
        warning("No exact linear implementation available. Use linear = FALSE instead.
          The approximative result will be returned.")
      }
    }
    erg <- erg / length(res2) + (1/2)^(K - 1)
    if (M != N) {
      erg <- (choose(M, K) * erg + choose(N - M, K)) / choose(N, K)
    }
    if (transform) erg <- length(res) * (erg - (1/2)^(K - 1))
  }
  return(erg)
}
