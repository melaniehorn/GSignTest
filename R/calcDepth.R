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

  if (!linear) {
    if (K %in% 3:5)
      warning("It is recommended to use the linear implementation!")
    erg <- calcDepth_def(res, K)
    if (transform) erg <- length(res) * (erg - (1/2)^(K - 1))
  } else {
    erg <- asymp_K_depth(res, K)
    if (exact) {
      if (K == 4) {
        N <- length(res)
        S <- cumsum(sign(res))
        erg <- erg + N * linearprod(res, 4) / (choose(N, 4) * 8) +
          N^3 / ((N - 1) * (N - 2) * (N - 3)) * 3 * (1 / (4 * N^2) - 1 / (2 * N^3)) *
          (S[N]^2 - N)
      } else if (K == 5) {
        N <- length(res)
        S <- cumsum(sign(res))
        z1 <- sum(S^2)
        erg <- erg + N / (16 * choose(N, 5)) * (N * linearprod(res, 4) -
            2 * prod_one_factor(res, 4, 4) + 2 * prod_one_factor(res, 4, 3) -
            2 * prod_one_factor(res, 4, 2) + 2 * prod_one_factor(res, 4, 1)) +
          N^4 / (32 * choose(N, 5)) * ((1 / (2 * N) * (S[N]^2 - N) -
              4 / (3 * N^2) * (S[N]^2 - N) - 1/N * S[N]^2 + 1/N^2 * (z1 + sum((S[N] - S)^2)) +
              8 / (3 * N^2) * S[N]^2 - 8 / (3 * N^3) * (z1 + sum((S[N] - S)^2))))
      } else if (K >= 6) {
        warning("No exact linear implementation available. Use linear = FALSE instead.
          The approximative result will be returned.")
      }
    }
    if (!transform) erg <- erg / length(res) + (1/2)^(K - 1)
  }
  return(erg)
}
