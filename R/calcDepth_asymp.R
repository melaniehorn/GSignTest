#' @title Calculating the K-sign depth (asymptotically)
#'
#' @description
#' \code{calcDepth_asymp} calculates the K-sign depth of a given
#' vector of residuals (asymptotically).
#'
#' @param res [\code{numeric}]\cr
#'   numeric vector of residuals
#' @param K [\code{integer(1)}]\cr
#'   the parameter K for the sign depth
#' @param transform [\code{logical(1)}]\cr
#'   Should the depth be transformed by the formula \eqn{N * (res - 0.5^(K-1))},
#'   so that the distribution of the K-sign depth converges for large number
#'   of data points N? Default is \code{FALSE}.
#' @param exact [\code{logical(1)}]\cr
#'   Should the depth be calculated exactly? If \code{FALSE}, terms of scale
#'   o(1) will ne ignored. At the moment, for exact
#'   results this method is available only for K = 2, 3, 4, 5. Default is
#'   \code{TRUE}.
#' @return [\code{numeric(1)}] the calculated K-sign depth of \code{res}.
#'
#' @note The implementation is based on a work of Dennis Malcherczyk.
#'
#' @references
#' Malcherczyk, D. (2021+). K-sign depth: Asymptotic distribution, efficient
#' computation and applications. Dissertation in preparation.
#'
#'
#' @examples
#' calcDepth_asymp(rnorm(10), 3)
#' calcDepth_asymp(runif(100, -1, 1), 4, transform = TRUE)
#' calcDepth_asymp(runif(100, -1, 1), 4, transform = TRUE, exact = FALSE)
#' @export

calcDepth_asymp <- function(res, K, transform = FALSE, exact = TRUE) {
  assert_numeric(res, min.len = K, any.missing = FALSE)
  assert_integerish(K, lower = 2, len = 1, any.missing = FALSE)
  assert_logical(transform, any.missing = FALSE, len = 1)
  assert_logical(exact, any.missing = FALSE, len = 1)

  res2 <- res[res != 0]
  N <- length(res)
  M <- length(res2)

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
    erg <- erg / length(res2) + (1/2)^(K - 1)
    if (M != N) {
      erg <- (choose(M, K) * erg + choose(N - M, K)) / choose(N, K)
    }
    if (transform) erg <- length(res) * (erg - (1/2)^(K - 1))
  }
  return(erg)
}
