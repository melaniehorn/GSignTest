
# linearprod: helper-function: Will be needed for calculation of the
#             exact depth in linear runtime
# input: residuals: Vector of residuals; K: Parameter of the K Sign Depth
# output: Needed result
linearprod <- function(residuals, K = 4){
    M <- sum(residuals < 0)
    N <- length(residuals)
    C <- sum(choose(M, 0:K) * choose(N - M, K - 0:K) * (-1)^(0:K))
    return(C)
  }
