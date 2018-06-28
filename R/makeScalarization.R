
## makeScalarization

makeScalarization <- function(data, weights) {
  assert_numeric(weights, lower = 0, len = ncol(data), any.missing = FALSE)
  tmp <- rowSums(sapply(1:ncol(data), function(x) data[, x] * weights[x]))
  inds <- order(tmp)
  return(inds)
}
