
## makeNorm

makeNorm <- function(data, degree) {
  assert_numeric(degree, len = 1, lower = 1, upper = Inf, any.missing = FALSE)
  if(is.infinite(degree)) {
    res <- apply(data, 1, max)
  } else {
    res <- apply(data, 1, function(x) sum(abs(x)^degree)^(1/degree))
  }
  return(order(res))
}
