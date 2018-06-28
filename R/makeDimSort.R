
## makeDimSort

makeDimSort <- function(data, order) {
  assert_integerish(order, any.missing = FALSE, len = ncol(data), lower = 1,
    upper = ncol(data), unique = TRUE)

  return(do.call("order", data[order]))
}
