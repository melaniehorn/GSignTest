
makeProjection <- function(data, vec.pos, vec.dir) {
  assert_numeric(vec.pos, any.missing = FALSE, len = ncol(data))
  assert_numeric(vec.dir, any.missing = FALSE, len = ncol(data))

  denom <- sum(vec.dir^2)
  erg <- apply(data, 1, function(x) {
    vec.pos + sum((x - vec.pos) * vec.dir) / denom * vec.dir
  })
  erg <- if(ncol(data) == 1) erg else t(erg)
  erg <- as.data.frame(erg)

  return(makeDimSort(erg, 1:ncol(data)))
}
