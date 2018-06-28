
## makeMedian

makeMedian <- function(data) {
  meds <- apply(data, 1, median)
  inds <- order(meds)
  return(inds)
}
