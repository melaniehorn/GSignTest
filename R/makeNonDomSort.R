
## makeNonDomSort

makeNonDomSort <- function(data, tie.breaking = "identity", control) {
  assert_choice(tie.breaking, choices = c("identity", "random", "norm",
    "dimSort", "projection", "nn", "shp"))

  res <- fastNonDominatedSorting(data)
  if(tie.breaking == "random") {
    res <- lapply(res, function(x) if(length(x) == 1) x else sample(x))
  } else if(tie.breaking != "identity") {
    res <- lapply(res, function(x) {
      if(length(x) == 1){
        x
      } else {
        x[multiSorting(data[x, ], method = tie.breaking, control)$inds]
      }
    })
  }
  return(unlist(res))
}
