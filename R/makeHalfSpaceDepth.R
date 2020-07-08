
### makeHalfSpaceDepth

makeHalfSpaceDepth <- function(data, path) {
  tryMethod <- try(capture.output(concorde_path(path)), silent = TRUE)
  method <- if (class(tryMethod) == "try-error") "repetitive_nn" else "concorde"

  if(ncol(data) == 1) {
    data2 <- cbind(data, 0)
    depths <- ddalpha::depth.halfspace(data2, data2, num.directions = 100000)
  } else {
    depths <- ddalpha::depth.halfspace(data, data, num.directions = 100000)
  }

  erg <- split(data, depths)

  calcDist <- function(mat1, mat2) {
    dists <- apply(mat1, 1, function(x) {
      apply(mat2, 1, function(y) sqrt(sum((x - y)^2)))
    })
    return(which(dists == min(dists), arr.ind = TRUE))
  }

  inds <- c()
  for(i in seq_along(erg)) {
    if(nrow(erg[[i]]) > 2) {
      tsp <- TSP(dist(erg[[i]]))
      tour <- suppressWarnings(solve_TSP(tsp, method = method, verbose = FALSE))
      names(tour) <- NULL
      inds <- c(inds, calcDist(erg[[i]], data)[tour, ][, 1])
    } else {
      inds <- c(inds, calcDist(erg[[i]], data)[1:nrow(erg[[i]]), , drop = FALSE][, 1])
      names(inds) <- NULL
    }
  }
  return(inds)
}
