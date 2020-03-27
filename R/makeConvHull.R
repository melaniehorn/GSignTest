
### makeConvHull

makeConvHull <- function(data, path) {
  tryMethod <- try(capture.output(concorde_path(path)), silent = TRUE)
  method <- if (class(tryMethod) == "try-error") "nn" else "shp"

  if (ncol(data) == 1) {
    return(multiSorting(data, method = method, control = list(path = path))$inds)
  }

  method <- ifelse(method == "nn", "repetitive_nn", "concorde")

  tmp <- data
  erg <- list()
  i <- 1
  repeat {
    if(nrow(tmp) <= ncol(tmp)) {
      if(nrow(tmp) > 0) erg[[i]] <- tmp
      break
    }
    ind <- try(as.numeric(geometry::convhulln(tmp)), silent = TRUE)
    if(class(ind) == "try-error") {
      erg[[i]] <- tmp
      tmp <- tmp[-(1:nrow(tmp)), , drop = FALSE]
    } else {
      erg[[i]] <- tmp[unique(ind), ]
      tmp <- tmp[-unique(ind), , drop = FALSE]
    }
    i <- i + 1
  }

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
