
## makeTSP

makeTSP <- function(data, method, path, dist.method, dist.p) {
  assert_character(path, len = 1, any.missing = FALSE)
  assert_choice(dist.method, choices = c("euclidean", "maximum", "manhattan",
    "canberra", "binary", "minkowski"))
  assert_numeric(dist.p, lower = 1, len = 1, any.missing = FALSE)
  if(method == "nn") {
    tsp.method <- "repetitive_nn"
  } else {
    invisible(capture.output(concorde_path(path)))
    tsp.method <- "concorde"
  }
  tsp <- TSP(dist(data, method = dist.method, p = dist.p))
  tsp <- insert_dummy(tsp)
  tour <- suppressWarnings(solve_TSP(tsp, method = tsp.method, verbose = FALSE))
  tour_tmp <- cut_tour(tour, max(tour))
  names(tour_tmp) <- NULL
  return(tour_tmp)
}
