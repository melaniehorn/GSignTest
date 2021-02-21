#' @title Multidimensional Sorting/Ordering
#'
#' @description Performs a multidimensional sorting of a \code{data.frame}.
#' One can choose between various approaches for the sorting, see the Details.
#'
#' @param data [\code{data.frame}]\cr \code{data.frame} where each row consists
#' of one numeric data point vector.
#' @param method [\code{character(1)}]\cr Method for sorting, see Details.
#' @param control [\code{list}]\cr Further Arguments for the methods, see Details.
#'
#' @return Named list with the entries:
#' \describe{
#'   \item{\code{sortedData} [\code{data.frame}]}{The sorted data.}
#'   \item{\code{inds} [\code{integer}]}{The sorted indices of the rows of the input data.}
#' }
#'
#' @details The following methods for sorting are supported:
#' \describe{
#'   \item{\code{"identity"}}{order in the data set (i.e. the function does no further sorting).}
#'   \item{\code{"random"}}{a random order.}
#'   \item{\code{"norm"}}{ordering via a vector norm.}
#'   \item{\code{"median"}}{ordering via the median value per data point vector.}
#'   \item{\code{"dimSort"}}{sorting the data one-dimensional. Ties are broken by the next dimension and so on.}
#'   \item{\code{"weightedSum"}}{ordering the data via the value of a weighted sum per data point vector.}
#'   \item{\code{"projection"}}{projects the points of the data set orthogonal on a line.}
#'   \item{\code{"nonDomSort"}}{Nondominated Sorting. Ties can be broken via one of the other methods.}
#'   \item{\code{"convhull"}}{Computing repeatedly convex hulls for ordering the data.}
#'   \item{\code{"halfspace"}}{ordering according to the (approximate) values of Tukey's halfspace depth.}
#'   \item{\code{"clust"}}{Take the order gotten from a hierarchical clustering by \code{\link[stats]{hclust}}}
#'   \item{\code{"nn"}}{Nearest-Neighbor-Heuristic.
#'       Tries to find a shortest path trough the data by always taking the next nearest not-so-far-chosen point.}
#'   \item{\code{"shp"}}{Computes the Shortest Hamiltonian Path through the data.}
#' }
#'
#' The \code{control} argument is a list that can supply any of the following components:
#' \describe{
#'   \item{\code{degree} [\code{numeric(1)}]}{For method \code{"norm"}:\cr
#'     Non-negative (possibly infinity) parameter of a vector-norm. Default: \code{2}.}
#'   \item{\code{order} [\code{integer}]}{For method \code{"dimSort"}:\cr
#'     Vector of dimension indices which determine the order of sorting. Default: \code{1:ncol(data)}.}
#'   \item{\code{weights} [\code{numeric}]}{For method \code{"weightedSum"}:\cr
#'     Non-negative vector of length \code{ncol(data)} which give the weight for each dimension.
#'     Default: \code{rep(1, ncol(data))}.}
#'   \item{\code{vec.pos} [\code{numeric}]}{For method \code{"projection"}:\cr
#'     Vector of length \code{ncol(data)} which gives the position vector of the line to project on.
#'     Default: \code{rep(0, ncol(data))}.}
#'   \item{\code{vec.dir} [\code{numeric}]}{For method \code{"projection"}:\cr
#'     Vector of length \code{ncol(data)} which gives the direction vector of the line to project on.
#'     Default: \code{rep(1, ncol(data))}.}
#'   \item{\code{tie.breaking} [\code{character(1)}]}{For method  \code{"nonDomSort"}:\cr
#'     A character string determining one of the other methods for tie-breaking of the domination ranks.
#'     Default: \code{"identity"}. Note, that if the chosen method has control arguments itsself,
#'     they can/should be also given via \code{control}.}
#'   \item{\code{clust.method} [\code{character(1)}]}{For method \code{"clust"}:\cr
#'     The \code{method}-argument of \code{\link[stats]{hclust}}. Default: \code{"complete"}.}
#'   \item{\code{dist.method} [\code{character(1)}]}{For method \code{"clust"}, \code{"nn"} and \code{"shp"}:\cr
#'     The \code{method}-argument of \code{\link[stats]{dist}}. Default: \code{"euclidean"}.}
#'   \item{\code{p} [\code{numeric(1)}]}{For method \code{"clust"}, \code{"nn"} and \code{"shp"}:\cr
#'     The \code{p}-argument of \code{\link[stats]{dist}}. Default: \code{2}.}
#'   \item{\code{path} [\code{character(1)}]}{For method \code{"convhull"}, \code{"halfspace"} and \code{"shp"}:\cr
#'     Path to the executable code of the Concorde solver. Default: The working directory.}
#'   \item{\code{prec} [\code{numeric(1)}]}{For method \code{"shp"}: \cr
#'     The precision, i.e. the number of decimal places using for the distances between two points.
#'     Default: \code{6}}.
#' }
#'
#' @note
#' The code of the Concorde TSP package is not included in this package and has to be obtained separately
#' from the Concorde web site (see references). Either download the precompiled executables
#' and place them in a suitable directory and make them executable, or you can get the source code and
#' compile it on your own. TSP needs to know where the executables are.
#'
#' @references
#' \itemize{
#'   \item Horn M. Sign depth for parameter tests in multiple regression. TU Dortmund University; 2021.
#'   \item Barber, C. B., Dobkin, D. P., and Huhdanpaa, H. (1996). “The Quickhull Algorithm for Convex Hulls”. In: ACM Trans. Math. Softw. 22.4, pp. 469–483.
#'   \item Tukey, J. W. (1975). “Mathematics and the picturing of data”. In: Proceedings of the International Congress of Mathematicians (Vancouver, BC, 1974), Volume 2. Montreal, Quebec, Canada: Canadian Mathematical Congress, pp. 523–531.
#'   \item Deb, K., Agrawal, S., Pratap, A., and Meyarivan, T. (2000). “A Fast Elitist Non- dominated Sorting Genetic Algorithm for Multi-objective Optimization: NSGA-II”. In: Parallel Problem Solving from Nature PPSN VI. Berlin, Heidelberg: Springer Berlin Heidelberg, pp. 849–858.
#'   \item Rockafellar, R. T. (1970). Convex analysis. Vol. 28. Princeton mathematical series. Princton University Press.
#'   \item Concorde home page, \url{http://www.tsp.gatech.edu/concorde/}
#'   \item Concorde download page, \url{http://www.tsp.gatech.edu/concorde/downloads/downloads.htm}
#'   \item David Appletgate, Robert Bixby, Vasek Chvatal, William Cook (2001): TSP cuts which do not
#'         conform to the template paradigm, Computational Combinatorial Optimization, M. Junger and D.
#'         Naddef (editors), Springer.
#' }
#'
#' @seealso \code{\link{plotMultiSorting}}
#'
#' @examples
#' multiSorting(iris[, 1:4], method = "identity")
#' multiSorting(iris[, 1:4], method = "random")
#' multiSorting(iris[, 1:4], method = "norm")
#' multiSorting(iris[, 1:4], method = "median")
#' multiSorting(iris[, 1:4], method = "dimSort")
#' multiSorting(iris[, 1:4], method = "weightedSum")
#' multiSorting(iris[, 1:4], method = "projection")
#' multiSorting(iris[, 1:4], method = "nonDomSort")
#' multiSorting(iris[, 1:4], method = "convhull")
#' multiSorting(iris[, 1:4], method = "halfspace")
#' multiSorting(iris[, 1:4], method = "clust")
#' multiSorting(iris[, 1:4], method = "nn")
#' \dontrun{multiSorting(iris[, 1:4], method = "shp")}
#'
#' @export
multiSorting <- function(data, method = "identity", control = list()) {
  assert_data_frame(data, types = "numeric", any.missing = FALSE, min.rows = 2,
    min.cols = 1)
  assert_choice(method, choices = c("identity", "random", "norm", "median",
    "dimSort", "weightedSum", "projection", "nonDomSort", "convhull", "halfspace",
    "clust", "nn", "shp"))
  assert_list(control, any.missing = FALSE)

  if(method == "identity") {
    inds <- 1:nrow(data)
  } else if(method == "random") {
    inds <- sample(nrow(data))
  } else if(method == "norm") {
    degree <- if(is.null(control[["degree"]])) 2 else control[["degree"]]
    inds <- makeNorm(data, degree)
  } else if(method == "median") {
    inds <- makeMedian(data)
  } else if(method == "dimSort") {
    order <- if(is.null(control[["order"]])) 1:ncol(data) else control[["order"]]
    inds <- makeDimSort(data, order)
  } else if(method == "weightedSum") {
    weights <- if(is.null(control[["weights"]])) rep(1, ncol(data)) else control[["weights"]]
    inds <- makeWeightedSum(data, weights)
  } else if(method == "projection") {
    vec.pos <- if(is.null(control[["vec.pos"]])) rep(0, ncol(data)) else control[["vec.pos"]]
    vec.dir <- if(is.null(control[["vec.dir"]])) rep(1, ncol(data)) else control[["vec.dir"]]
    inds <- makeProjection(data, vec.pos, vec.dir)
  } else if(method == "nonDomSort") {
    tie.breaking <- if(is.null(control[["tie.breaking"]])) "identity" else control[["tie.breaking"]]
    inds <- makeNonDomSort(data, tie.breaking, control)
  } else if(method == "convhull") {
    path <- if(is.null(control[["path"]])) getwd() else control[["path"]]
    inds <- makeConvHull(data, path)
  } else if(method == "halfspace") {
    path <- if(is.null(control[["path"]])) getwd() else control[["path"]]
    inds <- makeHalfSpaceDepth(data, path)
  } else if(method == "clust") {
    dist.method <- if(is.null(control[["dist.method"]])) "euclidean" else control[["dist.method"]]
    p <- if(is.null(control[["p"]])) 2 else control[["p"]]
    clust.method <- if(is.null(control[["clust.method"]])) "complete" else control[["clust.method"]]
    inds <- makeClust(data, clust.method, dist.method, p)
  } else {
    path <- if(is.null(control[["path"]])) getwd() else control[["path"]]
    dist.method <- if(is.null(control[["dist.method"]])) "euclidean" else control[["dist.method"]]
    p <- if(is.null(control[["p"]])) 2 else control[["p"]]
    prec <- if(is.null(control[["prec"]])) 6 else control[["prec"]]
    inds <- makeTSP(data, method, path, dist.method, p, prec)
  }

  res <- list(sortedData = data[inds, ], inds = inds)
  class(res) <- c(method, "multiSorting")
  return(res)
}
