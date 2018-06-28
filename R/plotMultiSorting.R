
#' @title Plot the sorted data (two-dimensional)
#'
#' @description
#' Plots two dimensions of the data as a scatter plot and the sorting of the
#' data as lines, which connect the points. Therefore the function ggplot from
#' the ggplot2-package is used.
#'
#' @param x [\code{multiSorting}]\cr An object returned from \code{\link{multiSorting}}.
#' @param xlab [\code{character(1)}]\cr Optional x axis label.
#' @param ylab [\code{character(1)}]\cr Optional y axis label.
#' @param main [\code{character(1)}]\cr Optional title of the plot.
#' @param dims [\code{integer(2)}]\cr The dimensions to be plottet. Default is the first two.
#' @param point.size [\code{numeric(1)}]\cr Non-negative size of the points in the plot.
#' @param path.size [\code{numeric(1)}]\cr Non-negative size of the path in the plot.
#'
#' @return The ggplot object
#'
#' @seealso \code{\link{multiSorting}}
#' @examples
#' x <- data.frame(x = sort(rnorm(20)), y = sort(rnorm(20)))
#' multSort <- multiSorting(x)
#' plotMultiSorting(multSort, point.size = 3, path.size = 2)
#'
#' @export
plotMultiSorting <- function(x, xlab, ylab, main, dims = 1:2, point.size = 1, path.size = 1) {
  assert_class(x, "multiSorting")
  assert_integerish(dims, lower = 1, upper = ncol(x$sortedData),
    any.missing = FALSE, len = 2)
  assert_numeric(point.size, lower = 1, any.missing = FALSE, len = 1)
  assert_numeric(path.size, lower = 1, any.missing = FALSE, len = 1)
  if(ncol(x$sortedData) < 2) {
    stop("Plots for one-dimensional data not supported!")
  }

  if(missing(xlab)) {
    xlab <- "x"
  }
  if(missing(ylab)) {
    ylab <- "y"
  }
  if(missing(main)) {
    main <- ""
  }

  dat <- x$sortedData[, dims]
  names(dat) <- c("x", "y")

  p <- ggplot(mapping = aes(x = dat$x, y = dat$y)) + geom_path(size = path.size) +
    geom_point(color = "red", size = point.size) + xlab(xlab) + ylab(ylab) +
    ggtitle(main)

  return(p)
}
