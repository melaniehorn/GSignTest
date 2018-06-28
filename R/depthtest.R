
#' @title Sign Depth Test
#'
#' @description
#' Performs a robust test whether a statistical model fits to some data.
#'
#' @param x [\code{numeric}]\cr vector of residuals.
#' @param k [\code{integer(1)}]\cr parameter of the Sign Depth.
#' @param formula [\code{formula}]\cr a formula of the form \code{lhs ~ rhs} where
#' lhs is a numeric giving the target variable and rhs are the influencing variables .
#' @param data [\code{data.frame}]\cr a data frame containing the variables in the formula.
#' @param params [\code{numeric}]\cr a parameter vector for the null-hypotheses.
#' @param ... [any] further arguments in the formula-method. Arguments are passed
#' to the \code{multiSorting}-function.
#'
#' @details
#' In the default method a residual vector is given and it is tested whether this
#' residual vector comes from a model that describes the data sufficiently.
#' With the formula interface it is possible to give a formula and a parameter vector,
#' so that a linear model with the given parameters is fitted and the residuals are
#' calculated with this model.
#'
#' At the moment only linear models for the formula interface are supported!
#'
#' The parameter k of the Sign Depth has to be a natural number. Currently
#' only \eqn{k = 2, 3, 4, 5} is supported!
#'
#' The quantiles used for calculating the p-value of the test are simulation based.
#' For more information see \code{\link{qdepth}}.
#'
#' @return A list with class \code{"htest"} containing the following components:
#'\describe{
#'  \item{\code{statistic}}{the value of the depth.}
#'  \item{\code{parameter}}{the k of the Sign Depth.}
#'  \item{\code{p.value}}{the p-value for the test.}
#'  \item{\code{alternative}}{a character string describing the alternative hypothesis.}
#'  \item{\code{data.name}}{a character string giving the name of the data.}
#'  \item{\code{method}}{a character string describing the method.}
#'}
#'
#' @seealso \code{\link{qdepth}}, \code{\link{calcDepth}}, \code{\link{multiSorting}},
#' \code{\link{sign.test}}, \code{\link{f.test}}
#'
#' @examples
#' depth.test(rnorm(30), 3)
#' depth.test(y ~ ., data = data.frame(x = rnorm(30), y = rnorm(30)),
#'   params = c(1, 1), k = 3)
#'
#' @rdname depthtest
#' @export
depth.test <- function(x, ...) {
  UseMethod("depth.test")
}


#' @rdname depthtest
#' @method depth.test default
#' @export
depth.test.default <- function(x, k, ...) {
  assert_numeric(x, min.len = 10, any.missing = FALSE)
  assert_integerish(k, lower = 2, upper = 5, len = 1, any.missing = FALSE)

  dname <- deparse(substitute(x))
  method <- "Sign Depth Test (via given residuals)"
  alternative <- "Model where the residuals come from fits not to the data"
  parameter <- k

  stat <- calcDepth(x, k)

  assign("dat", get(paste0("quants", k)))
  n <- min(length(x), 100)
  tmp <- which(dat[, as.character(n)] > stat)
  if(length(tmp) == 0) {
    p.value <- 1
  } else {
    ind <- min(tmp)
    p.value <- as.numeric(rownames(dat[ind, ]))
  }

  names(parameter) <- "k"
  names(stat) <- "depth"

  rval <- list(statistic = stat, parameter = parameter, p.value = p.value,
   alternative = alternative, data.name = dname, method = method)
  class(rval) <- "htest"
  return(rval)
}


#' @rdname depthtest
#' @method depth.test formula
#' @export
depth.test.formula <- function(formula, data, params, k, ...) {
  assert_data_frame(data, any.missing = FALSE, min.cols = 1)
  assert_numeric(params, min.len = 1, any.missing = FALSE)
  assert_integerish(k, lower = 2, upper = 5, len = 1, any.missing = FALSE)

  mm <- model.matrix(formula, data)
  mr <- model.response(model.frame(formula, data))

  if(ncol(mm) != length(params)) {
    stop(paste("Model matrix has", ncol(mm), "columns and params has length",
      length(params)))
  }

  sorted <- multiSorting(as.data.frame(mm), ...)$inds
  mm <- mm[sorted, ]
  mr <- mr[sorted]

  res <- mm %*% params - mr

  dname <- deparse(substitute(data))
  method <- "Sign Depth Test (via given parameters and data)"
  alternative <- paste("True parameter vector is not equal to", paste(params, collapse = " "))
  parameter <- k

  stat <- calcDepth(res, k)

  assign("dat", get(paste0("quants", k)))
  n <- min(length(res), 100)
  tmp <- which(dat[, as.character(n)] > stat)
  if(length(tmp) == 0) {
    p.value <- 1
  } else {
    ind <- min(tmp)
    p.value <- as.numeric(rownames(dat[ind, ]))
  }

  names(parameter) <- "k"
  names(stat) <- "depth"

  rval <- list(statistic = stat, parameter = parameter, p.value = p.value,
    alternative = alternative, data.name = dname, method = method)
  class(rval) <- "htest"
  return(rval)
}