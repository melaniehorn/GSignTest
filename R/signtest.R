
#' @title Sign Test
#'
#' @description
#' Performs a test whether a statistical model fits to some data via a binomial
#' test on the signs of residuals.
#'
#' @param x [\code{numeric}]\cr vector of residuals.
#' @param formula [\code{formula}]\cr a formula of the form \code{lhs ~ rhs} where
#' lhs is a numeric giving the target variable and rhs are the influencing variables .
#' @param data [\code{data.frame}]\cr a data frame containing the variables in the formula.
#' @param params [\code{numeric}]\cr a parameter vector for the null-hypotheses.
#' @param ... [\code{any}]\cr further arguments.
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
#' @return A list with class \code{"htest"} containing the following components:
#'\describe{
#'  \item{\code{statistic}}{the number of positive residuals.}
#'  \item{\code{parameter}}{the number of observations.}
#'  \item{\code{p.value}}{the p-value for the test.}
#'  \item{\code{alternative}}{a character string describing the alternative hypothesis.}
#'  \item{\code{data.name}}{a character string giving the name of the data.}
#'  \item{\code{method}}{a character string describing the method.}
#'}
#'
#' @seealso \code{\link{depth.test}}, \code{\link{f.test}}
#'
#' @examples
#' sign.test(rnorm(30))
#' sign.test(y ~ ., data = data.frame(x = rnorm(30), y = rnorm(30)), params = c(1, 1))
#'
#' @rdname signtest
#' @export
sign.test <- function(x, ...) {
  UseMethod("sign.test")
}

#' @rdname signtest
#' @method sign.test default
#' @export
sign.test.default <- function(x, ...) {
  assert_numeric(x, min.len = 2, any.missing = FALSE)

  dname <- deparse(substitute(x))
  method <- "Sign Test (via given residuals)"
  alternative <- "Model where the residuals come from fits not to the data"

  stat <- sum(x > 0)
  n <- length(x)
  names(n) <- "number of observations"
  p.value <- binom.test(stat, n)$p.value
  names(p.value) <- NULL
  names(stat) <- "positive residuals"
  rval <- list(statistic = stat, parameter = n, p.value = p.value,
    alternative = alternative, data.name = dname, method = method)
  class(rval) <- "htest"
  return(rval)
}

#' @rdname signtest
#' @method sign.test formula
#' @export
sign.test.formula <- function(formula, data, params, ...) {
  assert_data_frame(data, any.missing = FALSE, min.cols = 1)
  assert_numeric(params, min.len = 1, any.missing = FALSE)

  mm <- model.matrix(formula, data)
  mr <- model.response(model.frame(formula, data))
  if (ncol(mm) != length(params)) {
    stop(paste("Model matrix has", ncol(mm), "columns and params has length",
      length(params)))
  }

  res <- mm %*% params - mr
  dname <- deparse(substitute(data))
  method <- "Sign Test (via given parameters and data)"
  alternative <- paste("True parameter vector is not equal to",
    paste(params, collapse = " "))
  stat <- sum(res > 0)
  n <- length(res)
  names(n) <- "number of observations"
  p.value <- binom.test(stat, n)$p.value
  names(p.value) <- NULL
  names(stat) <- "positive residuals"
  rval <- list(statistic = stat, parameter = n, p.value = p.value,
    alternative = alternative, data.name = dname, method = method)
  class(rval) <- "htest"
  return(rval)
}
