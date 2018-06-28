
#' @title F Test
#'
#' @description
#' Performs a test whether a linear model fits to some data via a regression F Test.
#'
#' @param formula [\code{formula}]\cr a formula of the form \code{lhs ~ rhs} where
#' lhs is a numeric giving the target variable and rhs are the influencing variables .
#' @param data [\code{data.frame}]\cr a data frame containing the variables in the formula.
#' @param params [\code{numeric}]\cr a parameter vector for the null-hypotheses.
#'
#' @details
#' Here, a classical F Test for testing the parameter vector of a linear regression
#' model is done. Only the two sided test is currently supported.
#'
#' @return A list with class \code{"htest"} containing the following components:
#'\describe{
#'  \item{\code{statistic}}{the F Test statistic.}
#'  \item{\code{parameter}}{the degrees of freedom of the appropriate F distribution.}
#'  \item{\code{p.value}}{the p-value for the test.}
#'  \item{\code{alternative}}{a character string describing the alternative hypothesis.}
#'  \item{\code{data.name}}{a character string giving the name of the data.}
#'  \item{\code{method}}{a character string describing the method.}
#'}
#'
#' @seealso \code{\link{depth.test}}, \code{\link{sign.test}}
#'
#' @examples
#' f.test(y ~ ., data = data.frame(x = rnorm(30), y = rnorm(30)), params = c(1, 1))
#'
#' @rdname ftest
#' @export
f.test <- function(formula, data, params) {
  assert_data_frame(data, any.missing = FALSE, min.cols = 1)
  assert_numeric(params, min.len = 1, any.missing = FALSE)

  mm <- model.matrix(formula, data)
  if (ncol(mm) != length(params)) {
    stop(paste("Model matrix has", ncol(mm), "columns and params has length",
      length(params)))
  }

  mod <- lm(formula, data)
  res <- mod$residuals
  coef <- mod$coefficients

  dname <- deparse(substitute(data))
  method <- "F Test (via given parameters and data)"
  alternative <- paste("True parameter vector is not equal to",
    paste(params, collapse = " "))
  n <- length(res)
  k <- length(params)

  num <- t(coef - params) %*% (t(mm) %*% mm) %*% (coef - params) * (n - k)
  denom <- sum(res^2) * k
  stat <- num / denom

  q <- pf(stat, k, n - k)
  p.value <- 2 * min(q, 1 - q)
  names(stat) <- "F"
  parameter <- c(num.df = k, denom.df = n - k)
  rval <- list(statistic = stat, parameter = parameter, p.value = p.value,
    alternative = alternative, data.name = dname, method = method)
  class(rval) <- "htest"
  return(rval)
}
