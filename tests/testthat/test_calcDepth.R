library(GSignTest)
context("calcDepth")

set.seed(2306)

test_that("general", {
  expect_lt(calcDepth(rnorm(10), 2), 1)
  expect_gte(calcDepth(rnorm(10), 2), 0)
  expect_lt(calcDepth(rnorm(10), 3), 1)
  expect_gte(calcDepth(rnorm(10), 3), 0)
  expect_lt(calcDepth(rnorm(10), 4), 1)
  expect_gte(calcDepth(rnorm(10), 4), 0)
  expect_lt(calcDepth(rnorm(10), 5), 1)
  expect_gte(calcDepth(rnorm(10), 5), 0)
})

test_that("worst", {
  expect_equal(calcDepth(rep(c(-1, 1), 10), 2), qdepth(1, 20, K = 2))
  expect_equal(calcDepth(rep(c(-1, 1), 10), 3), qdepth(1, 20, K = 3))
  expect_equal(calcDepth(rep(c(-1, 1), 10), 4), qdepth(1, 20, K = 4))
  expect_equal(calcDepth(rep(c(-1, 1), 10), 5), qdepth(1, 20, K = 5))
})

test_that("error", {
  expect_error(calcDepth(1, 2))
  expect_error(calcDepth(c(-1, 1), 3))
  expect_error(calcDepth(c(-1, 1, -1), 4))
  expect_error(calcDepth(c(-1, 1, -1, 1), 5))
  expect_error(calcDepth("a", 3))
})

test_that("equality", {
  res <- rnorm(10)
  expect_equal(calcDepth(res, 2),
    suppressWarnings(calcDepth(res, 2, linear = FALSE)))
  expect_equal(calcDepth(res, 3),
    suppressWarnings(calcDepth(res, 3, linear = FALSE)))
  expect_equal(calcDepth(res, 4),
    suppressWarnings(calcDepth(res, 4, linear = FALSE)))
  expect_equal(calcDepth(res, 5),
    suppressWarnings(calcDepth(res, 5, linear = FALSE)))
})

test_that("warning", {
  expect_warning(calcDepth(rnorm(10), 3, linear = FALSE))
  expect_warning(calcDepth(rnorm(10), 6))
})
