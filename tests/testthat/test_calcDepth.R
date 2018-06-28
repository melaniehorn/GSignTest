library(RobRegTest)
context("calcDepth")

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
  expect_equal(calcDepth(rep(c(-1, 1), 10), 2), qdepth(1, 20, k = 2))
  expect_equal(calcDepth(rep(c(-1, 1), 10), 3), qdepth(1, 20, k = 3))
  expect_equal(calcDepth(rep(c(-1, 1), 10), 4), qdepth(1, 20, k = 4))
  expect_equal(calcDepth(rep(c(-1, 1), 10), 5), qdepth(1, 20, k = 5))
})

test_that("error", {
  expect_error(calcDepth(1, 2))
  expect_error(calcDepth(c(-1, 1), 3))
  expect_error(calcDepth(c(-1, 1, -1), 4))
  expect_error(calcDepth(c(-1, 1, -1, 1), 5))
  expect_error(calcDepth("a", 3))
})
