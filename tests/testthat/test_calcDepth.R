library(GSignTest)
context("calcDepth")

set.seed(0904)

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

test_that("general_old", {
  suppressWarnings(expect_lt(calcDepth_old(rnorm(10), 2), 1))
  suppressWarnings(expect_gte(calcDepth_old(rnorm(10), 2), 0))
  suppressWarnings(expect_lt(calcDepth_old(rnorm(10), 3), 1))
  suppressWarnings(expect_gte(calcDepth_old(rnorm(10), 3), 0))
  suppressWarnings(expect_lt(calcDepth_old(rnorm(10), 4), 1))
  suppressWarnings(expect_gte(calcDepth_old(rnorm(10), 4), 0))
  suppressWarnings(expect_lt(calcDepth_old(rnorm(10), 5), 1))
  suppressWarnings(expect_gte(calcDepth_old(rnorm(10), 5), 0))
})

test_that("general_asymp", {
  expect_lt(calcDepth_asymp(rnorm(10), 2), 1)
  expect_gte(calcDepth_asymp(rnorm(10), 2), 0)
  expect_lt(calcDepth_asymp(rnorm(10), 3), 1)
  expect_gte(calcDepth_asymp(rnorm(10), 3), 0)
  expect_lt(calcDepth_asymp(rnorm(10), 4), 1)
  expect_gte(calcDepth_asymp(rnorm(10), 4), 0)
  expect_lt(calcDepth_asymp(rnorm(10), 5), 1)
  expect_gte(calcDepth_asymp(rnorm(10), 5), 0)
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
  expect_error(calcDepth_old("a", 3))
  expect_error(calcDepth_asymp("a", 3))
})

test_that("equality", {
  res <- rnorm(10)
  expect_equal(calcDepth(res, 2),
    suppressWarnings(calcDepth_old(res, 2)))
  expect_equal(calcDepth(res, 2),
    suppressWarnings(calcDepth_asymp(res, 2)))
  expect_equal(calcDepth(res, 3),
    suppressWarnings(calcDepth_old(res, 3)))
  expect_equal(calcDepth(res, 3),
    suppressWarnings(calcDepth_asymp(res, 3)))
  expect_equal(calcDepth(res, 4),
    suppressWarnings(calcDepth_old(res, 4)))
  expect_equal(calcDepth(res, 5),
    suppressWarnings(calcDepth_old(res, 5)))
})

test_that("warning", {
  expect_warning(calcDepth_old(rnorm(10), 3))
  expect_warning(calcDepth_asymp(rnorm(10), 6))
})
