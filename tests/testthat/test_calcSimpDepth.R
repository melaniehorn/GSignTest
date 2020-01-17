library(GSignTest)
context("calcSimpDepth")

test_that("general", {
  expect_lte(calcSimpDepth(rnorm(10), 2), 1)
  expect_gte(calcSimpDepth(rnorm(10), 2), 0)
  expect_lte(calcSimpDepth(rnorm(10), 3), 1)
  expect_gte(calcSimpDepth(rnorm(10), 3), 0)
  expect_lte(calcSimpDepth(rnorm(10), 4), 1)
  expect_gte(calcSimpDepth(rnorm(10), 4), 0)
  expect_lte(calcSimpDepth(rnorm(10), 5), 1)
  expect_gte(calcSimpDepth(rnorm(10), 5), 0)
})

test_that("worst", {
  expect_equal(calcSimpDepth(rep(c(-1, 1), 10), 2), 1)
  expect_equal(calcSimpDepth(rep(c(-1, 1), 10), 3), 1)
  expect_equal(calcSimpDepth(rep(c(-1, 1), 10), 4), 1)
  expect_equal(calcSimpDepth(rep(c(-1, 1), 10), 5), 1)
})

test_that("error", {
  expect_error(calcSimpDepth(1, 2))
  expect_error(calcSimpDepth(c(-1, 1), 3))
  expect_error(calcSimpDepth(c(-1, 1, -1), 4))
  expect_error(calcSimpDepth(c(-1, 1, -1, 1), 5))
  #expect_error(calcSimpDepth("a", 3))
})
