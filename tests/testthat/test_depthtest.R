library(RobRegTest)
context("depthtest")

test_that("default", {
  res1 <- rep(c(-1, 1), 10)
  res2 <- rnorm(20)

  expect_s3_class(depth.test(res1, 3), "htest")
  expect_s3_class(depth.test(res2, 3), "htest")

  expect_equal(depth.test(res1, 3)$p.value, 1)
  expect_lte(depth.test(res2, 3)$p.value, 1)
  expect_gte(depth.test(res2, 3)$p.value, 0)

  expect_error(depth.test(res1, 21))
  expect_error(depth.test(res2, 21))
})

test_that("formula", {
  dat1 <- data.frame(x = rnorm(20), y = rnorm(20), z = rnorm(20))
  expect_s3_class(depth.test(z ~ ., data = dat1, params = c(1, 1, 1), k = 3), "htest")
  expect_lte(depth.test(z ~ ., data = dat1, params = c(1, 1, 1), k = 3)$p.value, 1)
  expect_gte(depth.test(z ~ ., data = dat1, params = c(1, 1, 1), k = 3)$p.value, 0)
})
