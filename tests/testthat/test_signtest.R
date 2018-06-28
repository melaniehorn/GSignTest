library(RobRegTest)
context("signtest")

test_that("default", {
  res1 <- rep(c(-1, 1), 10)
  res2 <- rnorm(20)

  expect_s3_class(sign.test(res1), "htest")
  expect_s3_class(sign.test(res2), "htest")

  expect_equal(sign.test(res1)$p.value, 1)
  expect_lte(sign.test(res2)$p.value, 1)
  expect_gte(sign.test(res2)$p.value, 0)

  expect_error(depth.test("a"))
})

test_that("formula", {
  dat1 <- data.frame(x = rnorm(20), y = rnorm(20), z = rnorm(20))
  expect_s3_class(sign.test(z ~ ., data = dat1, params = c(1, 1, 1)), "htest")
  expect_lte(sign.test(z ~ ., data = dat1, params = c(1, 1, 1))$p.value, 1)
  expect_gte(sign.test(z ~ ., data = dat1, params = c(1, 1, 1))$p.value, 0)
})
