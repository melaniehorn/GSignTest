library(GSignTest)
context("ftest")

test_that("formula", {
  dat1 <- data.frame(x = rnorm(20), y = rnorm(20), z = rnorm(20))
  expect_lte(f.test(z ~ ., data = dat1, params = c(1, 1, 1))$p.value, 1)
  expect_gte(f.test(z ~ ., data = dat1, params = c(1, 1, 1))$p.value, 0)

  expect_error(f.test(z ~ ., data = dat1, params = c(1, 1)))
})
