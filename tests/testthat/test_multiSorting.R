library(GSignTest)
context("multiSorting")

test_that("one-dimensional", {
  dat1 <- data.frame(x = 1:10)
  dat2 <- data.frame(x = rnorm(10))

  methods <- c("identity", "random", "norm", "median", "dimSort",
    "weightedSum", "projection", "nonDomSort", "clust", "nn")

  sapply(methods, function(x)
    expect_equal(sort(multiSorting(dat1, method = x)$inds), 1:10))
  sapply(methods, function(x)
    expect_equal(sort(multiSorting(dat2, method = x)$inds), 1:10))
})

test_that("more-dimensional", {
  dat1 <- data.frame(x = 1:10, y = 1:10)
  dat2 <- data.frame(x = rnorm(10), y = rnorm(10))
  dat3 <- as.data.frame(matrix(rnorm(100), ncol = 10))

  methods <- c("identity", "random", "norm", "median", "dimSort",
    "weightedSum", "projection", "nonDomSort", "clust", "nn")

  sapply(methods, function(x)
    expect_equal(sort(multiSorting(dat1, method = x)$inds), 1:10))
  sapply(methods, function(x)
    expect_equal(sort(multiSorting(dat2, method = x)$inds), 1:10))
  sapply(methods, function(x)
    expect_equal(sort(multiSorting(dat3, method = x)$inds), 1:10))
})

test_that("non-default-values", {
  dat1 <- data.frame(x = 1:10, y = 1:10)
  dat2 <- data.frame(x = rnorm(10), y = rnorm(10))
  dat3 <- as.data.frame(matrix(rnorm(100), ncol = 10))

  methods <- c("identity", "random", "norm", "median", "dimSort",
    "weightedSum", "projection", "nonDomSort", "clust", "nn")
  control12 <- list(degree = 1, order = 2:1, weights = c(3, 5), vec.pos = c(1, 2),
    vec.dir = c(2, 1), tie.breaking = "random", clust.method = "single",
    dist.method = "minkowski", p = 4)
  control3 <- control12
  control3$order <- 10:1
  control3$weights <- 1:10
  control3$vec.pos <- 1:10
  control3$vec.dir <- 1:10

  sapply(methods, function(x)
    expect_equal(sort(multiSorting(dat1, method = x, control = control12)$inds), 1:10))
  sapply(methods, function(x)
    expect_equal(sort(multiSorting(dat2, method = x, control = control12)$inds), 1:10))
  sapply(methods, function(x)
    expect_equal(sort(multiSorting(dat3, method = x, control = control3)$inds), 1:10))

  expect_equivalent(multiSorting(dat1, method = "projection", control = list(vec.dir = c(0, 1))),
    multiSorting(dat1, method = "dimSort", control = list(order = 2:1)))
  expect_equivalent(multiSorting(dat2, method = "projection", control = list(vec.dir = c(0, 1))),
    multiSorting(dat2, method = "dimSort", control = list(order = 2:1)))
  expect_equivalent(multiSorting(dat3, method = "projection", control = list(vec.dir = c(rep(0, 9), 1))),
    multiSorting(dat3, method = "dimSort", control = list(order = 10:1)))
  expect_equivalent(multiSorting(dat1, method = "weightedSum", control = list(weights = c(0, 1))),
    multiSorting(dat1, method = "dimSort", control = list(order = 2:1)))
  expect_equivalent(multiSorting(dat2, method = "weightedSum", control = list(weights = c(0, 1))),
    multiSorting(dat2, method = "dimSort", control = list(order = 2:1)))
  expect_equivalent(multiSorting(dat3, method = "weightedSum", control = list(weights = c(rep(0, 9), 1))),
    multiSorting(dat3, method = "dimSort", control = list(order = 10:1)))
})

test_that("error", {
  dat1 <- data.frame(x = 1:10, y = 1:10)
  expect_error(multiSorting(data.frame(x = 1, y = 2)))
  expect_error(multiSorting(dat1, method = "bla"))
  expect_error(multiSorting(dat1, method = "norm", control = list(degree = -1)))
  expect_error(multiSorting(dat1, method = "dimSort", control = list(order = 1)))
  expect_error(multiSorting(dat1, method = "weightedSum", control = list(weights = c(0, 0, 0))))
  expect_error(multiSorting(dat1, method = "projection", control = list(vec.pos = c(0, 0, 0))))
  expect_error(multiSorting(dat1, method = "nonDomSort", control = list(tie.breaking = "nonDomSort")))
  expect_error(multiSorting(dat1, method = "clust", control = list(clust.method = "bla")))
  expect_error(multiSorting(dat1, method = "nn", control = list(dist.method = "bla")))
  expect_error(suppressWarnings(multiSorting(dat1, method = "shp", control = list(path = "bla"))))
})
