test_that("removes unadjustable numeric values", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  # make first feature unadjustable in batch 1
  y[1:2,1] = NA
  y_formatted <- format_DF(y)
  expect_true(is.na(y_formatted[3,1]))
})

test_that("removes empty columns", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  # make first feature unadjustable in batch 1
  y[,1] = NA
  y_formatted <- format_DF(y)
  expect_equal(dim(y_formatted)[2], 10)
})

test_that("removes empty columns", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  # make first feature unadjustable in batch 1
  y[1,] = NA
  y_formatted <- format_DF(y)
  expect_equal(dim(y_formatted)[1], 8)
})

test_that("works with categorical batch", {
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c("A", "A", "A", "B", "B", "B", "C", "C", "C")
  expect_error(format_DF(y), NA)
})

test_that("works with categorical label", {
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y["Label"] <- c("A", "A", "A", "B", "B", "B", "C", "C", "C")
  expect_error(format_DF(y), NA)
})

test_that("works with categorical covariate", {
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y["Cov_3"] <- c("A", "A", "A", "B", "B", "B", "C", "C", "C")
  expect_error(format_DF(y), NA)
})

test_that("ordinal encoding works", {
  col = c("A", "A", "B", "B")
  enc <- ordinal_encode(col)
  expected <- c(1,1,2,2)
  expect_true(all.equal(expected, enc))
})

