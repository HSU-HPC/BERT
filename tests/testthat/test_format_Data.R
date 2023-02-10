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