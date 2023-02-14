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
  expect_true(all.equal(format_DF(y)[["Batch"]], ordinal_encode(y[["Batch"]])))
})

test_that("works with categorical label", {
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y["Label"] <- c("A", "A", "A", "B", "B", "B", "C", "C", "C")
  expect_error(format_DF(y), NA)
  expect_true(all.equal(format_DF(y)[["Label"]], ordinal_encode(y[["Label"]])))
})

test_that("works with categorical covariate", {
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y["Cov_3"] <- c("A", "A", "A", "B", "B", "B", "C", "C", "C")
  expect_error(format_DF(y), NA)
  expect_true(all.equal(format_DF(y)[["Cov_3"]], ordinal_encode(y[["Cov_3"]])))
})

test_that("ordinal encoding works", {
  col = c("A", "A", "B", "B")
  enc <- ordinal_encode(col)
  expected <- c(1,1,2,2)
  expect_true(all.equal(expected, enc))
})

test_that("NaNs are replaced correctly", {
  mat <- matrix(rnorm(5*20), nrow = 5, ncol=20)
  mat <- data.frame(mat)
  mat[5, 7] <- NaN
  mat_rep <- replace_missing(mat)
  expect_true(is.na(mat_rep[5,7]))
  expect_true(!is.nan(mat_rep[5,7]))
  expect_equal(sum(is.na(mat_rep)), 1)
})

test_that("Formatting calls NaN replacement function", {
  mat <- matrix(rnorm(5*20), nrow = 5, ncol=20)
  mat <- data.frame(mat)
  mat["Batch"] <- c(1,1,2,2,2)
  mat[5, 7] <- NaN
  mat_rep <- format_DF(mat)
  expect_true(is.na(mat_rep[5,7]))
  expect_true(!is.nan(mat_rep[5,7]))
  expect_equal(sum(is.na(mat_rep)), 1)
})


test_that("Formatting will use get_adjustable_features_with_mod, if covariables are present", {
  mat <- matrix(rnorm(5*5), nrow=5, ncol=5)
  mat <- data.frame(mat)
  mat["Batch"] <- c(1,1,1,1,1)
  mat["Cov_1"] <- c(1,1,1,2,2)
  mat[1,4] <- NA
  mat[4,1] <- NA
  formatted_df <- format_DF(mat)
  expect_equal(all(is.na(formatted_df[,1])), TRUE)
})