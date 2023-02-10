test_that("dataset not corrupted without adjustment", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y_adjusted <- hierarchical_adjustment(y, method="None")
  expect_true(all.equal(y, y_adjusted[rownames(y), colnames(y)]))
})

test_that("works equally with dataframes and matrices", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y_adjusted_df <- hierarchical_adjustment(y, method="None")
  y_adjusted_mat <- hierarchical_adjustment(as.matrix(y), method="None")
  expect_true(all.equal(y_adjusted_df, y_adjusted_mat[rownames(y_adjusted_df), colnames(y_adjusted_df)]))
})