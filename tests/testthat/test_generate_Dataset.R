test_that("correctly counts numeric values", {
  # generate dataset, 3 samples, 5 features
  y <- matrix(rnorm(10*9),3,5)
  y[1,1] <- NA
  y[1:3, 2] <- NA
  y[1:2, 3] <- NA
  # BERT typecasts to dataframe
  y <- data.frame(y)
  expect_equal(count_existing(y), 9)
})

test_that("Ignores Batch, Label, Sample and covariable columns", {
  # generate dataset, 3 samples, 5 features
  y <- matrix(rnorm(10*9),3,5)
  y[1,1] <- NA
  y[1:3, 2] <- NA
  y[1:2, 3] <- NA
  # BERT typecasts to dataframe
  y <- data.frame(y)
  y["Batch"] = 1
  y["Label"] = 2
  y["Cov_1"] = 4
  y["Cov_2"] = 20
  expect_equal(count_existing(y), 9)
})

