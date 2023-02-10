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

test_that("ASW ignores Sample and covariable columns", {
  # generate dataset, 3 samples, 5 features
  y <- matrix(rnorm(10*9),3,5)
  y[1,1] <- NA
  y[1:3, 2] <- NA
  y[1:2, 3] <- NA
  # BERT typecasts to dataframe
  y <- data.frame(y)
  y["Batch"] = c(1,1,2)
  y["Label"] = c(1,1,2)
  y2 <- data.frame(y)
  y2["Cov_1"] = 4
  y2["Cov_2"] = 20
  
  no_cat <- compute_asw(y)
  with_cat <- compute_asw(y2)
  
  expect_equal(no_cat$Label, with_cat$Label)
  expect_equal(no_cat$Batch, with_cat$Batch)
})

test_that("Computation of ASW works the same for batch and label", {
  # generate dataset, 3 samples, 5 features
  y <- matrix(rnorm(20*20),20,20)
  # BERT typecasts to dataframe
  y <- data.frame(y)
  y["Batch"] = sample(c(1,2), size=20, replace = TRUE)
  y["Label"] = y[["Batch"]]
  
  asw <- compute_asw(y)
  
  expect_equal(asw$Label, asw$Label)
})


test_that("ASW can also be computed for Batch alone", {
  # generate dataset, 3 samples, 5 features
  y <- matrix(rnorm(20*20),20,20)
  # BERT typecasts to dataframe
  y <- data.frame(y)
  y["Batch"] = sample(c(1,2), size=20, replace = TRUE)
  
  asw <- compute_asw(y)
  
  expect_true(is.na(asw$Label)&&!is.na(asw$Batch))
})

test_that("ASW can also be computed for Label alone", {
  # generate dataset, 3 samples, 5 features
  y <- matrix(rnorm(20*20),20,20)
  # BERT typecasts to dataframe
  y <- data.frame(y)
  y["Label"] = sample(c(1,2), size=20, replace = TRUE)
  
  asw <- compute_asw(y)
  
  expect_true(!is.na(asw$Label)&&is.na(asw$Batch))
})

test_that("Test strip covariable 1", {
  # generate dataset, 3 samples, 5 features
  y <- matrix(rnorm(20*20),20,20)
  # BERT typecasts to dataframe
  y <- data.frame(y)
  y["Cov_1"] = sample(c(1,2), size=20, replace = TRUE)
  
  
  y_nocov = strip_Covariable(y)
  
  expect_true(!("Cov_1" %in% names(y_nocov)))
})


