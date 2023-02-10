test_that("correctly identifies adjustable features", {
  # generate dataset, 3 samples, 5 features
  y <- matrix(rnorm(10*9),3,5)
  y[1,1] <- NA
  y[1:3, 2] <- NA
  y[1:2, 3] <- NA
  # BERT typecasts to dataframe
  y <- data.frame(y)
  adjustable = c(TRUE, FALSE, FALSE, TRUE, TRUE)
  expect_true(all.equal(adjustable, as.vector(get_adjustable_features(y))))
})

