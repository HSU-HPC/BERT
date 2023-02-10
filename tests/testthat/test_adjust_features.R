test_that("correctly identifies adjustable features", {
  # generate dataset, 3 samples, 5 features
  y <- matrix(rnorm(3*5),3,5)
  y[1,1] <- NA
  y[1:3, 2] <- NA
  y[1:2, 3] <- NA
  # BERT typecasts to dataframe
  y <- data.frame(y)
  adjustable = c(TRUE, FALSE, FALSE, TRUE, TRUE)
  expect_true(all.equal(adjustable, as.vector(get_adjustable_features(y))))
})

test_that("only combat modes 1-4 are allowed", {
  # generate dataset, 3 samples, 5 features
  y <- matrix(rnorm(10*6),6,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2)
  # no covariates
  mod = data.frame(matrix(NA, nrow=length(rownames(y)), ncol=0))
  modes = c(1,2,3,4)
  for(m in modes){
    adjust_node(y, 1, 2, mod, m, "ComBat")
  }
  # should crash
  expect_error(adjust_node(y, 1, 2, mod, "invalid mode", "ComBat"))
})

