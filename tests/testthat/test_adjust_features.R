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

test_that("selection of adjustable features works with covariables (1 covariate)", {
  mat <- matrix(rnorm(5*5), nrow=5, ncol=5)
  mat <- data.frame(mat)
  mat["Batch"] <- c(1,1,1,1,1)
  mat["Cov_1"] <- c(1,1,1,2,2)
  mat["Cov_2"] <- c(1,1,1,2,2)
  mat[1,4] <- NA
  mat[4,1] <- NA
  mod <- data.frame(mat [ , grepl( "Cov" , names( mat  ) ) ])
  mat <- mat [ , !grepl( "Cov" , names( mat  ) ) ]
  adjustable <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
  computed_adjustable <- as.vector(get_adjustable_features_with_mod(mat, mod))
  expect_true(all.equal(adjustable, computed_adjustable))
})

test_that("selection of adjustable features works with covariables (1 covariate)", {
  mat <- matrix(rnorm(5*5), nrow=5, ncol=5)
  mat <- data.frame(mat)
  mat["Batch"] <- c(1,1,1,1,1)
  mat["Cov_1"] <- c(1,1,1,2,2)
  mat[1,4] <- NA
  mat[4,1] <- NA
  mod <- data.frame(mat [ , grepl( "Cov" , names( mat  ) ) ])
  mat <- mat [ , !grepl( "Cov" , names( mat  ) ) ]
  adjustable <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
  computed_adjustable <- as.vector(get_adjustable_features_with_mod(mat, mod))
  expect_true(all.equal(adjustable, computed_adjustable))
})

test_that("BERT stops, if not enough samples from batch/covariate level", {
  mat <- matrix(rnorm(1*5), nrow=1, ncol=5)
  mat <- data.frame(mat)
  mat["Batch"] <- c(1)
  expect_error(get_adjustable_features(mat))
})