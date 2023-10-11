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

test_that("correctly formatted Reference column is recognized as such",{
  y <- matrix(rnorm(10*9),10,9)
  y <- data.frame(y)
  y["Reference"] <- c(1,1,2,2,0,1,1,2,2,0)
  #y["Batch"] <- c(1,1,1,1,1,0,0,0,0,0)
  expect_true(verify_references(y))
})

test_that("correctly missing Reference column is ok",{
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  expect_true(verify_references(y))
})

test_that("MVs in reference columns are recognized",{
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Reference"] <- c(1,NA,1,2,1,0,1,1,1)
  expect_true(!verify_references(y))
})

test_that("recognizes too low number of references",{
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Reference"] <- c(0,0,0,0,0,0,0,0,1)
  expect_true(!verify_references(y))
})

test_that("casts SummarizedExperiment to dataframe",{
  nrows <- 200
  ncols <- 8
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  colData <- data.frame(Batch=c(1,1,1,1,2,2,2,2))
  y = SummarizedExperiment::SummarizedExperiment(assays=list(counts=counts), colData=colData)
  y_df = format_DF(y, assayname = "counts")
  # dimension correct
  dimension = dim(y_df)
  expect_true(all(dimension==c(8,201)))
  
  # label, sample and references can be appended
  colData <- data.frame(Batch=c(1,1,1,1,2,2,2,2), "Reference"=c(1,1,0,0,1,1,0,0), "Label"=c(1,2,1,2,1,2,1,2), "Sample"=c(1,2,3,4,5,6,7,8))
  y = SummarizedExperiment::SummarizedExperiment(assays=list(counts=counts), colData=colData)
  y_df = format_DF(y, assayname = "counts")
  # dimension correct
  dimension = dim(y_df)
  expect_true(all(dimension==c(8,204)))
  
  # check for covariables
  nrows <- 10
  ncols <- 8
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  colData <- data.frame(Batch=c(1,1,1,1,2,2,2,2), "Label"=c(1,2,1,2,1,2,1,2), "Sample"=c(1,2,3,4,5,6,7,8), "Cov_1"=c(1,1,2,2,1,1,2,2))
  y = SummarizedExperiment::SummarizedExperiment(assays=list(counts=counts), colData=colData)
  y_df = format_DF(y, assayname = "counts")
  # dimension correct
  dimension = dim(y_df)
  expect_true(all(dimension==c(8,14)))
  
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

test_that("Formatting will use get_adjustable_features, if no covariables are present", {
  mat <- matrix(rnorm(5*5), nrow=5, ncol=5)
  mat <- data.frame(mat)
  mat["Batch"] <- c(1,1,1,1,1)
  mat[1,4] <- NA
  mat[4,1] <- NA
  formatted_df <- format_DF(mat)
  expect_true(all.equal(mat, formatted_df))
})

test_that("Formatting preserves names, if input was matrix", {
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y <- as.matrix(y)
  rownames(y) <- c("A","B","C","D","E","F","G","H","I")
  colnames(y) <- c("A1","B2","C3","D4","E5","F6","G7","H8","I9","J10", "Batch")
  formatted_df <- format_DF(y)
  expect_true(all.equal(rownames(y), rownames(formatted_df)))
  expect_true(all.equal(colnames(y), colnames(formatted_df)))
})