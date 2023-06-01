test_that("dataset chunking identifies correct number of chunks",{
  y <- matrix(rnorm(10*12),12,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3,4,4,4)
  chunks <- chunk_data(y, 2, "file")
  # correct number of chunks
  expect_equal(length(chunks),2)
  # max number of chunks is n_batches//2
  chunks <- chunk_data(y, 3)
  expect_equal(length(chunks),2)
})

test_that("dataset chunks are readable and have correct shape",{
  y <- matrix(rnorm(10*12),12,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3,4,4,4)
  chunks <- chunk_data(y, 2, "file")
  for(c in chunks){
    data = readRDS(c)
    shape = dim(data)
    expect_equal(shape[1],6)
    expect_equal(shape[2],11)
  }
})

test_that("dataset chunks are read as dataframes",{
  y <- matrix(rnorm(10*12),12,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3,4,4,4)
  chunks <- chunk_data(y, 2, "file")
  for(c in chunks){
    data = readRDS(c)
    expect_equal(typeof(data), "list")
  }
})

test_that("chunking can be configured to return dataframes",{
  y <- matrix(rnorm(10*12),12,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3,4,4,4)
  chunks <- chunk_data(y, 2, "default")
  for(c in chunks){
    expect_equal(typeof(c), "list")
  }
})

test_that("chunking the dataset to one split works",{
  y <- matrix(rnorm(10*12),12,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3,4,4,4)
  chunks <- chunk_data(y, 1, "file")
  for(c in chunks){
    data = readRDS(c)
    expect_equal(dim(data)[1],12)
    expect_equal(dim(data)[2],11)
  }
})

test_that("only one batch remains after subtree adjustment -- parallel_bert",{
  # generate data
  y <- matrix(rnorm(10*15),15,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5)
  chunks <- chunk_data(y, 1)
  cl <- parallel::makeCluster(1)
  doParallel::registerDoParallel(cl)
  adjusteddata <- parallel_bert(chunks, 1, method = "None")
  parallel::stopCluster(cl)
  num_batches = length(unique(adjusteddata$Batch))
  expect_equal(num_batches, 1)
})

test_that("parallel_bert works with file communication backend",{
  # generate data
  y <- matrix(rnorm(10*15),15,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5)
  chunks <- chunk_data(y, 1, "file")
  cl <- parallel::makeCluster(1)
  doParallel::registerDoParallel(cl)
  adjusteddata <- parallel_bert(chunks, 1, method = "None", "file")
  parallel::stopCluster(cl)
  num_batches = length(unique(adjusteddata$Batch))
  expect_equal(num_batches, 1)
})


test_that("dataset not corrupted without adjustment -- BERT", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y_adjusted <- BERT(y, method="None")
  expect_true(all.equal(y, y_adjusted[rownames(y), colnames(y)]))
})

test_that("BERT does not allow combination of covariable and reference columns", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y["Cov_1"] <- c(1,1,2,2,3,3,4,4,5)
  y["Reference"] <- c(1,1,2,2,3,3,4,4,5)
  # sequential
  testthat::expect_error(BERT(y, method="None"))
  # parallel
  testthat::expect_error(BERT(y, 2, method="None"))
})

test_that("BERT likes SummarizedExperiments",{
  nrows <- 200
  ncols <- 8
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  colData <- data.frame(Batch=c(1,1,1,1,2,2,2,2))
  y = SummarizedExperiment::SummarizedExperiment(assays=list(counts=counts), colData=colData)
  expect_error(BERT(y), NA)
})


test_that("BERT preserves order of samples and features", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  names(y) = c("A","B","C","D","E","F","G","H","I","J")
  rownames(y) = c("F1","F2","F3","F4","F5","F6","F7","F8","F9")
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  # sequential
  y_adj = BERT(y)
  expect_true(all(names(y)==names(y_adj)))
  rownames(all(rownames(y)==rownames(y_adj)))
})

test_that("bert works for simulated data without any formatting of the input -- BERT", {
  # generate dataset, 9 samples, 10 features
  y <- generateDataset(100,5,10,0.1,2)
  y_adjusted <- BERT(y, method="None", verify=FALSE)
  expect_true(all.equal(y, y_adjusted[rownames(y), colnames(y)]))
})


test_that("works equally with dataframes and matrices -- BERT", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y_adjusted_df <- BERT(y, method="None")
  y_adjusted_mat <- BERT(as.matrix(y), method="None")
  expect_true(all.equal(y_adjusted_df, y_adjusted_mat[rownames(y_adjusted_df), colnames(y_adjusted_df)]))
})

test_that("BERT preserves rownames and column names -- BERT", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y <- as.matrix(y)
  rownames(y) <- c("A","B","C","D","E","F","G","H","I")
  colnames(y) <- c("A1","B2","C3","D4","E5","F6","G7","H8","I9","J10", "Batch")
  y_adjusted_df <- BERT(y, method="None")
  expect_true(all.equal(rownames(y), rownames(y_adjusted_df)))
  expect_true(all.equal(colnames(y), colnames(y_adjusted_df)))
})


test_that("BERT prints ASW Batch if qualitycontrol is TRUE -- BERT", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y <- as.matrix(y)
  rownames(y) <- c("A","B","C","D","E","F","G","H","I")
  colnames(y) <- c("A1","B2","C3","D4","E5","F6","G7","H8","I9","J10", "Batch")
  expect_output(BERT(y, method="None", qualitycontrol = TRUE), "\\w*ASW Batch was\\w*")
})


test_that("BERT prints ASW Label if qualitycontrol is TRUE -- BERT", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y["Label"] <- c(0,1,2,0,1,2,0,1,2)
  y <- as.matrix(y)
  rownames(y) <- c("A","B","C","D","E","F","G","H","I")
  colnames(y) <- c("A1","B2","C3","D4","E5","F6","G7","H8","I9","J10", "Batch", "Label")
  expect_output(BERT(y, method="None", qualitycontrol = TRUE), "\\w*ASW Label was\\w*")
})

test_that("BERT works with references", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y[1:3,1:10] <- y[1:3,1:10] + 3
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y["Label"] <- c(0,1,2,0,1,2,0,1,2)
  y["Reference"] <- c(1,2,0,1,2,0,1,2,0)
  expect_error(BERT(y, method="ref"), NA)
})