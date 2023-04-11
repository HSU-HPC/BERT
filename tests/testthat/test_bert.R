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


test_that("dataset not corrupted without adjustment -- hierarchical_adjustment", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y_adjusted <- hierarchical_adjustment(y, method="None")
  expect_true(all.equal(y, y_adjusted[rownames(y), colnames(y)]))
})

test_that("dataset not corrupted without adjustment -- BERT", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y_adjusted <- hierarchical_adjustment(y, method="None")
  expect_true(all.equal(y, y_adjusted[rownames(y), colnames(y)]))
})

test_that("bert works for simulated data without any formatting of the input -- hierarchical_adjustment", {
  # generate dataset, 9 samples, 10 features
  y <- generateDataset(100,5,10,0.1,2)
  y_adjusted <- hierarchical_adjustment(y, method="None", verify=FALSE)
  expect_true(all.equal(y, y_adjusted[rownames(y), colnames(y)]))
})

test_that("bert works for simulated data without any formatting of the input -- BERT", {
  # generate dataset, 9 samples, 10 features
  y <- generateDataset(100,5,10,0.1,2)
  y_adjusted <- hierarchical_adjustment(y, method="None", verify=FALSE)
  expect_true(all.equal(y, y_adjusted[rownames(y), colnames(y)]))
})

test_that("works equally with dataframes and matrices -- hierarchical_adjustment", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y_adjusted_df <- hierarchical_adjustment(y, method="None")
  y_adjusted_mat <- hierarchical_adjustment(as.matrix(y), method="None")
  expect_true(all.equal(y_adjusted_df, y_adjusted_mat[rownames(y_adjusted_df), colnames(y_adjusted_df)]))
})

test_that("works equally with dataframes and matrices -- BERT", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y_adjusted_df <- hierarchical_adjustment(y, method="None")
  y_adjusted_mat <- hierarchical_adjustment(as.matrix(y), method="None")
  expect_true(all.equal(y_adjusted_df, y_adjusted_mat[rownames(y_adjusted_df), colnames(y_adjusted_df)]))
})

test_that("BERT preserves rownames and column names -- hierarchical adjustment", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y <- as.matrix(y)
  rownames(y) <- c("A","B","C","D","E","F","G","H","I")
  colnames(y) <- c("A1","B2","C3","D4","E5","F6","G7","H8","I9","J10", "Batch")
  y_adjusted_df <- hierarchical_adjustment(y, method="None")
  expect_true(all.equal(rownames(y), rownames(y_adjusted_df)))
  expect_true(all.equal(colnames(y), colnames(y_adjusted_df)))
})

test_that("BERT preserves rownames and column names -- BERT", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y <- as.matrix(y)
  rownames(y) <- c("A","B","C","D","E","F","G","H","I")
  colnames(y) <- c("A1","B2","C3","D4","E5","F6","G7","H8","I9","J10", "Batch")
  y_adjusted_df <- hierarchical_adjustment(y, method="None")
  expect_true(all.equal(rownames(y), rownames(y_adjusted_df)))
  expect_true(all.equal(colnames(y), colnames(y_adjusted_df)))
})

test_that("BERT prints ASW Batch if qualitycontrol is TRUE -- hierarchical_adjustment", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y <- as.matrix(y)
  rownames(y) <- c("A","B","C","D","E","F","G","H","I")
  colnames(y) <- c("A1","B2","C3","D4","E5","F6","G7","H8","I9","J10", "Batch")
  expect_output(hierarchical_adjustment(y, method="None", qualitycontrol = TRUE), "\\w*ASW Batch was\\w*")
})

test_that("BERT prints ASW Batch if qualitycontrol is TRUE -- BERT", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y <- as.matrix(y)
  rownames(y) <- c("A","B","C","D","E","F","G","H","I")
  colnames(y) <- c("A1","B2","C3","D4","E5","F6","G7","H8","I9","J10", "Batch")
  expect_output(hierarchical_adjustment(y, method="None", qualitycontrol = TRUE), "\\w*ASW Batch was\\w*")
})

test_that("BERT prints ASW Label if qualitycontrol is TRUE -- hierarchical adjustment", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y["Label"] <- c(0,1,2,0,1,2,0,1,2)
  y <- as.matrix(y)
  rownames(y) <- c("A","B","C","D","E","F","G","H","I")
  colnames(y) <- c("A1","B2","C3","D4","E5","F6","G7","H8","I9","J10", "Batch", "Label")
  expect_output(hierarchical_adjustment(y, method="None", qualitycontrol = TRUE), "\\w*ASW Label was\\w*")
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
  expect_output(hierarchical_adjustment(y, method="None", qualitycontrol = TRUE), "\\w*ASW Label was\\w*")
})