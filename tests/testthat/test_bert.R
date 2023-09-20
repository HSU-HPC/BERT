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

test_that("BERT allows string input to Reference and not to covariable columns", {
    # generate dataset, 9 samples, 10 features
    y <- matrix(rnorm(10*9),9,10)
    y <- data.frame(y)
    y["Batch"] <- c(1,1,1,1,2,2,2,2,2)
    y["Cov_1"] <- c(1,1,2,2,1,1,1,2,2)
    # numeric covariates
    expect_error(BERT(y, method="limma"), NA)
    # categorical covariates
    y2 = data.frame(y)
    y2$Cov_1 = sapply(y2$Cov_1, as.character)
    expect_error(BERT(y2, method="limma"))
    
    y <- matrix(rnorm(10*9),9,10)
    y <- data.frame(y)
    y["Batch"] <- c(1,1,1,1,2,2,2,2,2)
    y["Reference"] <- c(1,1,2,2,1,1,1,2,2)
    # numeric covariates
    y_cor1 = BERT(y, method="limma")
    y_cor1$Reference = sapply(y_cor1$Reference, as.character)
    # categorical covariates
    y2 = data.frame(y)
    y2$Reference = sapply(y2$Reference, as.character)
    y_cor2 = BERT(y2, method="limma")
    expect_equal(y_cor1, y_cor2)
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
  y <- generate_dataset(100,5,10,0.1,2)
  y_adjusted <- BERT(y, method="None", verify=FALSE)
  expect_true(all.equal(y, y_adjusted[rownames(y), colnames(y)]))
})

test_that("BERT allows the user to specify custom names for batch", {
    y <- generate_dataset(5,5,10,0.1,2)
    y_1 = BERT(y, method="limma")
    # rename
    names(y)[names(y)=="Batch"] = "X"
    y_2 = BERT(y, method="limma", batchname = "X")
    names(y_2)[names(y_2)=="X"] = "Batch"
    expect_true(all.equal(y_1, y_2[rownames(y_2), colnames(y_2)]))
})

test_that("BERT allows the user to specify custom names for label", {
    y <- generate_dataset(5,5,10,0.1,2)
    y_1 = BERT(y, method="limma")
    # rename
    names(y)[names(y)=="Label"] = "X"
    y_2 = BERT(y, method="limma", labelname = "X")
    names(y_2)[names(y_2)=="X"] = "Label"
    expect_true(all.equal(y_1, y_2[rownames(y_2), colnames(y_2)]))
})

test_that("BERT allows the user to specify custom names for references", {
    y <- generate_dataset(5,2,25,0.1,2)
    y["Reference"] = y$Label
    y_1 = BERT(y, method="limma")
    # rename
    names(y)[names(y)=="Reference"] = "X"
    y_2 = BERT(y, method="limma", referencename = "X")
    names(y_2)[names(y_2)=="X"] = "Reference"
    expect_true(all.equal(y_1, y_2[rownames(y_2), colnames(y_2)]))
})

test_that("BERT removes empty columns and still renames everything back", {
    nrows <- 8
    ncols <- 3
    counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
    y = data.frame(counts)
    y[,1] = NA
    y["t"] = c(1,1,1,1,2,2,2,2)
    y2 = BERT(y, batchname = "t")
    expect_true(all.equal(c("X2", "X3", "t"), colnames(y2)))
})

test_that("BERT allows the user to specify custom names for covariables", {
    y <- generate_dataset(5,2,25,0.1,2)
    y["Cov_1"] = y$Label
    y_1 = BERT(y, method="limma")
    # rename
    names(y)[names(y)=="Cov_1"] = "X"
    y_2 = BERT(y, method="limma", covariatename = c("X"))
    names(y_2)[names(y_2)=="X"] = "Cov_1"
    expect_true(all.equal(y_1, y_2[rownames(y_2), colnames(y_2)]))
})

test_that("bert validates all user input -- BERT", {
    # generate dataset, 9 samples, 10 features
    y <- generate_dataset(100,5,10,0.1,2)
    # this should work
    expect_error(validate_bert_input(y, 1, 1, TRUE, FALSE, FALSE, 1, 2,
                                     "file", "None", "X", "B", "R", "S",NULL), NA)
    expect_error(BERT(y, 1, 1,"ComBat", TRUE, FALSE, FALSE, 1, 2, "file", "X",
                      "B", "R", "S", NULL), NA)
    # this should crash
    expect_error(validate_bert_input("blubb", 1, 1, TRUE, FALSE, FALSE, 1, 2,
                                     "file", "None", "X", "B", "R", "S", NULL))
    expect_error(validate_bert_input(y, -1, 1, TRUE, FALSE, FALSE, 1, 2,
                                     "file", "None", "X", "B", "R", "S", NULL))
    expect_error(validate_bert_input(y, 1, 10, TRUE, FALSE, FALSE, 1, 2,
                                     "file", "None", "X", "B", "R", "S", NULL))
    expect_error(validate_bert_input(y, 1, 1, "TRUE", FALSE, FALSE, 1, 2,
                                     "file", "None", "X", "B", "R", "S", NULL))
    expect_error(validate_bert_input(y, 1, 1, TRUE, "", FALSE, 1, 2,
                                     "file", "None", "X", "B", "R", "S", NULL))
    expect_error(validate_bert_input(y, 1, 1, TRUE, FALSE, -10, 1, 2,
                                     "file", "None", "X", "B", "R", "S", NULL))
    expect_error(validate_bert_input(y, 1, 1, TRUE, FALSE, FALSE, FALSE, 2,
                                     "file", "None", "X", "B", "R", "S", NULL))
    expect_error(validate_bert_input(y, 1, 1, TRUE, FALSE, FALSE, 1, "",
                                     "file", "None", "X", "B", "R", "S", NULL))
    expect_error(validate_bert_input(y, 1, 1, TRUE, FALSE, FALSE, 1, 2,
                                     "f", "None", "X", "B", "R", "S", NULL))
    expect_error(validate_bert_input(y, 1, 1, TRUE, FALSE, FALSE, 1, 2,
                                     -1, "None", "X", "B", "R", "S", NULL))
    expect_error(validate_bert_input(y, 1, 1, TRUE, FALSE, FALSE, 1, 2,
                                     "file", FALSE, "X", "B", "R", "S", NULL))
    expect_error(validate_bert_input(y, 1, 1, TRUE, FALSE, FALSE, 1, 2,
                                     "file", "None", -1, "B", "R", "S", NULL))
    expect_error(validate_bert_input(y, 1, 1, TRUE, FALSE, FALSE, 1, 2,
                                     "file", "None", "X", FALSE, "R", "S", NULL))
    expect_error(validate_bert_input(y, 1, 1, TRUE, FALSE, FALSE, 1, 2,
                                     "file", "None", "X", FALSE, "R", "S", c()))
    expect_error(validate_bert_input(y, 1, 1, TRUE, FALSE, FALSE, 1, 2,
                                     "file", "None", "X", FALSE, "R",
                                     "S", c("c", 1)))
    
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