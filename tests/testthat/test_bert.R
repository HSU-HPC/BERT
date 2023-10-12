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
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  adjusteddata <- parallel_bert(chunks, method = "None")
  num_batches = length(unique(adjusteddata$Batch))
  expect_equal(num_batches, 1)
})

test_that("parallel_bert works equivalently in sequential or parallel mode",{
    # generate data
    y <- generate_dataset(50,10,10,0.1,2)
    chunks <- chunk_data(y, 2)
    BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
    adjusteddata_sequential <- parallel_bert(chunks, method = "None")
    BiocParallel::register(BiocParallel::SnowParam(workers = 2), default = TRUE)
    adjusteddata_parallel <- parallel_bert(chunks, method = "None")
    expect_true(all(adjusteddata_sequential==adjusteddata_parallel, na.rm = TRUE))
})

test_that("parallel_bert has correct rownames",{
    # generate data
    y <- generate_dataset(50,10,10,0.1,2)
    chunks <- chunk_data(y, 2)
    adjusteddata_sequential <- parallel_bert(chunks, method = "None")
    expect_true(all(y[colnames(y)[colnames(y)!="Batch"]]==adjusteddata_sequential[row.names(y),][colnames(y)[colnames(y)!="Batch"]], na.rm = TRUE))
})

test_that("parallel_bert works with file communication backend",{
  # generate data
  y <- matrix(rnorm(10*15),15,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5)
  chunks <- chunk_data(y, 1, "file")
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  adjusteddata <- parallel_bert(chunks, 1, method = "None", "file")
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
  expect_error(BERT(y, assayname = "counts"), NA)
})

test_that("BERT writes new assay to SummarizedExperiments",{
    nrows <- 200
    ncols <- 8
    counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
    colData <- data.frame(Batch=c(1,1,1,1,2,2,2,2))
    y = SummarizedExperiment::SummarizedExperiment(assays=list(counts=counts), colData=colData)
    y_corrected = BERT(y, assayname = "counts", method="None")
    raw_mat = SummarizedExperiment::assays(y)$counts
    proc_mat = SummarizedExperiment::assays(y)$counts_BERTcorrected
    expect_true(all(raw_mat=proc_mat))
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
    # this should not work
    expect_error(BERT("string"))
    expect_error(BERT(y, cores=-1))
    expect_error(BERT(y, cores="hui"))
    expect_error(BERT(y, combatmode = -1))
    expect_error(BERT(y, combatmode = "1"))
    expect_error(BERT(y, corereduction = "1"))
    expect_error(BERT(y, stopParBatches = "-1"))
    expect_error(BERT(y, backend = 1))
    expect_error(BERT(y, method = TRUE))
    expect_error(BERT(y, qualitycontrol = -1))
    expect_error(BERT(y, qualitycontrol = ""))
    expect_error(BERT(y, verify = -1))
    expect_error(BERT(y, verify = ""))
    expect_error(BERT(y, labelname = ""))
    expect_error(BERT(y, labelname = 1))
    expect_error(BERT(y, batchname = ""))
    expect_error(BERT(y, batchname = 1))
    expect_error(BERT(y, referencename = ""))
    expect_error(BERT(y, referencename = 1))
    expect_error(BERT(y, samplename = ""))
    expect_error(BERT(y, samplename = 1))
    expect_error(BERT(y, covariatename = ""))
    expect_error(BERT(y, covariatename = 1))
    expect_error(BERT(y, assayname = ""))
    expect_error(BERT(y, assayname = 1))
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

test_that("BERT works with both default backend (from bpparam) and with user-defined number of cores", {
    ds = generate_dataset(100,10,10,0.1,2)
    y_default = BERT(ds, method = "None")
    y_self = BERT(ds, method = "None", cores = 2)
    expect_true(all(y_default==y_self, na.rm = TRUE))
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