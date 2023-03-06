test_that("dataset not corrupted without adjustment", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y_adjusted <- hierarchical_adjustment(y, method="None")
  expect_true(all.equal(y, y_adjusted[rownames(y), colnames(y)]))
})

test_that("works equally with dataframes and matrices", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y_adjusted_df <- hierarchical_adjustment(y, method="None")
  y_adjusted_mat <- hierarchical_adjustment(as.matrix(y), method="None")
  expect_true(all.equal(y_adjusted_df, y_adjusted_mat[rownames(y_adjusted_df), colnames(y_adjusted_df)]))
})

test_that("BERT preserves rownames and column names", {
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

test_that("BERT prints ASW Batch if qualitycontrol is TRUE", {
  # generate dataset, 9 samples, 10 features
  y <- matrix(rnorm(10*9),9,10)
  y <- data.frame(y)
  y["Batch"] <- c(1,1,1,2,2,2,3,3,3)
  y <- as.matrix(y)
  rownames(y) <- c("A","B","C","D","E","F","G","H","I")
  colnames(y) <- c("A1","B2","C3","D4","E5","F6","G7","H8","I9","J10", "Batch")
  expect_output(hierarchical_adjustment(y, method="None", qualitycontrol = TRUE), "\\w*ASW Batch was\\w*")
})

test_that("BERT prints ASW Label if qualitycontrol is TRUE", {
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