test_that("correctly counts numeric values", {
  # generate dataset, 3 samples, 5 features
  y <- matrix(rnorm(3*5),3,5)
  y[1,1] <- NA
  y[1:3, 2] <- NA
  y[1:2, 3] <- NA
  # BERT typecasts to dataframe
  y <- data.frame(y)
  expect_equal(count_existing(y), 9)
})

test_that("validate input generate_dataset works", {
    # should work
    expect_error(validate_input_generate_dataset(features=1, batches=3,
                                                 samplesperbatch=4, mvstmt=0.1,
                                                 classes=2, housekeeping=NULL,
                                                 deterministic=FALSE), NA)
    expect_error(validate_input_generate_dataset(features=1, batches=3,
                                                 samplesperbatch=4, mvstmt=0.1,
                                                 classes=2, housekeeping=0.1,
                                                 deterministic=FALSE), NA)
    expect_error(validate_input_generate_dataset(features=1, batches=3,
                                                 samplesperbatch=2, mvstmt=0.0,
                                                 classes=2, housekeeping=NULL,
                                                 deterministic=FALSE), NA)
    # should crash
    expect_error(validate_input_generate_dataset(features=-1, batches=3,
                                                 samplesperbatch=4, mvstmt=0.1,
                                                 classes=2, housekeeping=NULL,
                                                 deterministic=FALSE))
    expect_error(validate_input_generate_dataset(features="test", batches=3,
                                                 samplesperbatch=4, mvstmt=0.1,
                                                 classes=2, housekeeping=NULL,
                                                 deterministic=FALSE))
    expect_error(validate_input_generate_dataset(features=1, batches=1,
                                                 samplesperbatch=4, mvstmt=0.1,
                                                 classes=2, housekeeping=NULL,
                                                 deterministic=FALSE))
    expect_error(validate_input_generate_dataset(features=1, batches="hi",
                                                 samplesperbatch=4, mvstmt=0.1,
                                                 classes=2, housekeeping=NULL,
                                                 deterministic=FALSE))
    expect_error(validate_input_generate_dataset(features=1, batches=3,
                                                 samplesperbatch=1, mvstmt=0.1,
                                                 classes=2, housekeeping=NULL,
                                                 deterministic=FALSE))
    expect_error(validate_input_generate_dataset(features=1, batches=3,
                                                 samplesperbatch="", mvstmt=0.1,
                                                 classes=2, housekeeping=NULL,
                                                 deterministic=FALSE))
    expect_error(validate_input_generate_dataset(features=1, batches=3,
                                                 samplesperbatch=1, mvstmt=-0.1,
                                                 classes=2, housekeeping=NULL,
                                                 deterministic=FALSE))
    expect_error(validate_input_generate_dataset(features=1, batches=3,
                                                 samplesperbatch=1, mvstmt="",
                                                 classes=2, housekeeping=NULL,
                                                 deterministic=FALSE))
    expect_error(validate_input_generate_dataset(features=1, batches=3,
                                                 samplesperbatch=4, mvstmt=0.1,
                                                 classes=0, housekeeping=NULL,
                                                 deterministic=FALSE))
    expect_error(validate_input_generate_dataset(features=1, batches=3,
                                                 samplesperbatch=4, mvstmt=0.1,
                                                 classes=FALSE, housekeeping=NULL,
                                                 deterministic=FALSE))
    expect_error(validate_input_generate_dataset(features=1, batches=3,
                                                 samplesperbatch=4, mvstmt=0.1,
                                                 classes=FALSE, housekeeping=0.1,
                                                 deterministic=FALSE))
    expect_error(validate_input_generate_dataset(features=1, batches=3,
                                                 samplesperbatch=4, mvstmt=0.3,
                                                 classes=FALSE, housekeeping=0.9,
                                                 deterministic=FALSE))
    expect_error(validate_input_generate_dataset(features=1, batches=3,
                                                 samplesperbatch=4, mvstmt=0.1,
                                                 classes=FALSE, housekeeping=NULL,
                                                 deterministic=""))
    # generate_dataset should call the check
    expect_error(generate_dataset(features=1, batches=3,
                                  samplesperbatch=4, mvstmt=0.1,
                                  classes=FALSE, housekeeping=NULL,
                                  deterministic=""))
    # same for generate_data_covariables
    expect_error(generate_data_covariables(features=1, batches=3,
                                           samplesperbatch=4, mvstmt=0.1,
                                           imbalcov = 0.1, housekeeping=""))
    expect_error(generate_data_covariables(features=1, batches=3,
                                           samplesperbatch=4, mvstmt=0.1,
                                           imbalcov = 1.1, housekeeping=NULL))
    
    
})

test_that("deterministic generation of datasets works", {
  ds <- generate_dataset(10, 4, 8, 0.1, 2, deterministic = TRUE)
  classes <- ds[["Label"]]
  expect_equal(max(classes), 2)
  expect_equal(min(classes), 1)
  expect_equal(sum(classes==1), 16) # class 2 also occurs 16 times then
  expect_equal(sum(classes[ds$Batch==1]==1), 4)
})



test_that("Ignores Batch, Label, Sample and covariable columns", {
  # generate dataset, 3 samples, 5 features
  y <- matrix(rnorm(3*5),3,5)
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

test_that("count_existing only takes dataframes as input", {
    # generate dataset, 3 samples, 5 features
    y <- matrix(rnorm(3*5),3,5)
    y[1,1] <- NA
    y[1:3, 2] <- NA
    y[1:2, 3] <- NA
    # this should crash
    expect_error(count_existing(y))
    y <- data.frame(y)
    y["Batch"] = 1
    y["Label"] = 2
    y["Cov_1"] = 4
    y["Cov_2"] = 20
    # and this should not
    expect_equal(count_existing(y), 9)
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

test_that("generate_truncated_dataset works", {
    # does not crash (bad style...)
    data <- generate_truncated_dataset(100, 3, 10, 0.1,0.1, 2)
    # should error if mvs truncated is string or in wrong interval
    expect_error(generate_truncated_dataset(100, 3, 10, 0.1,"test", 2))
    expect_error(generate_truncated_dataset(100, 3, 10, 0.1,-1, 2))
    expect_error(generate_truncated_dataset(100, 3, 10, 0.1,1.2, 2))
    # should error if absolute amount of missing values would be >100%
    expect_error(generate_truncated_dataset(100, 3, 10, 0.8,0.8, 2))
})


