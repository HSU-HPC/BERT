test_that("ordinal encoding works", {
    col = c("A", "A", "B", "B")
    enc <- ordinal_encode(col)
    expected <- c(1,1,2,2)
    expect_true(all.equal(expected, enc))
})

test_that("ASW is identical for string and integer representations of Batch and Label", {
    # generate dataset, 3 samples, 5 features
    y <- matrix(rnorm(20*20),20,20)
    # BERT typecasts to dataframe
    y <- data.frame(y)
    y["Batch"] = sample(c(1,2), size=20, replace = TRUE)
    y["Label"] = sample(c(1,2), size=20, replace = TRUE)
    
    asw <- compute_asw(y)
    integer_asw_batch = asw$Batch
    integer_asw_label = asw$Label
    
    y$Batch = sapply(y$Batch, as.character)
    y$Label = sapply(y$Label, as.character)
    asw <- compute_asw(y)
    str_asw_batch = asw$Batch
    str_asw_label = asw$Label
    expect_equal(integer_asw_label, str_asw_label)
    expect_equal(integer_asw_batch, str_asw_batch)
})

test_that("ASW can also be computed for Label alone", {
    # generate dataset, 3 samples, 5 features
    y <- matrix(rnorm(20*20),20,20)
    # BERT typecasts to dataframe
    y <- data.frame(y)
    y["Label"] = sample(c(1,2), size=20, replace = TRUE)
    
    asw <- compute_asw(y)
    
    expect_true(!is.na(asw$Label)&&is.na(asw$Batch))
})

test_that("Computation of ASW works the same for batch and label", {
    # generate dataset, 3 samples, 5 features
    y <- matrix(rnorm(20*20),20,20)
    # BERT typecasts to dataframe
    y <- data.frame(y)
    y["Batch"] = sample(c(1,2), size=20, replace = TRUE)
    y["Label"] = y[["Batch"]]
    
    asw <- compute_asw(y)
    
    expect_equal(asw$Label, asw$Label)
})


test_that("ASW can also be computed for Batch alone", {
    # generate dataset, 3 samples, 5 features
    y <- matrix(rnorm(20*20),20,20)
    # BERT typecasts to dataframe
    y <- data.frame(y)
    y["Batch"] = sample(c(1,2), size=20, replace = TRUE)
    
    asw <- compute_asw(y)
    
    expect_true(is.na(asw$Label)&&!is.na(asw$Batch))
})


test_that("ASW ignores Sample and covariable columns", {
    # generate dataset, 3 samples, 5 features
    y <- matrix(rnorm(3*5),3,5)
    y[1,1] <- NA
    y[1:3, 2] <- NA
    y[1:2, 3] <- NA
    # BERT typecasts to dataframe
    y <- data.frame(y)
    y["Batch"] = c(1,1,2)
    y["Label"] = c(1,1,2)
    y2 <- data.frame(y)
    y2["Cov_1"] = 4
    y2["Cov_2"] = 20
    
    no_cat <- compute_asw(y)
    with_cat <- compute_asw(y2)
    
    expect_equal(no_cat$Label, with_cat$Label)
    expect_equal(no_cat$Batch, with_cat$Batch)
})

test_that("compute_asw takes only dataframes as input", {
    # generate dataset, 3 samples, 5 features
    y <- matrix(rnorm(3*5),3,5)
    y[1,1] <- NA
    y[1:3, 2] <- NA
    y[1:2, 3] <- NA
    # this should crash
    expect_error(compute_asw(y))
    y <- data.frame(y)
    y["Batch"] = c(1,1,2)
    y["Label"] = c(1,1,2)
    # this should work
    expect_error(compute_asw(y), NA)
})