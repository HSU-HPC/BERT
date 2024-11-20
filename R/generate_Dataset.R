#' Validate the user input to the function generate_dataset.
#' Raises an error if and only if the input is malformatted.
#' 
#' @param features Integer indicating the number of features
#' (e.g. genes/proteins) in the dataset.
#' @param batches Integer indicating the number of batches in the
#'  dataset.
#' @param samplesperbatch Integer indicating the number of of samples
#' per batch.
#' @param mvstmt Float (in [0,1)) indicating the fraction of missing values
#' per batch. 
#' @param classes Integer indicating the number of classes in the dataset.
#' @param housekeeping If NULL, no huosekeeping features will be simulatd.
#' Else, housepeeping indicates the fraction of of housekeeping features.
#' @param deterministic Whether to assigns the classes deterministically,
#' instead of random sampling
#' @return None
validate_input_generate_dataset <- function(features, 
                                            batches, 
                                            samplesperbatch, 
                                            mvstmt, 
                                            classes, 
                                            housekeeping, 
                                            deterministic) {
    if(!is.numeric(features) || !features%%1==0 || features<=0){
        error_str <- paste("Please provide positive integer arguments to",
                           "features in function generate_dataset")
        stop(error_str)
    }
    if(!is.numeric(batches) || !batches%%1==0 || batches<=1){
        error_str <- paste("Please provide integer (>=2) arguments to",
                           "batches in function generate_dataset")
        stop(error_str)
    }
    if(!is.numeric(samplesperbatch) || !samplesperbatch%%1==0 ||
       samplesperbatch<=1){
        error_str <- paste("Please provide integer (>=2) arguments to",
                           "samplesperbatch in function generate_dataset")
        stop(error_str)
    }
    if(!is.numeric(mvstmt) || !(0<=mvstmt && mvstmt<1)){
        error_str <- paste("Parameter mvstmt in function generate_dataset",
                           "must be in [0,1)")
        stop(error_str)
    }
    if(!is.numeric(classes) || !classes%%1==0 || classes<=0){
        error_str <- paste("Please provide positive integer arguments to",
                           "classes in function generate_dataset")
        stop(error_str)
    }
    if(!is.null(housekeeping)){
        if(!is.numeric(housekeeping) || !(0<housekeeping && housekeeping<1)){
            error_str <- paste("Parameter housekeeping in function",
                               "generate_dataset must be in NULL or",
                               "(0,1)")
            stop(error_str)
        }
        if(mvstmt + housekeeping>1){
            error_str <- paste("Sum of parameters mvstmt and housekeeping in",
                               "generate_dataset must be <=1.")
            stop(error_str)
        }
    }
    if(!is.logical(deterministic)){
        error_str <- paste("Parameter deterministic in function",
                           "generate_dataset must be TRUE/FALSE.")
        stop(error_str)
    }
}


#' Strip column labelled Cov_1 from dataframe.
#'
#' @param dataset Dataframe in the shape (samples, features) with additional
#' column Cov_1
#' @return Dataset without column Cov_1.
strip_Covariable <- function(dataset){
    ds <- dataset[, !names(dataset) %in% c("Cov_1")] 
    return(ds)
}

#' Generate dataset with batch-effects and biological labels using
#' a simple LS model
#'
#'The data will be already in the correct format for BERT.
#'
#' @param features Integer indicating the number of features
#' (e.g. genes/proteins) in the dataset.
#' @param batches Integer indicating the number of batches in the
#'  dataset.
#' @param samplesperbatch Integer indicating the number of of samples
#' per batch.
#' @param mvstmt Float (in [0,1)) indicating the fraction of missing values
#' per batch. 
#' @param classes Integer indicating the number of classes in the dataset.
#' @param housekeeping If NULL, no huosekeeping features will be simulatd.
#' Else, housepeeping indicates the fraction of of housekeeping features.
#' @param deterministic Whether to assigns the classes deterministically,
#' instead of random sampling
#' @return A dataframe containing the simulated data.
#' @examples
#' # generate dataset with 1000 features, 5 batches, 10 samples per batch and
#' # two genotypes
#' data = generate_dataset(1000,5,10, 0.1, 2)
#' @export
generate_dataset <- function(
        features, 
        batches, 
        samplesperbatch, 
        mvstmt, 
        classes, 
        housekeeping = NULL, 
        deterministic = FALSE
        ){
    
    # validate input
    validate_input_generate_dataset(features, batches, samplesperbatch, mvstmt,
                                    classes, housekeeping, deterministic)
    
    # genewise offset
    a <- stats::rnorm(features, mean=0, sd=1)
    # condition-specific offset
    bix <- matrix(
        unlist(stats::rnorm(features*classes, mean=0, sd=1)), 
        nrow=features, 
        ncol=classes)
    # evenly distribute samples over batches
    batchvector <- comprehenr::to_vec(for(i in seq_len(
        batches*samplesperbatch)) (i %% batches)+1)
    # the class values we may have
    potential_classes <- seq_len(classes)
    if(deterministic){
        classvector <- rep(0, batches*samplesperbatch)
        for(b in unique(batchvector)){
            classvector[batchvector==b] <- seq_len(samplesperbatch)%%classes
        }
        classvector <- classvector+1
    }else{
        # randomly select the class labels for each sample, with equal 
        # probability!
        classvector <- sample(
            potential_classes, 
            batches*samplesperbatch, 
            replace = TRUE)
    }
    # make matrix for the numeric expression values
    values <- matrix(0, ncol=features, nrow=batches*samplesperbatch)
    
    # fill with data, based on condition
    for(i in seq_len(batches*samplesperbatch)){
        values[i,] <- a + bix[, classvector[i]]
    }
    
    # now add batch effects
    # add some normally distributed noise --> e.g. measurement error, epsilon in
    # L/S model
    noise <- matrix(
        unlist(stats::rnorm(features*batches*samplesperbatch, mean=0, sd=0.1)), 
        nrow=batches*samplesperbatch, 
        ncol=features)
    
    # iterate over batches
    for(b in unique(batchvector)){
        # additive batch effect, normally distributed
        proteinshift <- stats::rnorm(features, mean=0, sd=1)
        # multiplicative batch effect, inverse gamma
        proteinscale <- sqrt(invgamma::rinvgamma(features, shape=5, rate = 2))
        # for each sample in this batch
        for(index in which(batchvector==b)){
            # additive and multiplicative batch effect
            values[index, ] <- proteinshift + 
                values[index, ] + 
                proteinscale*noise[index, ]
        }
    }
    
    start_idx <- 1
    if(!is.null(housekeeping)){
        start_idx <- round(housekeeping*features, digits=0)
    }
    
    # introduce missing values for each batch --> TMT like
    for(b in unique(batchvector)){
        # randomly select features to be missing
        missingindices <- sample(
            start_idx:features, 
            round(mvstmt*features, digits = 0))
        # indices of samples from this batch
        batch_indices <- which(batchvector==b)
        # set values to NA
        values[batch_indices, missingindices] <- NA
    }
    
    # make data frame
    finaldf <- data.frame(values)
    # add column with batch
    finaldf["Batch"] <- batchvector
    # add column with label
    finaldf["Label"] <- classvector
    
    return(finaldf)
}

#' Generate dataset with batch-effects and 2 classes with a specified imbalance.
#'
#'The data will be already in the correct format for BERT.
#'
#' @param features Integer indicating the number of features
#' (e.g. genes/proteins) in the dataset.
#' @param batches Integer indicating the number of batches in the
#'  dataset.
#' @param samplesperbatch Integer indicating the number of of samples
#' per batch.
#' @param mvstmt Float (in [0,1)) indicating the fraction of missing values
#' per batch. 
#' @param imbalcov Float indicating the probability for one of the classes to be
#' drawn as class label for each sample. The second class will have 
#' probability of 1-imbalcov 
#' @param housekeeping If NULL, no huosekeeping features will be simulatd.
#' Else, housepeeping indicates the fraction of of housekeeping features.
#' @return A dataframe containing the simulated data. Column Cov_1 will contain
#' the simulated, imbalanced labels.
#' @examples
#' # generate dataset with 1000 features, 5 batches, 10 samples per batch and
#' # two genotypes. The class ratio will either be 7:3 or 3:7 per batch.
#' data = generate_data_covariables(1000,5,10, 0.1, 0.3)
#' @export
generate_data_covariables <- function(
        features, 
        batches, 
        samplesperbatch, 
        mvstmt, 
        imbalcov, 
        housekeeping = NULL){
    # validity check for input (except for imbalcov)
    validate_input_generate_dataset(features, batches, samplesperbatch, mvstmt,
                                    2, housekeeping, FALSE)
    if(!is.numeric(imbalcov) || !(0<imbalcov && imbalcov<1)){
        error_str <- paste("Parameter imbalcov should be in (0,1)",
                           "in function generate_data_covariables")
        stop(error_str)
    }
    
    # genewise offset
    a <- stats::rnorm(features, mean=0, sd=1)
    # condition-specific offset
    bix <- matrix(
        unlist(stats::rnorm(features*2, mean=0, sd=1)), 
        nrow=features, 
        ncol=2)
    # we only have two classes
    potential_classes <- seq_len(2)
    # randomly select the class labels for each sample, with equal 
    # probability (here!)
    classvector <- sample(
        potential_classes, 
        batches*samplesperbatch, 
        replace = TRUE)
    # evenly distribute samples over batches
    batchvector <- comprehenr::to_vec(for(i in seq_len(
        batches*samplesperbatch)) (i %% batches)+1)
    
    # make classes unbalanced
    for(b in unique(batchvector)){
        # for each batch, determine randomly, whether class 1 has probability 
        # imbalcov, of class 2
        if(stats::rbinom(1, 1, 0.5)>0){
            prob1 <- imbalcov
        }else{
            prob1 <- 1-imbalcov
        }
        prob2 <- 1-prob1
        # all samples from this batch
        indices <- which(batchvector==b)
        # now overwrite the class labels --> this time with unbalanced datasets
        classvector[indices] <- sample(
            c(1,2), 
            size = length(indices), 
            replace = TRUE, 
            prob = c(prob1, prob2))
    }
    
    # make matrix for the numeric expression values
    values <- matrix(0, ncol=features, nrow=batches*samplesperbatch)
    
    # fill with data, based on condition
    for(i in seq_len(batches*samplesperbatch)){
        values[i,] <- a + bix[, classvector[i]]
    }
    
    # now add batch effects
    # add some normally distributed noise --> e.g. measurement error, epsilon in
    # L/S model
    noise <- matrix(
        unlist(stats::rnorm(features*batches*samplesperbatch, mean=0, sd=0.1)), 
        nrow=batches*samplesperbatch, 
        ncol=features)
    
    # iterate over batches
    for(b in unique(batchvector)){
        # additive batch effect, normally distributed
        proteinshift <- stats::rnorm(features, mean=0, sd=1)
        # multiplicative batch effect, inverse gamma
        proteinscale <- sqrt(invgamma::rinvgamma(features, shape=5, rate = 2))
        # for each sample in this batch
        for(index in which(batchvector==b)){
            # additive and multiplicative batch effect
            values[index, ] <- proteinshift + values[index, ] + 
                proteinscale*noise[index, ]
        }
    }
    
    start_idx <- 1
    if(!is.null(housekeeping)){
        start_idx <- round(housekeeping*features, digits=0)
    }
    # introduce missing values for each batch --> TMT like
    for(b in unique(batchvector)){
        # randomly select features to be missing
        missingindices <- sample(
            start_idx:features, 
            round(mvstmt*features, digits = 0))
        # indices of samples from this batch
        batch_indices <- which(batchvector==b)
        # set values to NA
        values[batch_indices, missingindices] <- NA
    }
    
    # make dataframes and append columns for batch, label and covariate
    # covariate is equivalent to class label
    finaldf <- data.frame(values)
    finaldf["Batch"] <- batchvector
    finaldf["Label"] <- classvector
    finaldf["Cov_1"] <- classvector
    
    return(finaldf)
}

#' Generate dataset using simple L/S model as in ComBat
#' but introducing missing values in both TMT-like manner
#' and truncating a fixed amount of numerical values (to simulate
#' MNAR missing values from e.g. DDA acquisition in high-complexity
#' spectra)
#'
#'The data will be already in the correct format for BERT.
#'
#' @param features Integer indicating the number of features
#' (e.g. genes/proteins) in the dataset.
#' @param batches Integer indicating the number of batches in the
#'  dataset.
#' @param samplesperbatch Integer indicating the number of of samples
#' per batch.
#' @param mvstmt Float (in [0,1)) indicating the fraction of missing values
#' per batch. 
#' @param mvstruncated Float (in [0,1)) indicating the fraction of missing
#' values from truncation per batch. Sum of mvstruncated and mvstmt must
#' not be >1. Truncation happens by setting the corresponding fraction
#' of features with lowest mean expression values to NA.
#' @param classes Integer indicating the number of classes in the dataset.
#' @return A dataframe containing the simulated data.
#' @examples
#' # generate dataset with 1000 features, 5 batches, 10 samples per batch,
#' 10% tmt like missing values, 10% truncation NAs and
#' # two genotypes
#' data = generate_truncated_dataset(1000,5,10, 0.1, 0.1, 2)
#' @export
generate_truncated_dataset <- function(
        features, 
        batches, 
        samplesperbatch, 
        mvstmt,
        mvstruncated,
        classes){
    # generate complete data
    complete_data <- generate_dataset(features,
                                      batches,
                                      samplesperbatch,
                                      0,
                                      classes)
    # validate user input
    if(!is.numeric(mvstmt) || !(0<=mvstmt && mvstmt<1)){
        error_str <- paste("Parameter mvstmt in function",
                           " generate_truncated_dataset must be in [0,1)")
        stop(error_str)
    }
    if(!is.numeric(mvstruncated) || !(0<=mvstruncated && mvstruncated<1)){
        error_str <- paste("Parameter mvstruncated in function",
                           "generate_truncated_dataset must be in [0,1)")
        stop(error_str)
    }
    if(mvstmt + mvstruncated>1){
        error_str <- paste("Parameters mvstmt and mvstruncated",
                           " must have sum in [0,1]")
        stop(error_str)
    }
    # 1. get unique batches
    batches <- unique(complete_data$Batch)
    # 2. iterate over all batches
    for(b in batches){
        # select batch
        batch <- complete_data[complete_data$Batch==b, ]
        # introduce truncated data
        mean_expression <- colMeans(batch[,!names(batch) %in% c("Batch",
                                                                "Sample",
                                                                "Label",
                                                                "Reference")])
        lowest_indices <- order(mean_expression)[1:round(
            mvstruncated*features)]
        complete_data[complete_data$Batch==b, lowest_indices] <- NA
        # now introduce TMT-like missing values
        all_indices <- 1:features
        tmt_missing <- sample(all_indices[-lowest_indices],
                              round(mvstmt*features), replace = FALSE)
        complete_data[complete_data$Batch==b, tmt_missing] <- NA
    }
    # and return
    return(complete_data)
}

#' Count the number of numeric features in this dataset. Columns labeled 
#' "Batch", "Sample" or "Label" will be ignored.
#'
#' @param dataset Dataframe in the shape (samples, features) with optional
#' columns "Batch", "Sample" or "Label".
#' @return Integer indicating the number of numeric values
#' @examples
#' # generate dataset with 1000 features, 5 batches, 10 samples per batch and
#' # two genotypes
#' data = generate_dataset(1000,5,10, 0.1, 2)
#' count_existing(data)
#' @export
count_existing <- function(dataset){
    
    if(!is.data.frame(dataset)){
        error_str <- paste("Parameter dataset in function count_existing",
                           "must be of type dataframe.")
        stop(error_str)
    }
    
    dataset_nocov <- dataset [ , !grepl( "Cov" , names( dataset  ) ) ]
    # select only numeric columns
    num_values <- dataset_nocov[,!names(dataset_nocov) %in% c(
        "Batch", 
        "Sample", 
        "Label", 
        "Reference")]
    # sum up non-missing values
    return(sum(!is.na(num_values)))
}