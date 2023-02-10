#' Strip column labelled Cov_1 from dataframe.
#'
#' @param dataset Dataframe in the shape (samples, features) with additional
#' column Cov_1
#' @return Dataset without column Cov_1.
#' @export
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
#' @param mvstmt Float (in (0,1)) indicating the fraction of missing values
#' per batch. 
#' @param classes Integer indicating the number of classes in the dataset.
#' @param housekeeping If NULL, no huosekeeping features will be simulatd.
#' Else, housepeeping indicates the fraction of of housekeeping features.
#' @return A dataframe containing the simulated data.
#' @export
generateDataset <- function(features, batches, samplesperbatch, mvstmt, classes, housekeeping = NULL){
  # genewise offset
  a <- rnorm(features, mean=0, sd=1)
  # condition-specific offset
  bix <- matrix(unlist(rnorm(features*classes, mean=0, sd=1)), nrow=features, ncol=classes)
  # the class values we may have
  potential_classes <- 1:classes
  # randomly select the class labels for each sample, with equal probability!
  classvector <- sample(potential_classes, batches*samplesperbatch, replace = TRUE)
  # evenly distribute samples over batches
  batchvector <- comprehenr::to_vec(for(i in 1:(batches*samplesperbatch)) (i %% batches)+1)
  # make matrix for the numeric expression values
  values <- matrix(0, ncol=features, nrow=batches*samplesperbatch)
  
  # fill with data, based on condition
  for(i in 1:(batches*samplesperbatch)){
    values[i,] <- a + bix[, classvector[i]]
  }
  
  # now add batch effects
  # add some normally distributed noise --> e.g. measurement error, epsilon in
  # L/S model
  noise <- matrix(unlist(rnorm(features*batches*samplesperbatch, mean=0, sd=0.1)), nrow=batches*samplesperbatch, ncol=features)
  
  # iterate over batches
  for(b in unique(batchvector)){
    # additive batch effect, normally distributed
    proteinshift <- rnorm(features, mean=0, sd=1)
    # multiplicative batch effect, inverse gamma
    proteinscale <- sqrt(invgamma::rinvgamma(features, shape=5, rate = 2))
    # for each sample in this batch
    for(index in which(batchvector==b)){
      # additive and multiplicative batch effect
      values[index, ] <- proteinshift + values[index, ] + proteinscale*noise[index, ]
    }
  }
  
  start_idx <- 1
  if(!is.null(housekeeping)){
    start_idx <- round(housekeeping*features, digits=0)
  }
  
  # introduce missing values for each batch --> TMT like
  for(b in unique(batchvector)){
    # randomly select features to be missing
    missingindices = sample(start_idx:features, round(mvstmt*features, digits = 0))
    # indices of samples from this batch
    batch_indices = which(batchvector==b)
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
#' @param mvstmt Float (in (0,1)) indicating the fraction of missing values
#' per batch. 
#' @param imbalcov Float indicating the probability for one of the classes to be
#' drawn as class label for each sample. The second class will have 
#' probability of 1-imbalcov 
#' @param housekeeping If NULL, no huosekeeping features will be simulatd.
#' Else, housepeeping indicates the fraction of of housekeeping features.
#' @return A dataframe containing the simulated data. Column Cov_1 will contain
#' the simulated, imbalanced labels.
#' @export
generateDataCovariables <- function(features, batches, samplesperbatch, mvstmt, imbalcov, housekeeping = NULL){
  # genewise offset
  a <- rnorm(features, mean=0, sd=1)
  # condition-specific offset
  bix <- matrix(unlist(rnorm(features*2, mean=0, sd=1)), nrow=features, ncol=2)
  # we only have two classes
  potential_classes <- 1:2
  # randomly select the class labels for each sample, with equal probability (here!)
  classvector <- sample(potential_classes, batches*samplesperbatch, replace = TRUE)
  # evenly distribute samples over batches
  batchvector <- comprehenr::to_vec(for(i in 1:(batches*samplesperbatch)) (i %% batches)+1)
  
  # make classes unbalanced
  for(b in unique(batchvector)){
    # for each batch, determine randomly, whether class 1 has probability imbalcov,
    # of class 2
    if(rnorm(1)>0){
      prob1 <- imbalcov
    }else{
      prob1 <- 1-imbalcov
    }
    prob2 <- 1-prob1
    # all samples from this batch
    indices = which(batchvector==b)
    # now overwrite the class labels --> this time with unbalanced datasets
    classvector[indices] = sample(c(1,2), size = length(indices), replace = TRUE, prob = c(prob1, prob2))
  }
  
  # make matrix for the numeric expression values
  values <- matrix(0, ncol=features, nrow=batches*samplesperbatch)
  
  # fill with data, based on condition
  for(i in 1:(batches*samplesperbatch)){
    values[i,] <- a + bix[, classvector[i]]
  }
  
  # now add batch effects
  # add some normally distributed noise --> e.g. measurement error, epsilon in
  # L/S model
  noise <- matrix(unlist(rnorm(features*batches*samplesperbatch, mean=0, sd=0.1)), nrow=batches*samplesperbatch, ncol=features)
  
  # iterate over batches
  for(b in unique(batchvector)){
    # additive batch effect, normally distributed
    proteinshift <- rnorm(features, mean=0, sd=1)
    # multiplicative batch effect, inverse gamma
    proteinscale <- sqrt(invgamma::rinvgamma(features, shape=5, rate = 2))
    # for each sample in this batch
    for(index in which(batchvector==b)){
      # additive and multiplicative batch effect
      values[index, ] <- proteinshift + values[index, ] + proteinscale*noise[index, ]
    }
  }
  
  start_idx <- 1
  if(!is.null(housekeeping)){
    start_idx <- round(housekeeping*features, digits=0)
  }
  # introduce missing values for each batch --> TMT like
  for(b in unique(batchvector)){
    # randomly select features to be missing
    missingindices = sample(start_idx:features, round(mvstmt*features, digits = 0))
    # indices of samples from this batch
    batch_indices = which(batchvector==b)
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


#' Compute the average silhouette width (ASW) for the dataset with respect
#' to both label and batch.
#'
#'Columns labelled Sample and Cov_1 will be ignored.
#'
#' @param dataset Dataframe in the shape (samples, features) with additional
#' columns Batch and Label.
#' @return List with fields "Label" and "Batch" for the ASW with regards to Label
#' and Batch respectively.
#' @export
compute_asw <- function(dataset){
  # labels as vector
  labels <- as.vector(dataset[["Label"]])
  # batches as vector
  batches <- as.vector(dataset[["Batch"]])
  # numeric values in dataset only
  num_values <- dataset[,!names(dataset) %in% c("Batch", "Sample", "Label", "Cov_1")]
  # compute distance matrix based on euclidean distances, ignoring NAs
  distancematrix <- dist(num_values)
  # compute silhouette object wrt. labels
  sil_labels <- cluster::silhouette(labels, dist=distancematrix)
  # extract ASW wrt. labels
  asw_label <- summary(sil_labels)["avg.width"]$avg.width
  # compute silhouette object wrt. batches
  sil_batches <- cluster::silhouette(batches, dist=distancematrix)
  # extract ASW wrt. batches
  asw_batches <- summary(sil_batches)["avg.width"]$avg.width
  # create list object
  ret <- list("Label" = asw_label, "Batch" = asw_batches)
  return(ret)
}

#' Count the number of numeric features in this dataset. Columns labelled 
#' "Batch", "Sample" or "Label" will be ignored.
#'
#' @param dataset Dataframe in the shape (samples, features) with optional
#' columns "Batch", "Sample" or "Label".
#' @return Integer indicating the number of numeric values
#' @export
count_existing <- function(dataset){
  # select only numeric columns
  num_values <- dataset[,!names(dataset) %in% c("Batch", "Sample", "Label", "Cov_1")]
  # sum up non-missing values
  return(sum(!is.na(num_values)))
}