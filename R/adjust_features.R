#' Check, which features contain enough numeric data to be adjusted (at least
#' 2 numeric values)
#'
#' This function will be called automatically be BERT n data from each batch
#' independently.
#'
#' @param data Matrix or dataframe in the format (samples, features). 
#' Additional column names are "Batch", "Cov_X" (were X may be any number),
#' "Label" and "Sample".
#' @return A logical with TRUE for adjustable features and FALSE for features
#' with too many missing values.
#' @export
get_adjustable_features <- function(data_batch) {
  # get the numeric data (boolean)
  available_data <- !is.na(data_batch)
  # add the booleans per column (per feature)--> counts of numeric values
  available_data <- colSums(available_data)
  
  return(available_data > 1)
  
}

#' Adjust two batches to each other.
#'
#' This function is called by the BERT algorithm and should not be called by
#' the user directly.
#'
#' @param data Matrix or dataframe in the format (samples, features). 
#' Additional column names are "Batch", "Cov_X" (were X may be any number),
#' "Label" and "Sample".
#' @param b1 The first batch to adjust.
#' @param b2 The second batch to adjust.
#' @param mod Dataframe with potential covariables to use. May be emty.
#' @param combatmode Integer, encoding the parameters to use for ComBat.
#' 1 (default)    par.prior = TRUE, mean.only = FALSE
#' 2              par.prior = TRUE, mean.only = TRUE
#' 3              par.prior = FALSE, mean.only = FALSE
#' 4              par.prior = FALSE, mean.only = TRUE
#' Will be ignored, if method=="limma".
#' @param method Adjustment method to use. Should either be "ComBat" or "limma".
#' "None" is also allowed for testing purposes and will yield no batch effect
#' correction.
#' @return A matrix/dataframe mirroring the shape of the input. The data will
#' be batch-effect adjusted by the specified method.
adjust_node <- function(data, b1, b2, mod, combatmode, method) {
  # debug output
  log_str <- paste("Adjusting data from batch ", b1, " and ", b2)
  logging::logdebug(log_str)
  
  # data from the two respective batches
  data_b_1 <- data[data$Batch == b1,]
  data_b_2 <- data[data$Batch == b2,]
  # make joint matrix (by adding the samples fow-wise)
  total_data <- rbind(data_b_1, data_b_2)
  
  # select corresponding mod --> may still be empty
  mod_b_1 <- data.frame(mod[data$Batch == b1,])
  mod_b_2 <- data.frame(mod[data$Batch == b2,])
  
  # get adjustable features
  av_b1 <- get_adjustable_features(data_b_1)
  av_b2 <- get_adjustable_features(data_b_2)
  
  # we can only adjust features, which are adjustable for both batches
  # if something is not adjustable, that means the feature is missing
  # for the entire batch, since format_DF will have been called prior
  total_adjustable <- av_b1 & av_b2
  
  # if there is only one adjustable protein/gene, this is not enough for
  # ComBat. BERT will return both batches unadjusted!
  if (sum(total_adjustable) == 1) {
    logging::logerror("Singular overlap between batches. Returning fully unadjusted data.")
    return(total_data)
  }
  
  # select only the adjustable features from both batches and combine them
  # into one matrix/dataframe
  adjustable_data_b_1 <- data_b_1[, total_adjustable]
  adjustable_data_b_2 <- data_b_2[, total_adjustable]
  total_adjustable_data <- rbind(adjustable_data_b_1, adjustable_data_b_2)
  
  # same for the two dataframes with the covariates --> may stil be empty
  total_mod <- rbind(mod_b_1, stats::setNames(mod_b_2, names(mod_b_1)))
  
  # list with the respective batches
  batch_list <- total_adjustable_data[["Batch"]]
  
  # drop columns that we want to ignore. That is, all allowed columns that don't
  # contain numeric values
  total_adjustable_data <- total_adjustable_data[,!names(total_adjustable_data) %in% c("Batch", "Sample", "Label")]
  
  # set combat mode
  if(combatmode ==1){
    parprior <- TRUE
    meanonly <- FALSE
  }else if(combatmode == 2){
    parprior <- TRUE
    meanonly <- TRUE
  }else if(combatmode == 3){
    parprior <- FALSE
    meanonly <- FALSE
  }else if(combatmode == 4){
    parprior <- FALSE
    meanonly <- TRUE
  }else{
    stop("Unknown ComBat mode.")
  }
  
  # the following if-statements use the respective adjustment method and
  # fill them into the total_data matrix
  # if the covariate-dataframe is empty
  if(dim(mod)[2]==0){
    # if we use ComBat adjustment
    if(method=="ComBat"){
      total_data[, names(total_adjustable_data)] <- suppressMessages(t(sva::ComBat(dat = t(total_adjustable_data), batch = batch_list, par.prior = parprior, mean.only = meanonly)))
    }else if(method=="limma"){
      # if we use limma
      total_data[, names(total_adjustable_data)] <- suppressMessages(t(limma::removeBatchEffect(x = t(total_adjustable_data), batch = batch_list)))
    }else if (method=="None"){
      total_data[, names(total_adjustable_data)] <- total_adjustable_data
    }else{
      stop("Unknown adjustment method.")
    }
  }else{
    if(method=="ComBat"){
      total_data[, names(total_adjustable_data)] <-
        suppressMessages(t(sva::ComBat(dat = t(total_adjustable_data), batch = batch_list, mod=total_mod, par.prior = parprior, mean.only = meanonly)))
    }else if(method=="limma"){
      total_data[, names(total_adjustable_data)] <- suppressMessages(t(limma::removeBatchEffect(x = t(total_adjustable_data), batch = batch_list, design=total_mod)))
    }else if (method=="None"){
      total_data[, names(total_adjustable_data)] <- total_adjustable_data
    }else{
      stop("Unknown adjustment method.")
    }
    
  }
  
  # return adjusted data
  return(total_data)
}