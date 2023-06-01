#' Adjust a hierarchy level in parallel
#'
#' This function uses ComBat or limma to adjust an entire hierarchy level.
#'
#' @param data Matrix or dataframe in the format (samples, features). 
#' Additional column names are \"Batch\", \"Cov_X\" (were X may be any number),
#' \"Label\" and \"Sample\".
#' @param mod Dataframe with potential covariables to use. May be emty.
#' @param combatmode Integer, encoding the parameters to use for ComBat.
#' 1 (default)    par.prior = TRUE, mean.only = FALSE
#' 2              par.prior = TRUE, mean.only = TRUE
#' 3              par.prior = FALSE, mean.only = FALSE
#' 4              par.prior = FALSE, mean.only = TRUE
#' Will be ignored, if method=="limma".
#' @param method Adjustment method to use. Should either be \"ComBat\" or 
#' \"limma\".
#' @return A matrix/dataframe mirroring the shape of the input. The data will
#' be batch-effect adjusted by BERT.
adjustment_step_parallel <- function(data, mod, combatmode, method) {
    
    # get the unique batch values at this hierarchy level
    unique_batches <- unique(data$Batch)
    
    # indices for every second batch
    indices <- seq(1, length(unique_batches), by = 2)
    
    # define dopar in this namespace
    `%dopar%` <- foreach::`%dopar%`
    i <- NULL
    # adjust the data at this hierarchy level
    adjusted_data <-foreach::foreach(
        i = iterators::iter(indices),
        .combine = rbind, 
        .export = c("get_adjustable_features", "adjust_node")) %dopar% {
        # matrix / dataframe containing the adjusted data
        tempMatrix <- NULL
        if (i == length(unique_batches)) {
            # odd number of batches and this is the last one
            tempMatrix <- data[data$Batch == unique_batches[i],]
        } else{
            # can adjust a pair of batches (this and the last one)
            tempMatrix <- adjust_node(
                data, 
                unique_batches[i], 
                unique_batches[i + 1], 
                mod, 
                combatmode, 
                method)
        }
        # override batch description for the respectively adjusted batches
        tempMatrix["Batch"] <- unique_batches[i]
        
        # the adjusted data OR the single, unadjusted batch, if batch number is
        # odd
        tempMatrix
    }
    # adjusted data from this hierarchy level
    return(adjusted_data)
    
}


#' Adjust a hierarchy level sequentially.
#'
#' This function uses ComBat or limma to adjust an entire hierarchy level.
#'
#' @param data Matrix or dataframe in the format (samples, features). 
#' Additional column names are \"Batch\", \"Cov_X\" (were X may be any number),
#' \"Label\" and \"Sample\".
#' @param mod Dataframe with potential covariables to use. May be emty.
#' @param combatmode Integer, encoding the parameters to use for ComBat.
#' 1 (default)    par.prior = TRUE, mean.only = FALSE
#' 2              par.prior = TRUE, mean.only = TRUE
#' 3              par.prior = FALSE, mean.only = FALSE
#' 4              par.prior = FALSE, mean.only = TRUE
#' Will be ignored, if method=="limma".
#' @param method Adjustment method to use. Should either be \"ComBat\" or
#' \"limma\".
#' @return A matrix/dataframe mirroring the shape of the input. The data will
#' be batch-effect adjusted by BERT.
adjustment_step <- function(data, mod, combatmode, method) {
    
    # get the unique batch values at this hierarchy level
    unique_batches <- unique(data$Batch)
    
    # indices for every second batch
    indices <- seq(1, length(unique_batches), by = 2)
    
    # define do in this namespace
    `%do%` <- foreach::`%do%`
    i <- NULL
    # adjust the data at this hierarchy level
    adjusted_data <- foreach::foreach(
        i = iterators::iter(indices), .combine = rbind) %do% {
        # matrix / dataframe containing the adjusted data
        tempMatrix <- NULL
        if (i == length(unique_batches)) {
            # odd number of batches and this is the last one
            tempMatrix <- data[data$Batch == unique_batches[i],]
        } else{
            # can adjust a pair of batches (this and the last one)
            tempMatrix <- adjust_node(
                data, 
                unique_batches[i], 
                unique_batches[i + 1], 
                mod, 
                combatmode, 
                method)
        }
        # override batch description for the respectively adjusted batches
        tempMatrix["Batch"] <- unique_batches[i]
        # the adjusted data OR the single, unadjusted batch, if batch number is
        # odd
        tempMatrix
    }
    
    # adjusted data from this hierarchy level
    return(adjusted_data)
    
}