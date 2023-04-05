#' Ordinal encoding of a vector.
#' 
#'This function is usually called by BERT during formatting of the input.
#'The idea is, that Label, Batch and Covariables should only be integers
#' @param column The categorical vector
#' @return The encoded vector
#' @export
ordinal_encode <- function(column){
  temp <- as.integer(factor(column, levels=unique(column)))
  return(temp)
}

#' Replaces missing values (NaN) by NA, this appears to be faster
#' 
#' @param data The data as dataframe
#' @return The data with the replaced MVs
#' @export
replace_missing <- function(data){
  data[sapply(data, is.nan)] <- NA
  data[is.null(data)] <- NA
  return(data)
}


#' Format the data as expected by BERT.
#'
#'This function is called automatically by BERT. It removes empty columns
#'and removes a (usually very small) number of numeric values, if features are
#'unadjustable for lack of data.
#'
#' @param data Matrix or dataframe in the format (samples, features). 
#' Additional column names are "Batch", "Cov_X" (were X may be any number),
#' "Label" and "Sample".
#' @return The formatted matrix.
format_DF <- function(data){
  logging::loginfo("Formatting Data.")
  
  if(is.matrix(data)){
    logging::loginfo("Typecasting input to dataframe.")
    data <- data.frame(data)
  }
  
  logging::loginfo("Replacing NaNs with NAs.")
  data <- replace_missing(data)
  
  # get names of potential covariables
  cov_names <- names(data)[grepl( "Cov" , names( data  ) )]
  cat_names <- names(data)[names(data) %in% c("Label", "Batch")]
  all_names <- c(cov_names, cat_names)
  
  if(length(all_names)==1){
    if(!is.character(data[1, all_names])){
      all_names <- character(0)
    }
  }else{
    dtypes <- sapply(data[, all_names], typeof)
    all_names <- all_names[dtypes=="character"]
  }
  
  if (length(all_names>0)){
    logging::logwarn(paste("Identified", length(all_names),
                              "categorical variables among batch, label and all covariates. Note that BERT requires integer values there. Will apply ordinal encoding."))
    
    for(n in all_names){
      data[, n] <- ordinal_encode(data[[n]])
    }
  }
  
  
  logging::loginfo("Removing potential empty rows and columns")
  `%>%` <- janitor::`%>%`
  data <- data %>% janitor::remove_empty(c("rows", "cols"))
  
  # count number of missing values
  inital_mvs <- sum(is.na(data))
  
  logging::loginfo(paste("Found ", inital_mvs, " missing values."))
  
  # all unique batch levels
  unique_batches <- unique(data[["Batch"]])
  
  # select covariates
  mod <- data.frame(data [ , grepl( "Cov" , names( data  ) ) ])
  
  if(dim(mod)[2]!=0){
    logging::loginfo("BERT requires at least 2 numeric values per batch/covariate level. This may reduce the number of adjustable features considerably, depending on the quantification technique.")
  }
  
  # iterate over batches and remove numeric values, if a feature (e.g. protein)
  # does not contain at least 2 numeric values
  for(b in unique_batches){
    # data from batch b
    data_batch <- data[data["Batch"] == b,]
    mod_batch <- mod[data["Batch"] == b,]
    # logical with the features that can be adjusted (that is, contain more
    # than 2 numeric values in this batch/covariate level)
    if(dim(mod)[2]==0){
      adjustable_batch <- get_adjustable_features(data_batch)
    }else{
      adjustable_batch <- get_adjustable_features_with_mod(data_batch, data.frame(mod_batch))
    }
    # set features from this batch to missing, where adjustable_batch is FALSE
    data[data["Batch"] == b, !adjustable_batch] <- NA
    
  }
  # count missing values
  final_mvs <- sum(is.na(data))
  
  logging::loginfo(paste("Introduced ", final_mvs-inital_mvs, " missing values due to singular proteins at batch/covariate level."))
  logging::loginfo("Done")
  
  return(data)
}