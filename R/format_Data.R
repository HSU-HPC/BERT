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
  
  logging::loginfo("Removing potential empty rows and columns")
  `%>%` <- janitor::`%>%`
  data <- data %>% janitor::remove_empty(c("rows", "cols"))
  
  # count number of missing values
  inital_mvs <- sum(is.na(data))
  
  logging::loginfo(paste("Found ", inital_mvs, " missing values."))
  
  # all unique batch levels
  unique_batches <- unique(data[["Batch"]])
  
  # iterate over batches and remove numeric values, if a feature (e.g. protein)
  # does not contain at least 2 numeric values
  for(b in unique_batches){
    # data from batch b
    data_batch <- data[data["Batch"] == b,]
    # logical with the features that can be adjusted (that is, contain more
    # than 2 numeric values in this batch)
    adjustable_batch <- get_adjustable_features(data_batch)
    # set features from this batch to missing, where adjustable_batch is FALSE
    data[data["Batch"] == b, !adjustable_batch] <- NA
    
  }
  # count missing values
  final_mvs <- sum(is.na(data))
  
  logging::loginfo(paste("Introduced ", final_mvs-inital_mvs, " missing values due to singular proteins in batches."))
  logging::loginfo("Done")
  
  return(data)
}