#' Ordinal encoding of a vector.
#' 
#'This function is usually called by BERT during formatting of the input.
#'The idea is, that Label, Batch and Covariables should only be integers
#' @param column The categorical vector
#' @return The encoded vector
ordinal_encode <- function(column){
  temp <- as.integer(factor(column, levels=unique(column)))
  return(temp)
}

#' Verify that the Reference column of the data contains only zeros and ones
#' (if it is present at all)
#' @param batch the dataframe for this batch (samples in rows, samples in columns)
#' @return either TRUE (everything correct) or FALSE (something is not correct)
verify_references <- function(batch){
  if ("Reference" %in% names(batch)){
    # no missing values
    if(any(is.na(batch$Reference))){
      logging::loginfo("Found NA in Reference column")
      return(FALSE)
    }
    is_ref <- batch$Reference!=0
    if(sum(is_ref)>1){
      return(TRUE)
    }else{
      logging::logerror("Require at least two references per batch.")
      return(FALSE)
    }
  }
  return(TRUE)
}

#' Replaces missing values (NaN) by NA, this appears to be faster
#' 
#' @param data The data as dataframe
#' @return The data with the replaced MVs
replace_missing <- function(data){
  data[vapply(data, is.nan, logical(dim(data)[1]))] <- NA
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
  
  if(typeof(data)=="S4"){
    # Summarized Experiment
    logging::loginfo("Recognized input as S4 class - assuming SummarizedExperiment")
    if(length(SummarizedExperiment::assays(data))!=1){
      logging::logerror("BERT only supports batch effect correction for SummarizedExperiments with a single assay.")
      stop()
    }
    logging::loginfo("Typecasting input to dataframe.")
    # obtain raw data from assay with observations in rows and features in columns
    raw_data <- data.frame(t(SummarizedExperiment::assay(data)))
    # obtain batch/label/sample/reference column
    raw_data["Batch"] <- SummarizedExperiment::colData(data)$Batch
    if("Sample" %in% names(SummarizedExperiment::colData(data))){
      raw_data["Sample"] <- SummarizedExperiment::colData(data)$Sample
    }
    if("Label" %in% names(SummarizedExperiment::colData(data))){
      raw_data["Label"] <- SummarizedExperiment::colData(data)$Label
    }
    if("Reference" %in% names(SummarizedExperiment::colData(data))){
      raw_data["Reference"] <- SummarizedExperiment::colData(data)$Reference
    }
    # potential covariables
    cov_names <- names(SummarizedExperiment::colData(data))[grepl( "Cov" , names( SummarizedExperiment::colData(data)  ) )]
    if(length(cov_names)>0){
      for(n in cov_names){
        raw_data[n] <- SummarizedExperiment::colData(data)[n][,1]
      }
    }
    data <-  raw_data
  }
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
    dtypes <- vapply(data[, all_names], typeof, character(1))
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
    # require at least two references per batch
    if(!verify_references(data_batch)){
      logging::logerror(paste("Reference column error in batch", b))
      stop()
    }
  }
  # count missing values
  final_mvs <- sum(is.na(data))
  
  logging::loginfo(paste("Introduced ", final_mvs-inital_mvs, " missing values due to singular proteins at batch/covariate level."))
  
  logging::loginfo("Done")
  
  return(data)
}