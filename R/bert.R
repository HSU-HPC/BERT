#' Adjust data using the BERT algorithm.
#'
#' This function uses the hierarchical BERT algorithm to adjust data with 
#' batch effects. It assumes that the data is in the format
#' (samples, features) and that missing values are indicated by NA.
#' An additional column labelled "Batch" should indicate the batch. Furthermore
#' all columns named \"Cov_1\", \"Cov_2\", ... will be considered as covariate
#' for adjustment. Columns labelled "Label" and "Sample" will be ignored,
#' all other columns are assumed to contain data.
#'
#' @param data Matrix or dataframe in the format (samples, features). 
#' Additional column names are \"Batch\", \"Cov_X\" (were X may be any number),
#' "Label" and "Sample".
#' @param cores The number of cores to use for parallel adjustment. Increasing
#' this number leads to faster adjustment, especially on Linux machines. The
#' default is 1.
#' @param combatmode Integer, encoding the parameters to use for ComBat.
#' 1 (default)    par.prior = TRUE, mean.only = FALSE
#' 2              par.prior = TRUE, mean.only = TRUE
#' 3              par.prior = FALSE, mean.only = FALSE
#' 4              par.prior = FALSE, mean.only = TRUE
#' Will be ignored, if method=="limma".
#' @param method Adjustment method to use. Should either be \"ComBat\" or \"limma\".
#' @return A matrix/dataframe mirroring the shape of the input. The data will
#' be batch-effect adjusted by BERT.
#' @export
hierarchical_adjustment <- function(data, cores = 1, combatmode = 1, method="ComBat") {
  # measure starting time
  total_start <- Sys.time()
  
  # format dataframe
  data <- format_DF(data)
  
  # split covariates from remaining data
  mod <- data.frame(data [ , grepl( "Cov" , names( data  ) ) ])
  data <- data [ , !grepl( "Cov" , names( data  ) ) ]
  logging::loginfo(paste("Found ", dim(mod)[2], "covariates"))
  
  if (cores > 1) {
    # recommended only on linux
    logging::loginfo(paste("Setting up cluster with ", cores, " cores."))
    if (.Platform$OS.type == "windows") {
      # set up cluster
      cl <- parallel::makeCluster(cores)
      logging::loginfo("Identified OS as Windows. Using Parallel Socket Cluster (PSOCK).")
      logging::logwarn("Usage of cores>1 on Windows machines is not recommended.")
    } else{
      # set up cluster with forking.
      cl <- parallel::makeForkCluster(cores)
      logging::loginfo("Identified OS as UNIX (aka not windows). Using forking.")
    }
    # register parallel backend
    doParallel::registerDoParallel(cl)
    logging::loginfo("Done")
  }
  
  # store the original batches, because we need to manually set them again
  # after adjustment
  original_batches <- data[["Batch"]]
  original_rownames <- rownames(data)
  
  logging::loginfo("Starting hierarchical adjustment")
  
  # start timing
  adjustment_start <- Sys.time()
  
  # count the current number of batches. BERT will adjust on new
  # hierarchy levels, as long as there are at least 2 batches
  num_batches <- dim(unique(data["Batch"]))[1]
  logging::loginfo(paste("Found ", num_batches, " batches."))
  
  # counter for the hierarchy levels
  hierarchy_counter <- 1
  
  # repeat, until only one large, adjusted batch remains
  while (num_batches > 1) {
    logging::loginfo(paste("Adjusting Hierarchy Level ", hierarchy_counter, ": ", num_batches, " Batches"))
    # if we use parallelization
    if (cores > 1) {
      # do adjustment step in parallel
      data <- adjustment_step_parallel(data, mod, combatmode, method)
    } else{
      # do it sequentially
      data <- adjustment_step(data, mod, combatmode, method)
      
    }
    # re-count the batches (they will be altered due to the adjustment at the
    # current hierarchy level)
    num_batches <- dim(unique(data["Batch"]))[1]
    # increase hierarchy level
    hierarchy_counter <- hierarchy_counter + 1
  }
  # --- here, adjustment is finished
  
  # re-set the batch to the original values
  data[original_rownames, "Batch"] <- original_batches
  
  # append covariates again
  if(dim(mod)[2]>1){
    data <- cbind(data, mod)
  }else{
    # if we have only one covariable, its name will have changed --> do
    # it manually
    data["Cov_1"] <- mod[[colnames(mod)[1]]]
  }
  # end timing measurement for the adjustment
  adjustment_end <- Sys.time()
  logging::loginfo("Done")
  
  # stop cluster
  if (cores > 1) {
    logging::loginfo("Stopping cluster gracefully.")
    parallel::stopCluster(cl)
    logging::loginfo("Done")
  }
  # stop total timing measurement
  total_end <- Sys.time()
  
  # get execution time in seconds
  execution_time <- as.numeric(as.POSIXct(total_end,origin = "1970-01-01")) - as.numeric(as.POSIXct(total_start,origin = "1970-01-01"))
  # get adjustment time in seconds
  adjustment_time <- as.numeric(as.POSIXct(adjustment_end,origin = "1970-01-01")) - as.numeric(as.POSIXct(adjustment_start,origin = "1970-01-01"))
  
  
  logging::loginfo(paste("Total function execution time is ",execution_time," s and adjustment time is ",adjustment_time,"s (",round(100 * as.numeric(adjustment_time) / as.numeric(execution_time),digits = 2),")"))
  return(data)
  
}