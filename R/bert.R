#'Chunks data into n segments with (close-to) equivalent number of batches
#'and stores them in temporary RDS files
#'
#'@param data Dataframe with the data to adjust
#'@param n The number of chunks to create
#'@param backend The backend to choose for communicating the data, Valid choices
#'are "default" and "file". The latter will use temp files for communicating
#'data chunks between the processes.
#'@return Vector with the absolute paths to the temporary files, where the data
#'is stored
chunk_data <- function(data, n, backend="default"){
    unique_batches <- unique(data$"Batch")
    num_batches <- length(unique_batches)
    # we can have at most floor(n_batches/2) chunks
    n <- min(floor(num_batches/2), n)
    splits <- split(unique_batches, sort(unique_batches%%n))
    if(backend=="default"){
        chunks <- list()
    }else if (backend=="file"){
        chunks <- c()
    }else{
        logging::logerror("Unrecognized backend.")
        stop()
    }
    counter <- 1
    for(split in splits){
        # select sub-df
        data_sub <- data[(data$"Batch")%in%split,]
        if (backend=="file"){
            # save to file
            tempfile_name <- paste(tempfile(),".rds", sep="")
            saveRDS(data_sub, tempfile_name)
            chunks <- c(chunks, tempfile_name)
        }else{
            # mus be default backend
            chunks[[counter]] <- data_sub
        }
        counter <- counter+1
    }
    return(chunks)
}

#'
#'Adjusts all chunks of data (in parallel) as far as possible.
#'
#'@param chunks vector with the filenames to the temp files where the
#'sub-matrices are stored
#'@param method the BE-correction method to use. Possible choices are ComBat
#'and limma
#'@param combatmode The mode to use for combat (ignored if limma). Encoded options
#'are the same as for HarmonizR
#'@param backend The backend to choose for communicating the data, Valid choices
#'are "default" and "file". The latter will use temp files for communicating
#'data chunks between the processes.
#'@return dataframe with the adjusted matrix
parallel_bert <- function(chunks, method="ComBat", combatmode=1, backend="default"){
    `%dopar%` <- foreach::`%dopar%`
    chunk <- NULL
    # parallel adjustment as far as possible for this chunk
    adjusted_data <- foreach::foreach(chunk=iterators::iter(chunks), .combine = rbind, .export = "adjustment_step") %dopar% {
        if(backend=="file"){
            is_rank_1 <- (chunk==chunks[1])
            # read dataframe containing the adjusted data
            data <- readRDS(chunk)
        }else if (backend=="default"){
            data <- chunk
            is_rank_1 <- FALSE # deactivates logging on all processes
        }else{
            logging::logerror("Unrecognized communication backend.")
            stop()
        }
        
        # split data and covariates
        mod <- data.frame(data [ , grepl( "Cov" , names( data  ) ) ])
        data <- data [ , !grepl( "Cov" , names( data  ) ) ]
        # don't allow covariables AND references
        if((dim(mod)[2]) & ("Reference" %in% names(data))){
            logging::logerror("Covariable and reference columns should not exist simultanously.")
            stop()
        }
        
        # number of batches at current level
        num_batches <- length(unique(data$Batch))
        
        hierarchy_counter <- 1
        while(num_batches>1){
            if(is_rank_1){
                logging::loginfo(paste("Worker one is processing hierarchy level", hierarchy_counter, "with", num_batches))
            }
            data <- adjustment_step(data, mod, combatmode, method)
            
            # new number of batches
            num_batches <- length(unique(data$Batch))
            # increase hierarchy level
            hierarchy_counter <- hierarchy_counter + 1
        }
        
        # append covariates again
        if(dim(mod)[2]>1){
            data <- cbind(data, mod)
        }else{
            # if we have only one covariable, its name will have changed --> do
            # it manually
            data["Cov_1"] <- mod[[colnames(mod)[1]]]
        }
        # return the adjusted split
        data
    }
    
    # adjusted the assembled data with all the adjusted chunks
    return(adjusted_data)
}

#' Adjust data using the BERT algorithm.
#'
#' This function uses the hierarchical BERT algorithm to adjust data with 
#' batch effects. It assumes that the data is in the format
#' (samples, features) and that missing values are indicated by NA.
#' An additional column labelled "Batch" should indicate the batch. Furthermore
#' all columns named "Cov_1", "Cov_2", ... will be considered as covariate
#' for adjustment. Columns labelled "Label" and "Sample" will be ignored,
#' all other columns are assumed to contain data.
#'
#' @param data Matrix or dataframe in the format (samples, features). 
#' Additional column names are "Batch", "Cov_X" (were X may be any number),
#' "Label", "Sample" and "Reference".
#' @param cores The number of cores to use for parallel adjustment. Increasing
#' this number leads to faster adjustment, especially on Linux machines. The
#' default is 1.
#' @param combatmode Integer, encoding the parameters to use for ComBat.
#' 1 (default)    par.prior = TRUE, mean.only = FALSE
#' 2              par.prior = TRUE, mean.only = TRUE
#' 3              par.prior = FALSE, mean.only = FALSE
#' 4              par.prior = FALSE, mean.only = TRUE
#' Will be ignored, if method!="ComBat".
#' @param method Adjustment method to use. Should either be "ComBat", "limma"
#' or "ref". Also allows "None" for testing purposes, which will perform no BE adjustment
#' @param qualitycontrol Boolean indicating, whether ASWs should be computed before
#' and after batch effect adjustment. If TRUE, will compute ASW with respect to
#' the "Batch" and "Label" column (if existent).
#' @param verify Whether the input matrix/dataframe needs to be verified befire adjustment
#' (faster if FALSE)
#' @param mpi Whether to use MPI for parallelization.
#' @param stopParBatches The minimum number of batches required at a hierarchy level
#' to proceed with parallelized adjustment. If the number of batches
#' is smaller, adjustment will be performed sequentially to avoid overheads.
#' @param corereduction Reducing the number of workers by at least this number
#' @param backend The backend to choose for communicating the data, Valid choices
#' are "default" and "file". The latter will use temp files for communicating
#' data chunks between the processes.
#' after adjusting all sub-trees as far as possible with the previous number of cores.
#' @return A matrix/dataframe mirroring the shape of the input. The data will
#' be batch-effect adjusted by BERT.
#' @examples
#' # generate dataset wiith 1000 features, 5 batches, 10 samples per batch and
#' # two genotypes
#' data = generateDataset(1000,5,10,0.1, 2)
#' corrected = BERT(data)
#' @export
BERT <- function(data, cores = 1, combatmode = 1, method="ComBat", qualitycontrol=TRUE, verify=TRUE, mpi=FALSE, stopParBatches = 4, corereduction=2, backend="default"){
    # store original cores
    original_cores <- cores
    
    # measure starting time
    total_start <- Sys.time()
    
    # if SummarizedExperiment, we want to store the original input to preserve all metadata
    if(typeof(data)=="S4"){
        original_data <- data
    }
    
    # format dataframe
    if(verify){
        data <- format_DF(data)
    }else{
        logging::loginfo("Skipping initial DF formatting")
    }
    
    # compute ASWs, if required
    if(qualitycontrol){
        logging::loginfo("Acquiring quality metrics before batch effect correction.")
        asws_prior <- compute_asw(data)
    }
    
    if(mpi){
        logging::loginfo("Starting MPI cluster.")
        cl <- doMPI::startMPIcluster()
        doMPI::registerDoMPI(cl)
        logging::loginfo("Done")
    } else if (cores > 1) {
        # recommended only on linux
        logging::loginfo(paste("Setting up cluster with ", cores, " cores."))
        if (.Platform$OS.type == "windows") {
            # set up cluster
            cl <- parallel::makeCluster(cores)
            logging::loginfo("Identified OS as Windows. Using Parallel Socket Cluster (PSOCK).")
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
    
    sub_tree_counter <- 1
    while((num_batches>stopParBatches)&&(cores>1)){
        logging::loginfo(paste("Processing subtree level", sub_tree_counter,"with", num_batches,"batches using",cores,"cores."))
        chunks <- chunk_data(data, cores, backend = backend)
        data <- parallel_bert(chunks, method=method, combatmode=combatmode, backend = backend)
        
        # the number of batches that remain in the adjusted data
        num_batches <- dim(unique(data["Batch"]))[1]
        # we need to lower the number of cores, since the n_cores chunks have
        # already been adjusted as far as possible
        cores <- max(1, floor(cores/corereduction))
        
        sub_tree_counter <- sub_tree_counter+1
    }
    
    logging::loginfo(paste("Adjusting the last", num_batches, "batches sequentially"))
    # last few batches are adjusted sequentially to avoid overheads
    mod <- data.frame(data [ , grepl( "Cov" , names( data  ) ) ])
    data <- data [ , !grepl( "Cov" , names( data  ) ) ]
    # don't allow covariables AND references
    if((dim(mod)[2]) & ("Reference" %in% names(data))){
        logging::logerror("Covariable and reference columns should not exist simultanously.")
        stop()
    }
    hierarchy_level <- 1
    while (num_batches > 1) {
        logging::loginfo(paste("Adjusting sequential tree level", hierarchy_level, "with", num_batches, "batches"))
        data <- adjustment_step(data, mod, combatmode, method)
        # re-count the batches 
        num_batches <- dim(unique(data["Batch"]))[1]
        hierarchy_level <- hierarchy_level+1
    }
    
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
    if(mpi){
        doMPI::closeCluster(cl)
        Rmpi::mpi.finalize()
        logging::loginfo("Done")
    }else if (original_cores > 1) {
        logging::loginfo("Stopping cluster gracefully.")
        parallel::stopCluster(cl)
        logging::loginfo("Done")
    }
    # compute ASWs, if required
    if(qualitycontrol){
        logging::loginfo("Acquiring quality metrics after batch effect correction.")
        asws_after <- compute_asw(data)
        
        # batch information
        if(!is.na(asws_prior$Batch)){
            logging::loginfo(paste("ASW Batch was", asws_prior$Batch, "prior to batch effect correction and is now", asws_after$Batch,"."))
        }
        # label information
        if(!is.na(asws_prior$Label)){
            logging::loginfo(paste("ASW Label was", asws_prior$Label, "prior to batch effect correction and is now", asws_after$Label,"."))
        }
    }
    
    # if SummarizedExperiment, return as such object as well
    if(typeof(data)=="S4"){
        value <- t(as.matrix(data))
        rownames(value) <- NULL
        colnames(value) <- NULL
        SummarizedExperiment::assay(original_data) <- value
        data <- original_data
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