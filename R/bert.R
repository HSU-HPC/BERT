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
        stop("Unrecognized backend.")
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
#'@param combatmode The mode to use for combat (ignored if limma).
#'Encoded options 'are the same as for HarmonizR
#'@param backend The backend to choose for communicating the data, Valid choices
#'are "default" and "file". The latter will use temp files for communicating
#'data chunks between the processes.
#'@return dataframe with the adjusted matrix
parallel_bert <- function(
        chunks, 
        method="ComBat", 
        combatmode=1,
        backend = "default"){
    
    # parallel adjustment as far as possible for this chunk
    biocbackend <- BiocParallel::bpparam()
    result <- BiocParallel::bplapply(chunks, function(chunk){
        if(backend=="file"){
            is_rank_1 <- FALSE#(chunk==chunks[1])
            # read dataframe containing the adjusted data
            data <- readRDS(chunk)
        }else if (backend=="default"){
            data <- chunk
            is_rank_1 <- FALSE # deactivates logging on all processes
        }else{
            stop("Unrecognized communication backend.")
        }
        # split data and covariates
        mod <- data.frame(data [ , grepl( "Cov" , names( data  ) ) ])
        data <- data [ , !grepl( "Cov" , names( data  ) ) ]
        # don't allow covariables AND references
        if((ncol(mod)) & ("Reference" %in% names(data))){
            error_str <- paste("Covariable and reference columns should",
                               "not exist simultanously.")
            stop(error_str)
        }
        
        # number of batches at current level
        num_batches <- length(unique(data$Batch))
        
        hierarchy_counter <- 1
        while(num_batches>1){
            if(is_rank_1){
                logging::loginfo(paste(
                    "Worker one is processing hierarchy level",
                    hierarchy_counter, "with", num_batches))
            }
            data <- adjustment_step(data, mod, combatmode, method)
            
            num_batches <- length(unique(data$Batch))
            # increase hierarchy level
            hierarchy_counter <- hierarchy_counter + 1
        }
        
        # append covariates again
        if(ncol(mod)>1){
            data <- cbind(data, mod)
        }else{
            # if we have only one covariable, its name will have changed --> do
            # it manually
            data["Cov_1"] <- mod[[colnames(mod)[1]]]
        }
        # return the adjusted split
        data
    }, BPPARAM=biocbackend)
    
    adjusted_data <- do.call("rbind", result)
    
    # adjusted the assembled data with all the adjusted chunks
    return(adjusted_data)
}

#' @param data Matrix dataframe/SummarizedExperiment in the format (samples,
#' features). 
#' Additional column names are "Batch", "Cov_X" (were X may be any number),
#' "Label", "Sample" and "Reference". Must contain at least two features.
#' @param cores The number of cores to use for parallel adjustment. Increasing
#' this number leads to faster adjustment, especially on Linux machines. The
#' default is 1.
#' @param bpparameters Optional, default is NULL. If given, this should be
#' a BiocParallel BiocParallelParam instance. Note, that BERT will register
#' this instance as the default
#' @param combatmode Integer, encoding the parameters to use for ComBat.
#' 1 (default)    par.prior = TRUE, mean.only = FALSE
#' 2              par.prior = TRUE, mean.only = TRUE
#' 3              par.prior = FALSE, mean.only = FALSE
#' 4              par.prior = FALSE, mean.only = TRUE
#' Will be ignored, if method!="ComBat".
#' @param corereduction Reducing the number of workers by at least this number
#' @param stopParBatches The minimum number of batches required at a hierarchy
#' level to proceed with parallelized adjustment. If the number of batches
#' is smaller, adjustment will be performed sequentially to avoid overheads.
#' @param backend The backend to choose for communicating the data. 
#' Valid choices are "default" and "file". The latter will use temp files for
#' communicating data chunks between the processes. after adjusting all
#' sub-trees as far as possible with the previous number of cores.
#' @param method Adjustment method to use. Should either be "ComBat", "limma"
#' or "ref". Also allows "None" for testing purposes, which will perform no BE
#' adjustment
#' @param qualitycontrol Boolean indicating, whether ASWs should be computed
#' before and after batch effect adjustment. If TRUE, will compute ASW with
#' respect to the "Batch" and "Label" column (if existent).
#' @param verify Whether the input matrix/dataframe needs to be verified before
#' adjustment (faster if FALSE)
#' @param labelname A string containing the name of the column to use as class
#' labels. The default is "Label".
#' @param batchname A string containing the name of the column to use as batch
#' labels. The default is "Batch".
#' @param referencename A string containing the name of the column to use as ref.
#' labels. The default is "Reference".
#' @param samplename A string containing the name of the column to use as sample
#' name. The default is "Sample".
#' @param covariatename A vector containing the names of columns with
#' categorical covariables. The default is NULL, for which all columns with
#' the pattern "Cov" will be selected.
#' @param assayname User-defined string that specifies, which assay to select,
#' if the input data is a SummarizedExperiment. The default is NULL.
#' @return None. Will instead throw an error, if input is not as intended.
validate_bert_input <- function(data, cores, combatmode,
                                corereduction, stopParBatches, backend, method,
                                qualitycontrol, verify, labelname, batchname, 
                                referencename, samplename, covariatename,
                                assayname){
    
    # data should have at least two features
    stopifnot("data should be dataframe, matrix or SummarizedExperiment"={
        is.matrix(data) || is.data.frame(data) || 
            methods::is(data, "SummarizedExperiment")
    })
    stopifnot("cores should be integer >=1 or NULL" = {
        (cores%%1==0 && cores>0) || is.null(cores)
        })
    stopifnot("combatmode should be in c(1,2,3,4)"={
        combatmode %in% c(1,2,3,4)
    })
    stopifnot("combatmode should be in c(1,2,3,4)"={
        (combatmode %in% c(1,2,3,4))&&is.numeric(combatmode)
    })
    stopifnot("corereduction should be integer >=1"={
        corereduction%%1==0 && corereduction>0
    })
    stopifnot("stopParBatches should be integer >=1"={
        stopParBatches%%1==0 && stopParBatches>0
    })
    stopifnot("backend should be in c(\"default\", \"file\")"={
        backend %in% c("default", "file")
    })
    stopifnot("method should be in c(\"ComBat\", \"limma\", \"None\", \"ref\")"={
        method %in% c("ComBat", "limma", "None", "ref")
    })
    stopifnot("qualitycontrol should be either TRUE or FALSE" = {
        qualitycontrol %in% c(TRUE, FALSE)
    })
    stopifnot("verify should be either TRUE or FALSE" = {
        verify %in% c(TRUE, FALSE)
    })
    stopifnot("labelname must be string with length >=1" = {
        is.character(labelname) && nchar(labelname)>0
    })
    stopifnot("batchname must be string with length >=1" = {
        is.character(batchname) && nchar(batchname)>0
    })
    stopifnot("referencename must be string with length >=1" = {
        is.character(referencename) && nchar(referencename)>0
    })
    stopifnot("samplename must be string with length >=1" = {
        is.character(samplename) && nchar(samplename)>0
    })
    stopifnot("covariatename must be string with length >=1 or NULL" = {
        (is.character(covariatename) && nchar(covariatename)>0) ||
            is.null(covariatename)
    })
    stopifnot("assayname must be string with length >=1 or NULL" = {
        (is.character(assayname) && nchar(assayname)>0) ||
            is.null(assayname)
    })
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
#' @param data Matrix dataframe/SummarizedExperiment in the format (samples,
#' features). 
#' Additional column names are "Batch", "Cov_X" (were X may be any number),
#' "Label", "Sample" and "Reference". Must contain at least two features.
#' @param cores The number of cores to use for parallel adjustment. Increasing
#' this number leads to faster adjustment, especially on Linux machines. The
#' default is NULL, in which case the BiocParallel::bpparam() backend will be
#' used. If an integer is given, a backend with the corresponding number
#' of workers will be created and registered as default for usage.
#' @param combatmode Integer, encoding the parameters to use for ComBat.
#' 1 (default)    par.prior = TRUE, mean.only = FALSE
#' 2              par.prior = TRUE, mean.only = TRUE
#' 3              par.prior = FALSE, mean.only = FALSE
#' 4              par.prior = FALSE, mean.only = TRUE
#' Will be ignored, if method!="ComBat".
#' @param corereduction Reducing the number of workers by at least this number.
#' Only used if cores is an integer.
#' @param stopParBatches The minimum number of batches required at a hierarchy
#' level to proceed with parallelized adjustment. If the number of batches
#' is smaller, adjustment will be performed sequentially to avoid overheads.
#' @param backend The backend to choose for communicating the data. 
#' Valid choices are "default" and "file". The latter will use temp files for
#' communicating data chunks between the processes. after adjusting all
#' sub-trees as far as possible with the previous number of cores.
#' @param method Adjustment method to use. Should either be "ComBat", "limma"
#' or "ref". Also allows "None" for testing purposes, which will perform no BE
#' adjustment
#' @param qualitycontrol Boolean indicating, whether ASWs should be computed
#' before and after batch effect adjustment. If TRUE, will compute ASW with
#' respect to the "Batch" and "Label" column (if existent).
#' @param verify Whether the input matrix/dataframe needs to be verified before
#' adjustment (faster if FALSE)
#' @param labelname A string containing the name of the column to use as class
#' labels. The default is "Label".
#' @param batchname A string containing the name of the column to use as batch
#' labels. The default is "Batch".
#' @param referencename A string containing the name of the column to use as ref.
#' labels. The default is "Reference".
#' @param samplename A string containing the name of the column to use as sample
#' name. The default is "Sample".
#' @param covariatename A vector containing the names of columns with
#' categorical covariables. The default is NULL, for which all columns with
#' the pattern "Cov" will be selected.
#' @param assayname User-defined string that specifies, which assay to select,
#' if the input data is a SummarizedExperiment. The default is NULL.
#' @return A matrix/dataframe/SummarizedExperiment mirroring the shape of the
#' input. The data will be batch-effect adjusted by BERT.
#' @examples
#' # generate dataset with 1000 features, 5 batches, 10 samples per batch and
#' # two genotypes
#' data = generate_dataset(1000,5,10,0.1, 2)
#' corrected = BERT(data, cores=2)
#' @export
BERT <- function(
        data, 
        cores = NULL,
        combatmode = 1, 
        corereduction=4,
        stopParBatches=2,
        backend="default",
        method="ComBat",
        qualitycontrol=TRUE, 
        verify=TRUE,
        labelname="Label",
        batchname="Batch",
        referencename="Reference",
        samplename="Sample",
        covariatename=NULL,
        assayname=NULL){
    
    # dummy code to suppress bioccheck warning
    typeof(BiocStyle::html_document)
    
    validate_bert_input(data, cores, combatmode, corereduction,
                        stopParBatches, backend, method, qualitycontrol,
                        verify, labelname, batchname, referencename,
                        samplename, covariatename, assayname)
    
    
    
    # measure starting time
    total_start <- Sys.time()
    
    # if SummarizedExperiment, we want to store the original input to preserve
    #all metadata
    if(methods::is(data, "SummarizedExperiment") ){
        original_data <- data
    }
    
    # format dataframe
    if(verify){
        data <- format_DF(data, labelname, batchname, referencename, samplename,
                          covariatename, assayname)
    }else{
        logging::loginfo("Skipping initial DF formatting")
    }
    
    # compute ASWs, if required
    if(qualitycontrol){
        logging::loginfo(paste("Acquiring quality metrics before",
        "batch effect correction."))
        asws_prior <- compute_asw(data)
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
    num_batches <- nrow(unique(data["Batch"]))
    logging::loginfo(paste("Found ", num_batches, " batches."))
    
    user_defined_backend <- FALSE
    
    if(is.null(cores)){
        logging::loginfo(paste("Cores argument is not defined.",
                               "Defaulting to BiocParallel::bpparam().",
                               "Argumens corereduction",
                               "will not be used."))
        user_defined_backend <- TRUE
    }else{
        if(cores==1){
            bpparameters <- BiocParallel::SerialParam()
            logging::loginfo("Set up sequential backend")
        }else{
            if(.Platform$OS.type == "windows"){
                bpparameters <- BiocParallel::SnowParam(workers = cores)
            }else{
                bpparameters <- BiocParallel::MulticoreParam(workers = cores)
            }
            logging::loginfo(paste("Set up parallel execution backend with",
                                   cores, "workers"))
        }
        BiocParallel::register(bpparameters, default = TRUE)
        user_defined_backend <- FALSE
    }
    
    sub_tree_counter <- 1
    while(((num_batches>stopParBatches)&&is.null(cores))||
          ((num_batches>stopParBatches)&&(cores>1))){
        
        if(!is.null(cores)){
            logging::loginfo(paste(
                "Processing subtree level",
                sub_tree_counter,"with",
                num_batches,"batches using",cores,"cores."))
        }else{
            logging::loginfo(paste(
                "Processing subtree level",
                sub_tree_counter))
        }
        
        
        chunks <- chunk_data(data, cores, backend = backend)
        
        data <- parallel_bert(
            chunks, 
            method=method, 
            combatmode=combatmode,
            backend = backend)
        
        # the number of batches that remain in the adjusted data
        num_batches <- nrow(unique(data["Batch"]))
        # we need to lower the number of cores, since the n_cores chunks have
        # already been adjusted as far as possible
        
        if((!user_defined_backend)){
            cores <- max(1, floor(cores/corereduction))
            if(.Platform$OS.type == "windows"){
                bpparameters <- BiocParallel::SnowParam(workers = cores)
            }else{
                bpparameters <- BiocParallel::MulticoreParam(workers = cores)
            }
            BiocParallel::register(bpparameters, default = TRUE)
        }
        
        sub_tree_counter <- sub_tree_counter+1
    }
    
    logging::loginfo(paste(
        "Adjusting the last", num_batches,
        "batches sequentially"))
    # last few batches are adjusted sequentially to avoid overheads
    mod <- data.frame(data [ , grepl( "Cov" , names( data  ) ) ])
    data <- data [ , !grepl( "Cov" , names( data  ) ) ]
    # don't allow covariables AND references
    if((ncol(mod)) & ("Reference" %in% names(data))){
        error_str <- paste(
            "Covariable and reference columns should",
            "not exist simultanously.")
        stop(error_str)
    }
    hierarchy_level <- 1
    while (num_batches > 1) {
        logging::loginfo(paste(
            "Adjusting sequential tree level",
            hierarchy_level, "with", num_batches,
            "batches"))
        data <- adjustment_step(data, mod, combatmode, method)
        # re-count the batches 
        num_batches <- nrow(unique(data["Batch"]))
        hierarchy_level <- hierarchy_level+1
    }
    
    # re-set the batch to the original values
    data[original_rownames, "Batch"] <- original_batches
    
    # append covariates again
    if(ncol(mod)>1){
        data <- cbind(data, mod)
    }else{
        # if we have only one covariable, its name will have changed --> do
        # it manually
        data["Cov_1"] <- mod[[colnames(mod)[1]]]
    }
    
    # end timing measurement for the adjustment
    adjustment_end <- Sys.time()
    logging::loginfo("Done")
    
    # compute ASWs, if required
    if(qualitycontrol){
        logging::loginfo(paste(
            "Acquiring quality metrics after batch effect",
            "correction."))
        asws_after <- compute_asw(data)
        
        # batch information
        if(!is.na(asws_prior$Batch)){
            logging::loginfo(paste(
                "ASW Batch was", asws_prior$Batch,
                "prior to batch effect",
                "correction and is now",
                asws_after$Batch,"."))
        }
        # label information
        if(!is.na(asws_prior$Label)){
            logging::loginfo(paste(
                "ASW Label was", asws_prior$Label,
                "prior to batch effect correction",
                "and is now", asws_after$Label,"."))
        }
    }
    
    # rename again
    colnames(data)[colnames(data)=="Batch"] <- batchname
    colnames(data)[colnames(data)=="Label"] <- labelname
    colnames(data)[colnames(data)=="Sample"] <- samplename
    
    for(x in covariatename){
        colnames(data)[colnames(data)==paste("Cov_",x, sep="")] <- x
    }
    
    # if SummarizedExperiment, return as such object as well
    if(methods::is(data, "SummarizedExperiment") ){
        value <- t(as.matrix(data))
        rownames(value) <- NULL
        colnames(value) <- NULL
        SummarizedExperiment::assay(original_data, paste(assayname, "BERTcorrected", sep="_")) <- value
        data <- original_data
    }
    
    # stop total timing measurement
    total_end <- Sys.time()
    
    
    # get execution time in seconds
    a1 <- as.POSIXct(total_end,origin = "1970-01-01")
    e1 <- as.POSIXct(total_start,origin = "1970-01-01")
    execution_time <- as.numeric(a1) - as.numeric(e1)
    # get adjustment time in seconds
    a2 <- as.POSIXct(adjustment_end,origin = "1970-01-01")
    e2 <- as.POSIXct(adjustment_start,origin = "1970-01-01")
    adjustment_time <- as.numeric(a2) - as.numeric(e2)
    
    frac <- round(
        100 * as.numeric(adjustment_time) / as.numeric(execution_time),
        digits = 2)
    
    logging::loginfo(paste(
        "Total function execution time is ",
        execution_time," s and adjustment time is ",
        adjustment_time,"s (",frac,")"))
    
    return(data)
    
}