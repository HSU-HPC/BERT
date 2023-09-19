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


#' Compute the average silhouette width (ASW) for the dataset with respect
#' to both label and batch.
#'
#'Columns labelled Batch, Sample, Label, Reference and Cov_1 will be ignored.
#'
#' @param dataset Dataframe in the shape (samples, features) with additional
#' columns Batch and Label.
#' @return List with fields "Label" and "Batch" for the ASW with regards to 
#' Label and Batch respectively.
#' @examples
#' # generate dataset with 1000 features, 5 batches, 10 samples per batch and
#' # two genotypes
#' data = generate_dataset(1000,5,10,0.1, 2)
#' asw = compute_asw(data)
#' asw
#' @export
compute_asw <- function(dataset){
    # check that input is indeed a dataframe.
    if(!is.data.frame(dataset)){
        logging::logerror(paste("Parameter dataset in function compute_asw",
                                "must be of type dataframe."))
        stop()
    }
    
    # create copy in which we can apply ordinal-encoding to Batch and Label
    ds2 <- data.frame(dataset)
    if("Batch" %in% colnames(ds2)){
        ds2$Batch <- ordinal_encode(ds2$Batch)
    }
    if("Label" %in% colnames(ds2)){
        ds2$Label <- ordinal_encode(ds2$Label)
    }
    
    dataset_nocov <- ds2 [ , !grepl( "Cov" , names( ds2  ) ) ]
    # numeric values in dataset only
    num_values <- dataset_nocov[,!names(dataset_nocov) %in% c(
        "Batch", 
        "Sample", 
        "Label", 
        "Cov_1", 
        "Reference")]
    # compute distance matrix based on euclidean distances, ignoring NAs
    distancematrix <- stats::dist(num_values)
    
    if("Label" %in% names(dataset_nocov)){
        # labels as vector
        labels <- as.vector(dataset_nocov[["Label"]])
        # compute silhouette object wrt. labels
        sil_labels <- cluster::silhouette(labels, dist=distancematrix)
        # extract ASW wrt. labels
        asw_label <- summary(sil_labels)["avg.width"]$avg.width
    }else{
        asw_label <- NA
    }
    if("Batch" %in% names(dataset_nocov)){
        # batches as vector
        batches <- as.vector(dataset_nocov[["Batch"]])
        # compute silhouette object wrt. batches
        sil_batches <- cluster::silhouette(batches, dist=distancematrix)
        # extract ASW wrt. batches
        asw_batches <- summary(sil_batches)["avg.width"]$avg.width
    }else{
        asw_batches <- NA
    }
    # create list object
    ret <- list("Label" = asw_label, "Batch" = asw_batches)
    return(ret)
}