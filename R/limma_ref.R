#' Identifies the references to use for this specific batch effect adjustment
#' @param vector of batch numbers. Must contain 2 unique elements
#' @param references vector that contains 0, if the sample is to be c-adjusted
#' and a class otherwise
#' @return the indices of the reference samples
#' @export
identify_references <- function(batch, references){
  # numbers of the two batches <- must be two
  batch_vals <- unique(batch)
  if(length(batch_vals)!=2){
    logging::logerror("removeBatchEffectRefs takes only two batches")
    stop()
  }
  batch1 <- batch_vals[1]
  batch2 <- batch_vals[2]
  # assemble references
  possible_idx <- 1:length(batch)
  ref_idx <- c()
  
  for(cl in unique(references)){
    if(cl!=0){
      t1 <- (references==cl)&(batch==batch1)
      t2 <- (references==cl)&(batch==batch2)
      min_ref <- min(sum(t1), sum(t2))
      if(min_ref>1){
        idx_1 <- sample(possible_idx[t1], size = min_ref, replace = FALSE)
        idx_2 <- sample(possible_idx[t2], size = min_ref, replace = FALSE)
        ref_idx <- c(ref_idx, idx_1)
        ref_idx <- c(ref_idx, idx_2)
      }else if(min_ref==1){
        idx_1 <- possible_idx[t1]
        idx_2 <- possible_idx[t2]
        ref_idx <- c(ref_idx, idx_1)
        ref_idx <- c(ref_idx, idx_2)
      }
    }
  }
  if(length(ref_idx)<2){
    logging::logerror("Did not find enough references of common class.")
    stop()
  }
  # convert references to boolean
  isref <- (references==max(unique(references))+1)
  isref[ref_idx] <- TRUE
  return(isref)
}

#' A method to remove batch effects estimated from a subset (references)
#' per batch only.
#' Source code is heavily based on limma::removeBatchEffects by
#' Gordon Smyth and Carolyn de Graaf
#' @param  x the data matrix with samples in columns and features in rows
#' @param batch the batch list as vector.
#' @param references a vector of integers, indicating whether the corresponding
#' sample is to be co-adjusted (0) or may be used as a reference (>0)
#' @return the corrected data matrix
#' @export 
removeBatchEffectRefs <- function(x,batch,references)
{
  isref <- identify_references(batch, references)
  # set up initial design for entire data
  batch <- as.factor(batch)
  contrasts(batch) <- contr.sum(levels(batch))
  batch <- model.matrix(~batch)[,-1,drop=FALSE]
  X.batch <- cbind(batch,NULL,NULL)
  design <- matrix(1,nrow(X.batch),1)
  # select references
  Xref.batch <- as.matrix(X.batch[isref,])
  designref <- matrix(1,nrow(Xref.batch),1)
  xref <- x[,isref]
  fit <- limma::lmFit(xref, cbind(designref, Xref.batch))
  beta <- fit$coefficients[,-(1:ncol(design)),drop=FALSE]
  beta[is.na(beta)] <- 0
  # now on full data
  as.matrix(x) - beta %*% t(X.batch)
}