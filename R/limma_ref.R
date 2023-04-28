#' A method to remove batch effects estimated from a subset (references)
#' per batch only.
#' Source code is heavily based on limma::removeBatchEffects by
#' Gordon Smyth and Carolyn de Graaf
#' @param  x the data matrix with samples in columns and features in rows
#' @param batch the batch list as vector.
#' @param references a vector of 0/1, indicating whether the corresponding
#' sample is a reference (1) or not (0)
#' @return the corrected data matrix
#' @export 
removeBatchEffectRefs <- function(x,batch, references)
{
  # convert references to boolean
  isref = references == 1
  # set up initial design for entire data
  batch <- as.factor(batch)
  contrasts(batch) <- contr.sum(levels(batch))
  batch <- model.matrix(~batch)[,-1,drop=FALSE]
  X.batch <- cbind(batch,NULL,NULL)
  # select references
  Xref.batch <- X.batch[isref,]
  designref <- matrix(1,nrow(Xref.batch),1)
  xref <- x[,isref]
  fit <- limma::lmFit(xref, cbind(designref, Xref.batch))
  beta <- fit$coefficients[,-(1:ncol(design)),drop=FALSE]
  beta[is.na(beta)] <- 0
  # now on full data
  as.matrix(x) - beta %*% t(X.batch)
}