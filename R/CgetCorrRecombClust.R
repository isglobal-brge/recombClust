#' Recombination correlation matrix
#' 
#' Get correlation matrix from the probability of chromosome recombination
#' 
#' @export
#'
#' @param mat Matrix with model probabilities - recombClust output.
#' @param annot Matrix with initial and final coordinades of the model.
#' @param range range where the function is executed.
#' @param window window size.
#' @return A vector with chromosome probabilities to be in recomb population
CgetCorrRecombClust <- function( probmat, annot, window)
{ 

   # Get probability recombination   
   X <- cGetRecombProb(probmat, annot, window)

   # Gets correlation mattrix
   Xs <- BigDataStatMeth::Normalize_Data(X$res ,bcenter = TRUE,bscale = FALSE)
   covX <- BigDataStatMeth::bdcrossprod(Xs)/(nrow(X$res)-1)
   Sx <- sqrt(diag(covX))
   invSx <- diag(Sx^-1)
   
   corrX <-  BigDataStatMeth::blockmult(BigDataStatMeth::blockmult(invSx,covX),invSx)
   
   # Write_RHDF5
   
   return(list(recombmat = X, corrmat = corrX))
   
}
   