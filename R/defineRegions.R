#' Get regions
#'
#' Get detected regions with a minimum correlation and minimum length
#'
#' @export
#'
#' @param probTab Probability matrix by genetic windows
#' @param minR Minimum correlation between windows
#' @param minLength Minimum window length
#' @return Detected regions
#' 
define_Regions <- function(probTab, minR, minLength)
{
   
   # Filter positions with missing data > 5% 
   probTab <- probTab[, which(colMeans(!is.na(probTab)) > 0.05)]
   
   # Get correlation matrix ( Absolute values)
   corProbs <- getCorrelationMatrix(probTab, absval = TRUE);

   # Apply BigLD algorithm
   
}