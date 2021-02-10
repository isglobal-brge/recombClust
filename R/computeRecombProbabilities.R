#' Main function to run all the recombClust pipeline
#'
#' This function takes as input phased chromosomes and returns the responsibilities
#' PCA and the cluster classification of the samples.
#'
#' @export
#'
#' @param filename File with genotypic data (SNPs and samples) the data file type must be .vcf (Variant Call Format),
#' .gds (Genomic Data Structure) or .bed (Browser Extensible Data).
#' @param annot GenomicRanges with the SNPs annotation
#' @param range Region GenomicRanges
#' @param window only to debug, by default window = 500bp
#' @param ... Further arguments passed to runLDmixtureModel
#' @return probTab: Summarized probabilities by window
computeRecombProbabilities <- function(haplos, annot, range, window = 500, ...) 
{

  # Get models
  models <- runLDmixtureModel(haplos, annot, ...)

  ## Remove failed models
  goodModels <- vapply(models, class, character(1)) == "list"

  # Get matrix
  matProb <- createMatrixfromModels(models[goodModels], range)

  # # Summarize recombClust probabilities by windows
  matWindProb <- getRecombProb(matProb$mat, matProb$annot, range, window)

  #
  # # Return models and matrix with summarized probabilities by widow
  return(probTab = matWindProb)
}
