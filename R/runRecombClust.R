#' Main function to run all the LDclassifier pipeline
#'
#' This function takes as input phased chromosomes and returns the responsibilities
#' PCA and the cluster classification of the samples.
#'
#' @export
#'
#' @param haplos Matrix with the haplotypes (SNPs in rows, samples in columns)
#' @param annot GenomicRanges with the SNPs annotation
#' @param clusters Numeric with the clusters used in k-means
#' @param PCs Numeric with the number of PCA components used to make the clustering.
#' @param ... Further arguments passed to runLDmixtureModel
#' @return A list with two elements:
#' \itemize{
#'  \item{"class"}{Cluster classification of the chromosomes}
#'  \item{"pc"}{Responsibilities PCA}
#' }
runRecombClust <- function(haplos, annot, clusters = 2, PCs = 1, ...){
  # Get models
  models <- runLDmixtureModel(haplos, annot, ...)


  ## Remove failed models
  goodModels <- vapply(models, class, character(1)) == "list"
  ## Create matrix of chromosome responsibilities
  indsmat <- do.call(cbind, lapply(models[goodModels], `[[`, "r1"))

  ## Run PCA on individuals responsibilities
  pc <- stats::prcomp(indsmat)

  ## Get classification with k-means
  class <- stats::kmeans(pc$x[, seq_len(PCs)], centers = clusters, nstart = 1000)$cluster
  names(class) <- colnames(haplos)

  ## TO DO: create an object to encapsulate results
  return(list(class = class, pc = pc))
}
