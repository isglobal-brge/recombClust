#' Main function to run all the recombClust pipeline
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
#' @return A list with four elements:
#' \itemize{
#'  \item{class: Cluster classification of the chromosomes}
#'  \item{pc: Responsibilities PCA}
#'  \item{mat: Responsibilities matrix}
#'  \item{models: List of models}
#' }
runRecombClust <- function(haplos, annot, clusters = 2, PCs = 3, ...){

  # Tareas generales a mejorar
  ## Hacer los tests (http://r-pkgs.had.co.nz/tests.html)
  ### Comprobar un caso básico para que funcione
  ### Comprobar errores comunes (funciones altas)
  ## Aceptar missings (al menos 5% de individuos con missings)
  ### Individuos con missings, tienen un NA como resultados.
  ## Chequear missings al principio de la función (maximo 5% individuos y SNPs)
  ### Error y parar (stop)
  ## Comprobar que los SNPs en haplos son los mismos en annot

  # Get models
  models <- runLDmixtureModel(haplos, annot, ...)


  ## Remove failed models
  goodModels <- vapply(models, class, character(1)) == "list"
  ## Create matrix of chromosome responsibilities
  indsmat <- do.call(cbind, lapply(models[goodModels], `[[`, "r1"))
  rownames(indsmat) <- rownames(haplos)

  ## Run PCA on individuals responsibilities
  pc <- stats::prcomp(indsmat, rank. = PCs)
  rownames(pc$x) <- rownames(haplos)

  ## Get classification with k-means
  class <- stats::kmeans(pc$x[, seq_len(PCs)], centers = clusters, nstart = 1000)$cluster
  names(class) <- rownames(haplos)

  ## TO DO: create an object to encapsulate results
  return(list(class = class, pc = pc, mat = indsmat, models = models))
}
