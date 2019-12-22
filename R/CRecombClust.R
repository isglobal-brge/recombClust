#' Main function to run all the recombClust pipeline
#'
#' This function takes as input phased chromosomes and returns the responsibilities
#' PCA and the cluster classification of the samples.
#'
#' @export
#'
#' @param filename File with genotypic data (SNPs and samples) the data file type must be .vcf (Variant Call Format), 
#' .gds (Genomic Data Structure) or .bed (Browser Extensible Data).
#' @param start start position to be treated
#' @param end end position to be treated
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
CRecombClust <- function( filename, chromosome, gcstart = NULL, gcend = NULL, clusters = 2, PCs = 3, ...){ 
   
   # Convert datafile and gets the filename string route of new file
   gds_file <- CGetDatafromFile(filename)
   
   # Call CrunRecombclust
   res <- CrunRecombClust( gds_file, chromosome, gcstart, gcend, clusters, PCs)
   
   ## TO DO: create an object to encapsulate results
   return(res)
}




