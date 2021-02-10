#' Main function to run all the recombClust pipeline with data in a file
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
#' @return A list with four elements:
#' \itemize{
#'  \item{porbTab: Summarized probabilities by window}
#'  \item{models: List of models}
#' }
computeRecombProbabilities_inFile <- function(filename, chromosome, range, window = 500, minmaf = 0.1, ...) {
   
   # initialize values
   gcstart <- start(range)
   gcend <- end(range)
   
   # If necessary --> Converts vcf or bed file to gds file
   gds_file <- CGetDatafromFile(filename)
   
   # Read data 
   snpsData <- getData(gds_file, chromosome, range, minmaf)
   
   # get genomic coordinates position in file
   pos.gcoord <- CgetIndexfromGenCoord(gds_file, chromosome, start(range), end(range), 0.1)
   
   annot <- makeGRangesFromDataFrame(snpsData$map, start.field = "position", 
                                     end.field = "position")
   
   # Get models
   haplos <- snpsData$genotypes[, pos.gcoord$SNPs]

   results <- computeRecombProbabilities(haplos, annot, range, window = 500)

   # # Return models and matrix with summarized probabilities by widow
   return(results)
}
