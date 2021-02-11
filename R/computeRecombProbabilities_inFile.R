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
#' @param @param samples optional, character vector with sample ids to be used to compute models
#' @param window only to debug, by default window = 500bp
#' @param ... Further arguments passed to runLDmixtureModel
#' @return probTab: Summarized probabilities by window
computeRecombProbabilities_inFile <- function(filename, range, samples = NULL, window = 500, minmaf = 0.1, ...) {
   
   # initialize values
   gcstart <- start(range)
   gcend <- end(range)
   
   # If necessary --> Converts vcf or bed file to gds file
   gds_file <- CGetDatafromFile(filename)
   
   # Read data 
   snpsData <- getData(gds_file, range, samples, minmaf)
   
   # get genomic coordinates position in file
   pos.gcoord <- CgetIndexfromGenCoord(gds_file, seqlevels(range), start(range), end(range), 0.1)
   
   annot <- makeGRangesFromDataFrame(snpsData$map, start.field = "position", 
                                     end.field = "position")
   
   # Get models
   haplos <- snpsData$genotypes[, pos.gcoord$SNPs]

   results <- computeRecombProbabilities(haplos, annot, range, window = 500)

   # # Return models and matrix with summarized probabilities by widow
   return(results)
}
