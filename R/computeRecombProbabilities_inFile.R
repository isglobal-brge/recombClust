#' Main function to run all the recombClust pipeline with data in a file
#'
#' This function takes as input phased chromosomes and returns the responsibilities
#' PCA and the cluster classification of the samples.
#'
#' @export
#'
#' @param filename File with genotypic data (SNPs and samples) the data file 
#' type must be .vcf (Variant Call Format), .gds (Genomic Data Structure) or 
#' .bed (Browser Extensible Data).
#' @param annot GenomicRanges with the SNPs annotation
#' @param range Region GenomicRanges
#' @param samples optional, character vector with sample ids to be used to compute models
#' @param window only to debug, by default window = 500bp
#' @param minmaf optional, threshold to filter by minimum MAF, by default minmaf = 0.1
#' @param resfilename string, path and file name to store results. By default
#' this runLDmixture stores data in recombClustResults.hdf5, if resfilename 
#' is empty then function works on memory. To work with big regions it is recommended 
#' to work with data on disk. Only use resfilename = "" when regions are small in order to 
#' prevent memory overflows and low performance.
#' @param resgroup string, folder inside a file where we want to store results. By default
#' runLDmixture stores data inside group "AllModels" in resfilename file.
#' @param overwrite boolean, (optional) either a logical value indicating whether
#' the output file can be overwritten or not, by default files are not overwritten.
#' @param ... Further arguments passed to runLDmixtureModel
#' @return probTab: Summarized probabilities by window
computeRecombProbabilities_inFile <- function(filename, range, samples = NULL,
                                              window = 500, minmaf = 0.1,
                                              resfilename = "recombClustResults.hdf5", 
                                              resgroup = "AllModels",
                                              overwrite = FALSE, ...) 
{
   
   if(is.null(resfilename)){ resfilename <- ""; }
   
   # initialize values
   gcstart <- start(range)
   gcend <- end(range)
   
   # If necessary --> Converts vcf or bed file to gds file
   gds_file <- CGetDatafromFile(filename)
   
   # Read data 
   snpsData <- getData(gds_file, range, samples, minmaf)
   
   
   if( is.list(snpsData) ){
      
      # get genomic coordinates position in file
      pos.gcoord <- CgetIndexfromGenCoord(gds_file, seqlevels(range), start(range), end(range), 0.1)
      
      if(!is.na(pos.gcoord$SNPs[1])) {
      
         annot <- makeGRangesFromDataFrame(snpsData$map, start.field = "position", 
                                           end.field = "position")
         
         # Filter pos.gcoord$SNPs in snpsData$genotypes colnames before select
         # data to avoid subscript out of bounds.
         pos.gcoord$SNPs <- pos.gcoord$SNPs[ which(pos.gcoord$SNPs %in% colnames(snpsData$genotypes) )] 
         
         # Get intersection with pos.gcoord$SNPs and annot to filter common
         # data - subscript out of bounds in happlos and annot
         commonrs <- intersect(pos.gcoord$SNPs, colnames(snpsData$genotypes))
         
         # Get models
         # haplos <- snpsData$genotypes[, pos.gcoord$SNPs]
         haplos <- snpsData$genotypes[, commonrs]
      
         results <- computeRecombProbabilities(haplos, annot[which(names(annot) %in% commonrs),], range, resfilename, window = 500, ...)
      } else {
         results <- NA
      }
      
      
   } else {
      
      results <- NA
      
   }

   # # Return models and matrix with summarized probabilities by widow
   return(results)
}
