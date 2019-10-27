#' Main function to run all the recombClust pipeline
#'
#' This function takes a file as input, converts file to gds format
#' and get data ready to be processed
#'
#' @export
#'
#' @param filename character string with file path
#' @param minmaf minimum MAF, default 0.01
#' @return data ready to be processed
CgetData <- function(filename, minmaf = 0.1) {
## Utilitzar al treballar per parts per poder tenir en compte mides de finestra
#..# CgetData <- function(filename, start = 1, end = NULL, minmaf = 0.01) {
   
   file_output <- CGetDatafromFile(filename)
   gds<- SeqArray::seqOpen(file_output)
   
   # Get gds filtered data by minmaf
   SeqArray::seqSetFilterCond(gds, maf=minmaf, .progress = FALSE, verbose = TRUE)
   # Get genotype from filtered data
   geno.gds <- SeqArray::seqGetData(gds, "genotype")
   
   # getting sample and variable names
   allel <- c('_1', '_2')
   dimnames(geno.gds) <- list(allel, SeqArray::seqGetData(gds, "sample.id"), SeqArray::seqGetData(gds, "annotation/id"))
   
   # Transform list (3d-data) to matrix (2d-data)
   geno.snpmatrix <- CTransformtoSampleAlleles(geno.gds, allel, SeqArray::seqGetData(gds, "sample.id"), SeqArray::seqGetData(gds, "annotation/id"))

   # Getting positions from different variables
   snpsData <- list(genotypes = geno.snpmatrix, map = data.frame(name = colnames(geno.snpmatrix)))
   snpsData$map$position <- SeqArray::seqGetData(gds, "position")
   snpsData$map$chromosome <- SeqArray::seqGetData(gds, "chromosome")
   rownames(snpsData$map) <- SeqArray::seqGetData(gds, "annotation/id")
   
   SeqArray::seqClose(gds)
   
   return(snpsData)
   
}

