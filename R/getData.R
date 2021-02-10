#'
#' This function takes a file as input, converts file to gds format
#' and get data ready to be processed
#'
#' @export
#'
#' @param filename character string with file path
#' @param chromosome character or character vector with chromosomes
#' @param range genomic ranges with regions to be read
#' @param minmaf minimum MAF, default 0.01
#' @return data ready to be processed
getData <- function(filename, chromosome, range, minmaf = 0.1) {

   gds <- seqOpen(filename)
   #. 22/12/2019 -> Fi modif .# 
    
   # Filter - by regions
   seqSetFilterChrom(gds, chromosome, from.bp = start(range), to.bp = end(range))
   
   # Filter - by MAF
   minmaf = minmaf + 1e-9 # To force filter to be strictly greater
   seqSetFilterCond(gds, maf=rep(minmaf,length(start(range))), .progress = FALSE, verbose = TRUE)
   
   # Get genotype from filtered data
   geno.gds <- seqGetData(gds, "genotype")
   
   # getting sample and variable names
   allel <- c('_1', '_2')
   dimnames(geno.gds) <- list(allel, seqGetData(gds, "sample.id"), seqGetData(gds, "annotation/id"))
   
   # Transform list (3d-data) to matrix (2d-data)
   geno.snpmatrix <- CTransformtoSampleAlleles(geno.gds, allel, seqGetData(gds, "sample.id"), seqGetData(gds, "annotation/id"))

   # Getting positions from different variables
   snpsData <- list(genotypes = geno.snpmatrix, map = data.frame(name = colnames(geno.snpmatrix)))
   snpsData$map$position <- SeqArray::seqGetData(gds, "position")
   snpsData$map$chromosome <- SeqArray::seqGetData(gds, "chromosome")
   rownames(snpsData$map) <- seqGetData(gds, "annotation/id")
   
   seqClose(gds)
   
   return(snpsData)
   
}

