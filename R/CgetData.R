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
   
   #. 22/12/2019 -> No cal, el fitxer ja s'ha convertit prèviament, nomès cal obrir el fitxer  .# 
   #. 22/12/2019 .# file_output <- CGetDatafromFile(filename)
   #. 22/12/2019 .# gds<- seqOpen(file_output)
   gds<- seqOpen(filename)
   #. 22/12/2019 -> Fi modif .# 
   
   
   # Get gds filtered data by minmaf
   seqSetFilterCond(gds, maf=minmaf, .progress = FALSE, verbose = TRUE)
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


# CgetData <- function(filename, start = 1, end = NULL, minmaf = 0.01) {
#    
#    # filename = "inst/extdata/example.vcf"
#    file_output_vcf <- CGetDatafromFile(filename)
#    gds<- GdsGenotypeReader(file_output_vcf)
#    
#    # get scanAnnot
#    if ( hasVariable(gds, "sample.annot/sex") == TRUE ) {
#       sex <- getVariable(gds, "sample.annot/sex") 
#    }else {  sex <- rep(NA, nscan(gds)) }
#    
#    scanAnnot <- ScanAnnotationDataFrame(data.frame(getScanID(gds), sex, stringsAsFactors=FALSE))
#    
#    # get snpAnnot
#    snpAnnot <- SnpAnnotationDataFrame(data.frame(getSnpID(gds), as(getChromosome(gds), "integer"),  getPosition(gds),
#                                                  getVariable(gds, "snp.rs.id"), getAlleleA(gds), getAlleleB(gds),
#                                                  stringsAsFactors = FALSE))
#    
#    # Create genotype data
#    genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
#    
#    # get allele frequencies
#    snpfreq <- alleleFrequency(genoData)
#    
#    # get SNPMatrix
#    snpmat <- asSnpMatrix(genoData, snpNames="rsID", scanNames="scanID")
#    
#    # Filter data by minmaf
#    snpfiltrada <- as(snpmat, "character")[,snpfreq[,"MAF"] > minmaf]
#    
#    #get genotype matrix
#    geno.snpmatrix <- CgdsSNPpairMatrix(snpfiltrada)
#    
#    
#    snpsData <- list(genotypes = geno.snpmatrix, map = data.frame(name = colnames(geno.snpmatrix)))
#    snpsData$map$position <- getPosition(gds)[snpfreq[,"MAF"] > minmaf]
#    snpsData$map$chromosome <- getChromosome(gds)[snpfreq[,"MAF"] > minmaf]
#    rownames(snpsData$map) <- colnames(geno.snpmatrix)
#    
#    snpgdsClose(gds);
#    return(snpsData)
#    
# }



#..# showfile.gds(closeall=TRUE)
#..# CgetData("inst/extdata/example.vcf")

