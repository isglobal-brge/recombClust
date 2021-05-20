#'
#' This function takes a file as input, converts file to gds format
#' and get data ready to be processed
#'
#' @export
#'
#' @param filename character string with file path
#' @param range genomic ranges with regions to be read
#' @param samples optional, character vector with sample ids to be used to compute models
#' @param minmaf minimum MAF, default 0.01
#' @return data ready to be processed
getData <- function(filename, range, samples, minmaf = 0.1) {


   gds <- seqOpen(filename, readonly = TRUE)
   
   
   seqSetFilterChrom(gds, gsub("chr","",unique(as.character(seqnames(range)))), 
                     from.bp = start(range), to.bp = end(range), is.num = TRUE)
   
   # Filter by samples
   if(!is.null(samples)) {
      seqSetFilter(gds, sample.id =samples)
   }
   
   # Filter - by MAF
   minmaf = minmaf + 1e-9 # To force filter to be strictly greater
   
   # Filter by range
   if( sum(countOverlaps(granges(gds), range))>0 ){
      seqSetFilterCond(gds, maf=rep(minmaf,length(start(range))), .progress = FALSE, verbose = TRUE)
   } else {
      seqClose(gds)
      message("No SNP`s in selected range")
      return(NA)
   }
   
   # Filter indels
   variant.id <- seqGetData(gds, "variant.id")
   svn <- variant.id[which(isSNV(gds, biallelic=TRUE))]
   if(!is.null(svn)){
      seqSetFilter(gds, variant.id = svn)
   }else {
      seqClose(gds)
      message("No SNP`s in selected range")
      return(NA)
   }
   
   # Get genotype from filtered data
   geno.gds <- seqGetData(gds, "genotype")
   
   if(!is.null(geno.gds)){
      # getting sample and variable names
      allel <- c('_1', '_2')
      dimnames(geno.gds) <- list(allel, seqGetData(gds, "sample.id"), seqGetData(gds, "annotation/id"))
      
      # Transform list (3d-data) to matrix (2d-data)
      geno.snpmatrix <- CTransformtoSampleAlleles(geno.gds, allel, seqGetData(gds, "sample.id"), seqGetData(gds, "annotation/id"))
      
      # Getting positions from different variables
      snpsData <- list(genotypes = geno.snpmatrix, map = data.frame(name = colnames(geno.snpmatrix)))
      snpsData$map$position <- SeqArray::seqGetData(gds, "position")
      snpsData$map$chromosome <- SeqArray::seqGetData(gds, "chromosome")
      new_rownames <- sapply(seq(1,length(seqGetData(gds, "annotation/id"))), function(x, data) { name = data[x] 
                                                                                                  if(data[x] == "") 
                                                                                                     name = paste0("NA.",x)
                                                                                                  return(name) },
                             data = seqGetData(gds, "annotation/id") )
      rownames(snpsData$map) <- new_rownames
   } else {
      snpsData <- NA
   }

   seqClose(gds)
   return(snpsData)
   
}

