#'
#' This function get index position from genomic coordinates
#'
#'
#' @param filename character string with file path
#' @param chrom chromosome
#' @param gcstart genomic coordinates - start
#' @param gcend genomic coordinates - end
#' @param minmaf minimum MAF, default 0.01
#' @return data ready to be processed
#' @export
CgetIndexfromGenCoord <- function(filename, chrom, gcstart=NULL, gcend=NULL, minmaf = 0.1) 
{
   
   # Open file
   gds <- seqOpen(filename)

   # Get gds filtered data by minmaf
   seqSetFilterCond(gds, maf=minmaf, .progress = FALSE, verbose = TRUE)

   # Getting index position from Genomic coordinates 
   if( is.null(chrom))
      stop("chromosome not indicated!")
   
   
   ## Tornar els noms enllocs de la posició !!!
   if( is.null(gcstart)) 
      ipos <- min(which(SeqArray::seqGetData(gds, "chromosome")==chrom))
   else
      ipos <- which(SeqArray::seqGetData(gds, "chromosome")==chrom & SeqArray::seqGetData(gds, "position")==gcstart)
   
   if(is.null(gcend))
      epos <- max(which(SeqArray::seqGetData(gds, "chromosome")==chrom))
   else
      epos <- which(SeqArray::seqGetData(gds, "chromosome")==chrom & SeqArray::seqGetData(gds, "position")==gcend)
   
   workSNPs <- SeqArray::seqGetData(gds,"annotation/id")[ipos:epos]
      
   seqClose(gds)
   
   return(list(SNPs = workSNPs, istart = ipos, iend = epos))
}
