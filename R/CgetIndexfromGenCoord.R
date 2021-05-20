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
   gds <- seqOpen(filename, readonly = TRUE)

   # Get gds filtered data by minmaf
   seqSetFilterCond(gds, maf=minmaf, .progress = FALSE, verbose = TRUE)

   # Getting index position from Genomic coordinates 
   if( is.null(chrom))
      stop("chromosome not indicated!")
   
   
   ## Tornar els noms enllocs de la posició !!!
   if( is.null(gcstart)) 
      ipos <- min(which(SeqArray::seqGetData(gds, "chromosome")==chrom))
   else
      ipos <- min(intersect( which(SeqArray::seqGetData(gds, "chromosome")==chrom),
                             which(SeqArray::seqGetData(gds, "position")>=gcstart ))) 
   
   if(is.null(gcend))
      epos <- max(which(SeqArray::seqGetData(gds, "chromosome")==chrom))
   else
      epos <- max(intersect( which(SeqArray::seqGetData(gds, "chromosome")==chrom),
                             which(SeqArray::seqGetData(gds, "position")<=gcend ))) 
   
   if( is.null(ipos) || is.null(epos)) {
      message("No data found between start and end genomic coordinates")
      workSNPs <- NA
   } else {
      workSNPs <- SeqArray::seqGetData(gds,"annotation/id")[ipos:epos]
   }
      
   seqClose(gds)
   
   return(list(SNPs = workSNPs, istart = ipos, iend = epos))
}
