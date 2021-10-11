#' Get the windows within range
#'
#' Get the windows within range
#'
#' @export
#'
#' @param annot GenomicRanges with the SNPs annotation
#' @param range Range to look in
#' @param window windows size, by default 500bp
#' @return windows that overlap with the annotations within the range range
#' \itemize{
#'  \item{"chunk"}{Log-likelihood of mixture model}
#'  \item{"overLaps"}{Log-likelihood of linkage model}
#' }
getWindows <- function(annot, range, window = 500) {
   
   ## Define chunks
   starts <- seq(start(range), end(range), window)
   chunks <- GRanges(seqnames = seqnames(range), IRanges(start = starts, width = window))
   
   overLaps <- findOverlaps(chunks, annot, type = "within")
   
   return(list(chunks = chunks, overLaps = overLaps))
   
}
