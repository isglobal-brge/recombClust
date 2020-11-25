#' Get matrix from model
#'
#' Get matrix from models obtained by LDMixture
#'
#' @export
#'
#' @param recombRes Results from LDMixture
#' @param range Region GRanges
#' @param range Range to look in
#' @param window windows size, by default 500bp
#' @return Probability matrix by window
createMatrixfromModels <- function(recombRes, range){
   
   mat <- do.call(cbind, lapply(recombRes, `[[`, "r1"))
   allCoords <- data.frame(t(sapply(recombRes, function(x) x$annot)))
   allCoords$chr <- seqnames(range)
   grMods <- makeGRangesFromDataFrame(allCoords)
   list(mat = mat, annot = grMods)
}
