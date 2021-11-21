#' Get matrix probabilities by window size
#'
#' Get matrix probabilities by window size
#'
#' @export
#'
#' @param filename Results obtained with the different models LDMixture
#' @param group Results obtained with the different models LDMixture
#' @param annot Annotation obtained for the different models in LDMixture (A ELIMINAR !!)
#' @param range Range to look in
#' @param window windows size, by default 500bp
#' @return Probability matrix by window
getRecombProb_hdf5 <- function(filename, group, range, window = 500) {

  ## Define chunks
  starts <- seq(start(range), end(range), window)
  chunks <- GRanges(seqnames = seqnames(range), IRanges(start = starts, width = window))
  
  annot = getAnnotationDataHdf5(filename, group, unique(seqnames(range)))
  annot <- GRanges(seqnames = annot[,1], IRanges(start = annot[,2], end = annot[,3]))
  
  overLaps <- findOverlaps(chunks, annot, type = "within")
  
  ## Compute cluster Recomb freq by mean of voting .
  ## To avoid extra execution we start to compose matrix from first overlap avoiding NAs

  hdf5Dims <- NULL
  
  #. CONDICIÓ BONA !!! .# res <- sapply(seq( from = min(from(overLaps)), to = length(chunks)), function(chunk) {
  res <- sapply(seq( from = min(from(overLaps)), to = min(from(overLaps)) + 1 ), function(chunk) {
    sel <- to(overLaps)[from(overLaps) == chunk]
    
    if (length(sel) == 0) {
        
        if(is.null(hdf5Dims)){
            # Get dataset dimensions, we only get dimension once to avoid overload processes
            hdf5Dims <- get_dimHdf5(filename, paste0(group,"/", runValue(seqnames(chunks[chunk]))))
        }
        
        val <- rep(NA, hdf5Dims[1])
        
    } else {
        
        if(is.null(hdf5Dims)){
            # Get dataset dimensions, we only get dimension once to avoid overload processes
            hdf5Dims <- get_dimHdf5(filename, paste0(group,"/", runValue(seqnames(chunks[chunk]))))
        }
        
        datasetname <- runValue(seqnames(chunks[chunk]))
        getProbs_hdf5(filename, group, "22", sel, hdf5Dims[1])
        
        print("Tornem a estar al punt inicial")
        return(0)
    
      # val <- getProbs(mat, sel)
    }
  })

  # rownames(res) <- rownames(mat)
  # colnames(res) <- paste(start(chunks), end(chunks), sep = "-")
  
  #..# res
}
