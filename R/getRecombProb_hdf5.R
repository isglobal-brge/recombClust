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
#' @param samples Sample names
#' @param window windows size, by default 500bp
#' @return Probability matrix by window
getRecombProb_hdf5 <- function(filename, group, range, samples, window = 500) {

  ## Define chunks
  starts <- seq(start(range), end(range), window)
  chunks <- GRanges(seqnames = seqnames(range), IRanges(start = starts, width = window))
  
  annot = getAnnotationDataHdf5(filename, group, unique(seqnames(range)))
  annot <- GRanges(seqnames = annot[,1], IRanges(start = annot[,2], end = annot[,3]))
  
  overLaps <- findOverlaps(chunks, annot, type = "within")
  
  ## Compute cluster Recomb freq by mean of voting .
  ## To avoid extra execution we start to compose matrix from first overlap avoiding NAs

  hdf5Dims <- NULL
  
  chunkUsed <- sapply(seq( from = min(which(from(overLaps)== min(from(overLaps)))), 
                           to = length(chunks)), function(chunk) {
    sel <- to(overLaps)[from(overLaps) == chunk]
    
    if (length(sel) == 0) {
        
        if(is.null(hdf5Dims)){
            # Get dataset dimensions, we only get dimension once to avoid overload processes
            hdf5Dims <- get_dimHdf5(filename, paste0(group,"/", runValue(seqnames(chunks[chunk]))))
        }
        
        # val <- rep(NA, hdf5Dims[1])
        return(NA) 
    } else {
        
        if(is.null(hdf5Dims)){
            # Get dataset dimensions, we only get dimension once to avoid overload processes
            hdf5Dims <- get_dimHdf5(filename, paste0(group,"/", runValue(seqnames(chunks[chunk]))))
        }
        
        datasetname <- runValue(seqnames(chunks[chunk]))
        getProbs_hdf5(filename, group, as.character(unique(seqnames(range))), sel, hdf5Dims[1] - 4)
        
        return(chunk) 
        
      # val <- getProbs(mat, sel)
    }
  })

  columns_n <- paste(start(chunks[chunkUsed[which(!is.na(chunkUsed))]]) , 
                         end(chunks[chunkUsed[which(!is.na(chunkUsed))]]), sep = "-")
 
  write_dimNames( filename = filename, group = group, dataset = "ModelsProb", 
                   rownames = samples, 
                   colnames = columns_n, 
                   force = T )
}
