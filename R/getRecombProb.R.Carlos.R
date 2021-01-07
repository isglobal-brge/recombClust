## Functions to create probTab (matrix of probabilities per windows)

#' recombRes: resultado de CrunLDmixtureModel
#' range: GRanges de la región
createMatrix <- function(recombRes, range){
  
  mat <- do.call(cbind, lapply(recombRes, `[[`, "r1"))
  allCoords <- data.frame(t(sapply(recombRes, function(x) x$annot)))
  allCoords$chr <- seqnames(range)
  grMods <- makeGRangesFromDataFrame(allCoords)
  list(mat = mat, annot = grMods)
}

#' mat, annot: elementos de lista de createMatrix
#' range: ídem a range de createMatrix
#' window: tamaño de la ventana (500bp por defecto)
getRecombProb1 <- function(mat, annot, range, window){
  
  ## Define chunks
  starts <- seq(start(range), end(range), window)
  chunks <- GRanges(seqnames = seqnames(range), IRanges(start = starts, width = window))
  
  overLaps <- findOverlaps(chunks, annot, type = "within")
  
  ## Compute cluster Recomb freq by mean of voting
  res <- sapply(seq_len(length(chunks)), function(chunk){
    
    sel <- to(overLaps)[from(overLaps) == chunk]
    
    ## --- Des d'aquí
    ##          Aquí afegir la funció que retorni el val ¿? 
    
    if (length(sel) == 0){
      return(rep(NA, nrow(mat)))
    }
    
    col <- mat[, sel, drop = FALSE]
    val <- rowMeans(col)
    
    ## --- Fins aquí
  
  })
  colnames(res) <- paste(start(chunks), end(chunks), sep = "-")
  res
}



getRecombProb2 <- function(mat, annot, range, window){
  
  ## Define chunks
  starts <- seq(start(range), end(range), window)
  chunks <- GRanges(seqnames = seqnames(range), IRanges(start = starts, width = window))
  
  overLaps <- findOverlaps(chunks, annot, type = "within")
  
  ## Compute cluster Recomb freq by mean of voting
  res <- sapply(seq_len(length(chunks)), function(chunk){
    
    sel <- to(overLaps)[from(overLaps) == chunk]
    if (length(sel) == 0){
      return(rep(NA, nrow(mat)))
    }
    val <-  getProbs(mat, sel)

  })
  
  rownames(res) <- rownames(mat)
  colnames(res) <- paste(start(chunks), end(chunks), sep = "-")
  res
}
