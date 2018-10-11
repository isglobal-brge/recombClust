<<<<<<< HEAD
#' Reestimate linkage block frequencies
#'
#' @param Resp Numerical with the responsibilities for linkage population
#' @param Block Character with the blocks genotypes
#' @return Numerical with the new block frequencies for linkage population
newFreqLD <-function(Resp, Block)
{
  nSNP <- (nchar(Block[1]) - 1)/2

  ## Compute blocks frequency
  freqs <-tapply(Resp, Block, sum)

  # Get all possible left and right levels
  ## Make all combination with 0 and 1
  combs <- expand.grid(lapply(1:nSNP, function(x) c(0, 1)))
  ## Bind to get levels
  leftLevs <- rightLevs <- apply(combs, 1, paste, collapse = "")

  ## Find combination containing more samples
  sfreqs <- sort(freqs, decreasing = TRUE)
  leftLevsvec <- rightLevsvec <- leftLevs
  selAlleles <- character()
  while(length(sfreqs) > 0 & length(leftLevsvec) > 0){
    allele <- names(sfreqs)[1]
    leftp <- strsplit(allele, "+", fixed = TRUE)[[1]][1]
    rigthp <- strsplit(allele, "+", fixed = TRUE)[[1]][2]
    if (leftp %in% leftLevsvec & rigthp %in% rightLevsvec){
      selAlleles <- c(selAlleles, allele)
      leftLevsvec <- leftLevsvec[leftLevsvec != leftp]
      rightLevsvec <- rightLevsvec[rightLevsvec != rigthp]
    }
    sfreqs <- sfreqs[-1]
  }

  if (length(leftLevsvec) > 0){
    selAlleles <- c(selAlleles, paste(leftLevsvec,rightLevsvec, sep = "+" ))
  }

  props <- freqs[selAlleles]
  props[is.na(props)] <- 0
  names(props)[is.na(names(props))] <- selAlleles[!selAlleles %in% names(props)]
  props <- sort(props, decreasing = TRUE)

  #solve linear algebra to find new frequency values
  AA <- diag(2^nSNP)
  AA[,1] <- -props/max(props)
  AA[1, ] <- 1

  bb <- rep(0, nrow(AA))
  bb[1] <- 1

  res <- qr.solve(AA, bb, tol = 1e-10)
  names(res) <- names(props)

  ## Get list of all combinations
  namesVec <- expand.grid(leftLevs, rightLevs)
  namesVec <- paste(namesVec[, 1], namesVec[, 2], sep = "+")

  res <- c(res, rep(0, length(namesVec) - length(props)))
  res[res < 1e-5] <- 1e-5
  res <- res/sum(res)
  names(res)[(length(props) +1):length(res)] <-
    namesVec[!namesVec %in% names(res)]
  return(res)
}
