#' Reestimate recomb block frequencies
#'
#' @param Resp Numerical with the responsibilities for recomb population
#' @param Block Character with the blocks genotypes
#' @return Numerical with the new block frequencies for recomb population
newFreqBase <- function(Resp, Block)
{
  nSNP <- (nchar(Block[1]) - 1)/2

  ## Compute blocks frequency
  freqs <- tapply(Resp, Block, sum)

  # Get all possible left and right levels
  ## Make all combination with 0 and 1
  combs <- expand.grid(lapply(1:nSNP, function(x) c(0, 1)))
  ## Bind to get levels
  leftLevs <- rightLevs <- apply(combs, 1, paste, collapse = "")

  ## Compute frequencies partial blocks
  freqsLeft <- vapply(leftLevs, function(x)
    sum(freqs[grep(paste0(x, "+"), names(freqs), fixed = TRUE)]), numeric(1))
  freqsRigth <- vapply(rightLevs, function(x)
    sum(freqs[grep(paste0("+", x), names(freqs), fixed = TRUE)]), numeric(1))

  ## Sort combinations by frequency
  freqsLeft <- sort(freqsLeft, decreasing = TRUE)
  freqsRigth <- sort(freqsRigth, decreasing = TRUE)

  leftLevs <- names(freqsLeft)
  nLevsLeft <- length(leftLevs)

  rightLevs <- names(freqsRigth)
  nLevsRigth <- length(rightLevs)

  #solve linear algebra to find new frequency values
  AA<-diag(length(c(leftLevs, rightLevs)))
  AA[1:nLevsLeft,1]<- -freqsLeft/max(freqsLeft)
  AA[(nLevsLeft+1):nrow(AA), nLevsLeft+1]<- -freqsRigth/max(freqsRigth)

  AA[1,1:nLevsLeft]<-1
  AA[nLevsLeft+1, (nLevsLeft+1):ncol(AA)] <- 1

  bb<-rep(0, nrow(AA))
  bb[1] <- bb[nLevsLeft+1] <- 1

  props <- qr.solve(AA, bb, tol = 1e-10)
  names(props)<- c(leftLevs, c(rightLevs))

  vec1 <- rep(leftLevs, nLevsRigth)
  vec2 <- rep(rightLevs, each = nLevsLeft)

  ans <- props[1:nLevsLeft] %*% t(props[(nLevsLeft+1):length(props)])
  ans <- as.numeric(ans)

  names(ans) <- paste(vec1, vec2, sep = "+")
  return(ans)
}
