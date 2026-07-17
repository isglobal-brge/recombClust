## Legacy reference implementations (recombClust 1.0.0; Ruiz-Arenas et al. 2020).
##
## These are the original pure-R core functions, kept for INTERNAL use only: they
## are the golden reference against which the current C++ core is validated, and a
## fallback should a future check need the exact original behaviour. They are NOT
## exported and NOT part of the user API -- the going-forward functions are the
## unprefixed ones (recombFreq(), linkageFreq(), updateModel(), ...). Access these
## informally with recombClust:::recombFreq_legacy(), etc.

## recombFreq (recombClust 1.0.0): reestimate recomb-population block frequencies.
recombFreq_legacy <- function(Resp, Block) {
  nSNP <- (nchar(Block[1]) - 1) / 2
  freqs <- tapply(Resp, Block, sum)
  combs <- expand.grid(lapply(1:nSNP, function(x) c(0, 1)))
  leftLevs <- rightLevs <- apply(combs, 1, paste, collapse = "")
  freqsLeft <- vapply(leftLevs, function(x)
    sum(freqs[grep(paste0(x, "+"), names(freqs), fixed = TRUE)]), numeric(1))
  freqsRigth <- vapply(rightLevs, function(x)
    sum(freqs[grep(paste0("+", x), names(freqs), fixed = TRUE)]), numeric(1))
  freqsLeft <- sort(freqsLeft, decreasing = TRUE)
  freqsRigth <- sort(freqsRigth, decreasing = TRUE)
  leftLevs <- names(freqsLeft); nLevsLeft <- length(leftLevs)
  rightLevs <- names(freqsRigth); nLevsRigth <- length(rightLevs)
  AA <- diag(length(c(leftLevs, rightLevs)))
  AA[1:nLevsLeft, 1] <- -freqsLeft / max(freqsLeft)
  AA[(nLevsLeft + 1):nrow(AA), nLevsLeft + 1] <- -freqsRigth / max(freqsRigth)
  AA[1, 1:nLevsLeft] <- 1
  AA[nLevsLeft + 1, (nLevsLeft + 1):ncol(AA)] <- 1
  bb <- rep(0, nrow(AA)); bb[1] <- bb[nLevsLeft + 1] <- 1
  props <- qr.solve(AA, bb, tol = 1e-10)
  names(props) <- c(leftLevs, rightLevs)
  vec1 <- rep(leftLevs, nLevsRigth)
  vec2 <- rep(rightLevs, each = nLevsLeft)
  ans <- props[1:nLevsLeft] %*% t(props[(nLevsLeft + 1):length(props)])
  ans <- as.numeric(ans)
  names(ans) <- paste(vec1, vec2, sep = "+")
  ans
}

## linkageFreq (recombClust 1.0.0): reestimate linkage-population block frequencies.
linkageFreq_legacy <- function(Resp, Block) {
  nSNP <- (nchar(Block[1]) - 1) / 2
  freqs <- tapply(Resp, Block, sum)
  combs <- expand.grid(lapply(1:nSNP, function(x) c(0, 1)))
  leftLevs <- rightLevs <- apply(combs, 1, paste, collapse = "")
  sfreqs <- sort(freqs, decreasing = TRUE)
  leftLevsvec <- rightLevsvec <- leftLevs
  selAlleles <- character()
  while (length(sfreqs) > 0 & length(leftLevsvec) > 0) {
    allele <- names(sfreqs)[1]
    leftp <- strsplit(allele, "+", fixed = TRUE)[[1]][1]
    rigthp <- strsplit(allele, "+", fixed = TRUE)[[1]][2]
    if (leftp %in% leftLevsvec & rigthp %in% rightLevsvec) {
      selAlleles <- c(selAlleles, allele)
      leftLevsvec <- leftLevsvec[leftLevsvec != leftp]
      rightLevsvec <- rightLevsvec[rightLevsvec != rigthp]
    }
    sfreqs <- sfreqs[-1]
  }
  if (length(leftLevsvec) > 0)
    selAlleles <- c(selAlleles, paste(leftLevsvec, rightLevsvec, sep = "+"))
  props <- freqs[selAlleles]
  props[is.na(props)] <- 0
  names(props)[is.na(names(props))] <- selAlleles[!selAlleles %in% names(props)]
  props <- sort(props, decreasing = TRUE)
  AA <- diag(2^nSNP)
  AA[, 1] <- -props / max(props)
  AA[1, ] <- 1
  bb <- rep(0, nrow(AA)); bb[1] <- 1
  res <- qr.solve(AA, bb, tol = 1e-10)
  names(res) <- names(props)
  namesVec <- expand.grid(leftLevs, rightLevs)
  namesVec <- paste(namesVec[, 1], namesVec[, 2], sep = "+")
  res <- c(res, rep(0, length(namesVec) - length(props)))
  res[res < 1e-5] <- 1e-5
  res <- res / sum(res)
  names(res)[(length(props) + 1):length(res)] <-
    namesVec[!namesVec %in% names(res)]
  res
}

## updateModel (recombClust 1.0.0): one EM step (frequency functions passed in).
updateModel_legacy <- function(params, funcProps1, funcProps2) {
  newparams <- params
  newparams$r1 <- params$prob0 * params$props1[params$inds]
  newparams$r2 <- (1 - params$prob0) * params$props2[params$inds]
  R1 <- newparams$r1 / (newparams$r1 + newparams$r2)
  R2 <- newparams$r2 / (newparams$r1 + newparams$r2)
  newparams$prob0 <- mean(R1)
  newparams$props1 <- funcProps1(R1, newparams$inds)
  newparams$props2 <- funcProps2(R2, newparams$inds)
  newparams
}
