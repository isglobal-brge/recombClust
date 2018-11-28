#' Run LDmixture model to a pair of SNP-blocks
#'
#' @param dat Matrix with the genotype data
#' @param maxSteps Numerical with the maximum number of iterations run by the EM algorithm
#' @param prob0 Initial mixture probability.
#' @return A list with the LDmixture results
#' \itemize{
#'  \item{"logMix"}{Log-likelihood of mixture model}
#'  \item{"logLD"}{Log-likelihood of linkage model}
#'  \item{"logNoLD"}{Log-likelihood of recomb model}
#'  \item{"BIC"}{BIC of the mixture vs the base model}
#'  \item{"prob"}{Proportion of chromosomes belonging to recomb model}
#'  \item{"steps"}{Number of iterations until converge of the EM algorithm}
#'  \item{"pval"}{P-value of the Chi-square test}
#'  \item{"r1"}{Responsibilities for recomb population of each chromosomes. It is
#'  only available for selected models (BIC > 10, pval > 0.05)}
#' }
LDmixtureModel <- function(dat, maxSteps = 1000, prob0 = 0.5) {
  #reduce data to haplotype frequencies for block b1b2 and b3b4
  inds <- apply(dat, 1, paste, collapse = "+")
  r1 <- rep(1, length(inds))

  nSNP <- (nchar(inds[1]) - 1)/2

  propsRecomb <- recombFreq(r1, inds)
  propsLink <- linkageFreq(r1, inds)

  #compute likelihood of no inversion model first
  r1 <- propsRecomb[inds]
  LoglikeRecomb <-sum(log(r1))

  r2 <- propsLink[inds]
  LoglikeLinkage <- sum(log(r2))

  r1 <- prob0*propsRecomb[inds]
  r2<-(1-prob0)*propsLink[inds]


  params <- list(r1 = r1, r2 = r2, props1 = propsRecomb, props2 = propsLink,
                 prob0 = prob0, inds = inds)

  #EM loop
  #control while loop
  tol<- 1
  MINTOL<- .000000001
  steps<-1

  while(tol > MINTOL & steps <= maxSteps)
  {

    newparams <- updateModel(params, recombFreq, linkageFreq)

    tol<-sqrt((params$props1 - newparams$props1)%*%(params$props1 - newparams$props1) +
                (params$props2 - newparams$props2)%*%(params$props2 - newparams$props2) +
                abs(params$prob0 - newparams$prob0) )

    params <- newparams
    steps <- steps+1
  }

  #get last values to compute likelihood of the complete inversion model
  r1 <- params$prob0*params$props1[params$inds]
  r2 <- (1 - params$prob0)*params$props2[params$inds]
  LoglikeMix <- sum(log(r1 + r2))
  R1 <- r1/(r1 + r2)

  bicLD <- -2*LoglikeLinkage + log(length(params$inds))*(nSNP-1)
  bicNoLD <- -2*LoglikeRecomb + log(length(params$inds))*(nSNP-1)^2
  bicMix <- -2*LoglikeMix + log(length(params$inds))*((nSNP-1) + (nSNP-1)^2 + 1)
  bicDiff <- min(c(bicLD, bicNoLD)) - bicMix

  haplos <- names(params$props2)
  modelProps <- ((1 - params$prob0)*params$props2)[haplos] + (params$prob0*params$props1)[haplos]
  datProps <- table(inds)[names(modelProps)]
  datProps[is.na(datProps)] <- 0
  names(datProps) <- names(modelProps)

  chp.val <- stats::chisq.test(datProps, p = modelProps)$p.value


  ans <- list(logMix = LoglikeMix, logLD = LoglikeLinkage,
              logNoLD = LoglikeRecomb,
              bic = bicDiff, prob = params$prob0, steps = steps,
              pval = chp.val, r1 = R1)

  ans
}
