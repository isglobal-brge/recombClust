#' Run LDmixtureModel to a dataset
#'
#' This function makes SNP-blocks and pair them. Then, it applies LDmixtureModel
#' to all SNP-block pairs closer than a selected distance.
#'
#' @export
#'
#' @param haplos Matrix with the haplotypes (SNPs in rows, samples in columns)
#' @param annot GenomicRanges with the SNPs annotation
#' @param blockSize Numeric with the size of the SNP block (Default: 2)
#' @param distance Numeric with the maximum distance in bases to pair two blocks. (Default: 1e5)
#' @param BPPARAM  An object from \code{BiocParallelParam}. It allow running the models
#' in parallel. By default, the models are run in serial.
#' @return A list with the results of the LDmixture models. Each element is a model
#' and contains the following items:
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
runLDmixtureModel <- function(haplos, annot, blockSize = 2, distance = 1e5, BPPARAM = BiocParallel::SerialParam(1)){

  # Make list of SNP pairs to test
  GRblocks <- GenomicRanges::GRanges(
    paste0(GenomicRanges::seqnames(annot)[1], ":",
           GenomicRanges::start(annot)[1:(length(annot) - blockSize + 1)], "-",
           GenomicRanges::end(annot)[blockSize:(length(annot))]))
  GRblocks$Ind <- lapply(1:(length(annot) - blockSize + 1),
                         function(x) seq(x, x + blockSize -1, 1))

  # Select pairs of blocks that are closer than distance
  overlaps <- GenomicRanges::findOverlaps(GenomicRanges::resize(GRblocks, distance), GRblocks)
  autoOverlaps <- GenomicRanges::findOverlaps(GRblocks)

  ## Remove overlapping pairs of blocks
  overlaps <- overlaps[!S4Vectors::`%in%`(overlaps, autoOverlaps)]

  # Run model over the SNP-block pairs
  models <- BiocParallel::bplapply(seq_len(length(overlaps)), function(ind){
    bl1 <- S4Vectors::from(overlaps)[ind]
    bl2 <- S4Vectors::to(overlaps)[ind]
    ind1 <- GRblocks$Ind[[bl1]]
    ind2 <- GRblocks$Ind[[bl2]]
    res <- inversionModel(cbind(
      apply(haplos[, ind1], 1, function(x) paste(x, collapse = "")),
      apply(haplos[, ind2], 1, function(x) paste(x, collapse = ""))
    ))
    res$annot <- c(start = GenomicRanges::start(GRblocks[bl1]),
                   end = GenomicRanges::start(GRblocks[bl2]))
    if (res$bic < 10 | res$pval < 0.05) {
      res$r1 <- NULL
    }
    res
  }, BPPARAM = BPPARAM)
  models
}
