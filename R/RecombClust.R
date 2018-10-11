#' RecombClust: Group chromosomes by recombination history
#'
#' RecombClust infers recombination history based on a mixture model that, given a pair
#' of SNP-blocks, separates chromosomes in two populations, one with high Linkage Disequilibrium (LD)
#' and low recombination (linkage) and another with low LD and high recombination. The
#' method use the classification of several SNP-block pairs in a region to group chromosomes
#' in clusters with different recombination history. This package takes as input genotype phased data.
#'
#' @docType package
#' @name RecombClust
#'
#' @importFrom BiocParallel bplapply
#' @importFrom GenomicRanges seqnames start end findOverlaps resize
#' @importFrom gtools permutations
#' @importFrom S4Vectors from to
#' @importFrom stats chisq.test kmeans prcomp
NULL
