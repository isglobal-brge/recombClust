## hdf5omics imputation — fills missing entries of a 0/1 phased haplotype matrix,
## out-of-core, without a reference panel. Panel-based tools (Beagle, IMPUTE2) need a
## reference and assume 0/1/2 dosage; they do not apply here. Default is LD-kNNi
## (LinkImpute; Money et al. 2015), the best panel-free method in the literature; a
## per-SNP major-allele "mode" baseline is also available. See docs/DESIGN.md 6.3.

#' Impute missing genotypes of a 0/1 haplotype matrix (LD-kNNi or mode)
#'
#' Thin R wrapper over the C++ \code{hdf5omics_impute_hdf5}. Fills missing entries
#' (\code{NA}) of an HDF5 haplotype matrix and writes the completed matrix to a new
#' dataset. Imputation is reference-panel-free (in-sample), so it works on any HDF5
#' haplotype matrix without external data.
#'
#' \code{method = "ldknni"} (default) is LD-aware k-nearest-neighbour imputation
#' (LinkImpute): for each SNP with missing data it uses the \code{L} most-correlated
#' SNPs (linkage disequilibrium) and imputes each missing allele from the \code{k}
#' nearest haplotypes by inverse-distance-weighted majority vote. Because it uses LD,
#' it preserves the local haplotype structure that recombClust models, unlike the
#' \code{method = "mode"} baseline (per-SNP major allele), which biases allele
#' frequencies. Phased 1000 Genomes data are usually complete, so this is a safety
#' net for data that do contain missing calls.
#'
#' @param file HDF5 file.
#' @param in_group,in_dataset Source haplotype matrix (chromosomes x SNP, 0/1).
#' @param out_group,out_dataset Destination for the completed matrix.
#' @param method "ldknni" (default) or "mode".
#' @param k Number of nearest neighbours (LD-kNNi; default 5).
#' @param L Number of most-correlated SNPs used for the distance (LD-kNNi; default 20).
#' @param positions Optional SNP positions; returned unchanged (imputation keeps all
#'   SNPs) for convenience when chaining with the rest of the pipeline.
#' @return A list descriptor: \code{filename}, \code{group}, \code{dataset},
#'   \code{nrow}, \code{ncol}, \code{n_missing} (number of imputed entries),
#'   \code{method} and, if supplied, \code{positions}.
#' @export
hdf5omics_impute <- function(file, in_group = "haplotypes", in_dataset = "H",
                             out_group = in_group,
                             out_dataset = paste0(in_dataset, "_imp"),
                             method = c("ldknni", "mode"), k = 5L, L = 20L,
                             positions = NULL) {
  method <- match.arg(method)
  desc <- hdf5omics_impute_hdf5(file, in_group, in_dataset, out_group, out_dataset,
                                method = method, k = as.integer(k), L = as.integer(L))
  if (!is.null(positions)) desc$positions <- positions
  desc
}
