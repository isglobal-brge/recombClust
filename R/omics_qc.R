## hdf5omics QC — operates directly on an HDF5 haplotype matrix, out-of-core, so it
## works whether the data came from a VCF, a text file, or was already in HDF5 (no
## SeqArray required). Correct 0/1 haplotype MAF (BigDataStatMeth's filter_maf assumes
## 0/1/2 dosage coding; see docs/DESIGN.md 6.3). Reuses BDSM's out-of-core colMeans.

#' Filter SNPs of an HDF5 haplotype matrix by minor allele frequency
#'
#' Thin R wrapper over the C++ \code{hdf5omics_filter_maf_hdf5}: all computation (the
#' out-of-core column means via BigDataStatMeth's C++ \code{get_HDF5_colMeans}, the
#' phased-haplotype MAF \code{min(freq1, 1 - freq1)}, and the column subset) happens on
#' the BigDataStatMeth C++ API. This R layer only forwards arguments and, if positions
#' are supplied, keeps them in sync with the surviving SNPs.
#'
#' Works on any HDF5 matrix (no VCF/SeqArray needed) and never materialises the full
#' matrix in memory. Unlike \code{BigDataStatMeth::filter_maf} (0/1/2 dosage coding),
#' the MAF here is the phased-haplotype definition; see docs/DESIGN.md 6.3.
#'
#' @param file HDF5 file.
#' @param in_group,in_dataset Source haplotype matrix.
#' @param out_group,out_dataset Destination for the filtered matrix.
#' @param maf_threshold Minimum minor allele frequency to keep a SNP (default 0.05).
#' @param positions Optional SNP positions; if given, the returned descriptor carries
#'   the positions of the surviving SNPs (kept in sync with the filtered columns).
#' @return A list descriptor: \code{filename}, \code{group}, \code{dataset},
#'   \code{nrow}, \code{ncol}, \code{keep} (1-based kept column indices), \code{maf}
#'   (MAF of kept SNPs) and, if supplied, \code{positions}.
#' @export
hdf5omics_filter_maf <- function(file, in_group = "haplotypes", in_dataset = "H",
                                 out_group = in_group,
                                 out_dataset = paste0(in_dataset, "_maf"),
                                 maf_threshold = 0.05, positions = NULL) {
  desc <- hdf5omics_filter_maf_hdf5(file, in_group, in_dataset,
                                    out_group, out_dataset, maf_threshold)
  if (!is.null(positions)) desc$positions <- positions[desc$keep]
  desc
}
