#' Cluster chromosomes by recombination history (out-of-core, Path A)
#'
#' Runs the recombClust pipeline on a haplotype matrix stored in an HDF5 file: the
#' C++ region orchestrator fits the LD mixture model to every SNP-block pair and
#' streams the per-chromosome recomb responsibilities into an on-disk models matrix;
#' a principal component analysis of that matrix is then computed out-of-core with
#' BigDataStatMeth, and the chromosomes are grouped with k-means on the leading PCs.
#'
#' The heavy computation runs on the BigDataStatMeth C++/HDF5 backend; only the small
#' final k-means (chromosomes x PCs) runs in R.
#'
#' @param filename Path to the HDF5 file holding the haplotype matrix. The models
#'   matrix and PCA are written into the same file.
#' @param positions Numeric vector of SNP positions (bp), sorted ascending, one per
#'   column of the haplotype matrix.
#' @param hap_group,hap_dataset Location of the haplotype matrix (chromosomes x SNP,
#'   0/1) inside \code{filename}.
#' @param out_group Group under which the models matrix is written.
#' @param clusters Number of k-means clusters (default 2).
#' @param PCs Number of principal components used for clustering (default 3).
#' @param blockSize SNPs per block (default 2).
#' @param distance Maximum bp between paired blocks' starts (default 1e4).
#' @param maxSteps,prob0 EM parameters passed to the mixture model.
#' @param seed Optional integer seed for reproducible k-means.
#' @return An object of class \code{"recombClust"}: a list with \code{class}
#'   (chromosome cluster assignments), \code{sdev} (PCA standard deviations),
#'   \code{scores} (chromosomes x PCs), \code{models} (descriptor of the on-disk
#'   models matrix) and \code{params}.
#' @export
recombClust <- function(filename, positions,
                            hap_group = "haplotypes", hap_dataset = "H",
                            out_group = "recombClust",
                            clusters = 2, PCs = 3,
                            blockSize = 2, distance = 1e4,
                            maxSteps = 100, prob0 = 0.5, seed = NULL) {

  ## 1. Region orchestrator (C++): block pairing + per-pair EM + streamed models.
  desc <- run_region(filename, hap_group, hap_dataset, positions,
                        out_group, "models", blockSize, distance, maxSteps, prob0)

  ## 2. PCA of the models matrix via BigDataStatMeth's C++ SVD (out-of-core; not the
  ##    R prcomp S3 method — DESIGN 4.6). Returns scores (chromosomes x PCs) and sdev.
  BigDataStatMeth::hdf5_close_all()
  pc <- pca_hdf5(filename, out_group, "models", PCs)
  scores <- pc$scores

  ## 3. k-means on the leading PCs (small; stays in R — the one compute step in R).
  if (!is.null(seed)) set.seed(seed)
  cl <- stats::kmeans(scores, centers = clusters, nstart = 1000)$cluster

  structure(list(class = cl, sdev = pc$sdev, scores = scores, models = desc,
                 params = list(clusters = clusters, PCs = PCs,
                               blockSize = blockSize, distance = distance)),
            class = "recombClust")
}
