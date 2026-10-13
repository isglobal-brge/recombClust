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
#' @param compression gzip level (0-9) for the on-disk models matrix (the one
#'   dataset this function creates); \code{NULL} defers to BigDataStatMeth's own
#'   default (6). The package default is 0 (no compression): the responsibility
#'   matrix compresses to a measured ratio of ~0.96 (\code{docs/BENCHMARK.md}
#'   \S13.3), so deflating it costs CPU on the hot write path for essentially no
#'   space — with the rewritten EM core that cost would dominate the pipeline.
#'   Datasets that do compress well (the imported haplotypes) have their own
#'   \code{compression} argument in \code{\link{hdf5omics_import}} (default 1).
#' @param svd_k,svd_q Blocking of the hierarchical SVD used for the PCA: the models
#'   matrix is split into \code{svd_k^svd_q} blocks. \code{svd_q = NULL} (default)
#'   selects it automatically from the matrix dimensions and
#'   \code{svd_target_block}. Keep \code{svd_k = 2}: raising the branching factor
#'   costs both memory and accuracy (BENCHMARK \S10.5).
#' @param threads Number of OpenMP threads; \code{NULL} (default) uses all
#'   available cores (\code{OMP_NUM_THREADS} is respected). **Scope:** affects
#'   the two compute-heavy stages — the per-pair EM fits (each pair is an
#'   independent, deterministic fit, so the result is bit-identical for any
#'   thread count) and the PCA/SVD on the BigDataStatMeth backend. The final
#'   k-means runs in R, single-threaded. All published benchmark figures were
#'   measured with 1 thread.
#' @param svd_target_block Target number of matrix elements per SVD block
#'   (default 1e7). The automatic selection picks the *smallest* block count whose
#'   blocks fit this target, i.e. the most accurate configuration that stays within
#'   the memory budget. The target is an absolute, machine-independent constant on
#'   purpose: making it depend on detected RAM would make results differ between
#'   machines, since the trailing components' accuracy depends on the block count.
#' @section Accuracy under blocking:
#'   Above roughly 5.4e7 matrix elements BigDataStatMeth switches to an approximate
#'   hierarchical SVD, and the approximation grows with the block count. The leading
#'   component stays accurate (measured error <= 1.3\%) and recombClust's clustering
#'   is invariant to this (BENCHMARK \S10.5.1), but \code{sdev} and \code{scores}
#'   beyond the first component should not be interpreted quantitatively without
#'   checking. Set \code{svd_q = 1} to get the most accurate (and most
#'   memory-hungry) configuration.
#' @return An object of class \code{"recombClust"}: a list with \code{class}
#'   (chromosome cluster assignments), \code{sdev} (PCA standard deviations),
#'   \code{scores} (chromosomes x PCs), \code{models} (descriptor of the on-disk
#'   models matrix) and \code{params} (including the effective \code{svd_k},
#'   \code{svd_q}, \code{svd_M} and \code{compression}).
#' @export
recombClust <- function(filename, positions,
                            hap_group = "haplotypes", hap_dataset = "H",
                            out_group = "recombClust",
                            clusters = 2, PCs = 3,
                            blockSize = 2, distance = 1e4,
                            maxSteps = 100, prob0 = 0.5, seed = NULL,
                            compression = 0L, svd_k = 2L, svd_q = NULL,
                            svd_target_block = 1e7, threads = NULL) {

  ## 1. Region orchestrator (C++): block pairing + per-pair EM (parallel over
  ## pairs; identical result for any thread count) + batched models streaming.
  ## compression: gzip level (0-9) for the on-disk models matrix; NULL keeps
  ## BigDataStatMeth's default (6). Default 0: the matrix compresses to ~0.96.
  comp <- if (is.null(compression)) -1L else as.integer(compression)
  thr <- if (is.null(threads)) 0L else as.integer(threads)
  desc <- run_region(filename, hap_group, hap_dataset, positions,
                        out_group, "models", blockSize, distance, maxSteps, prob0,
                        comp, thr)

  ## 2. PCA of the models matrix via BigDataStatMeth's C++ SVD (out-of-core; not the
  ##    R prcomp S3 method — DESIGN 4.6). Returns scores (chromosomes x PCs) and sdev.
  ##    svd_k/svd_q set the hierarchical blocking (svd_k^svd_q blocks): more blocks
  ##    bound peak memory at the cost of accuracy in the trailing components.
  ## Auto-selection of the blocking: smallest block count whose blocks fit the
  ## target, i.e. the most accurate configuration within the memory budget. Bounded
  ## by BigDataStatMeth's requirement that svd_k^svd_q <= number of columns.
  if (is.null(svd_q)) svd_q <- .svd_auto_q(desc$nrow, desc$ncol, svd_k,
                                           svd_target_block)
  svd_k <- as.integer(svd_k); svd_q <- as.integer(svd_q)

  ## Same threads setting for the SVD (BDSM interprets NULL as all cores).
  thr_svd <- if (is.null(threads)) NULL else as.integer(threads)
  BigDataStatMeth::hdf5_close_all()
  pc <- pca_hdf5(filename, out_group, "models", PCs, svd_k, svd_q, thr_svd)
  scores <- pc$scores

  ## 3. k-means on the leading PCs (small; stays in R — the one compute step in R).
  if (!is.null(seed)) set.seed(seed)
  cl <- stats::kmeans(scores, centers = clusters, nstart = 1000)$cluster

  structure(list(class = cl, sdev = pc$sdev, scores = scores, models = desc,
                 params = list(clusters = clusters, PCs = PCs,
                               blockSize = blockSize, distance = distance,
                               svd_k = svd_k, svd_q = svd_q,
                               svd_M = svd_k^svd_q,
                               threads = if (is.null(threads)) "all-cores"
                                         else as.integer(threads),
                               compression = if (is.null(compression)) "bdsm-default"
                                             else compression)),
            class = "recombClust")
}


## Smallest q (>= 1) such that a block of the n x p matrix fits `target` elements,
## subject to BigDataStatMeth's constraint k^q <= p. q = 1 (two blocks) is both the
## default of the backend and the most accurate setting, so the search starts there
## and only grows when the block would exceed the target.
.svd_auto_q <- function(n, p, k = 2L, target = 1e7) {
  q <- 1L
  while (n * (p / k^q) > target && k^(q + 1L) <= p) q <- q + 1L
  q
}
