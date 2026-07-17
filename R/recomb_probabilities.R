#' Windowed recombination probabilities (Path B, out-of-core)
#'
#' Runs the region orchestrator to fit the LD mixture model to every SNP-block pair
#' and then summarises the per-chromosome recomb responsibilities into a windowed
#' probability matrix \code{ModelsProb} (chromosomes x windows): each window averages
#' the models whose genomic range spans it. All computation is on the BigDataStatMeth
#' C++/HDF5 backend; the result stays on disk (returned as a descriptor). This is the
#' Path B counterpart of \code{\link{recombClust}} (which clusters instead), used
#' downstream for region detection.
#'
#' @param filename HDF5 file holding the haplotype matrix; results are written here.
#' @param positions Numeric SNP positions (bp), sorted ascending.
#' @param hap_group,hap_dataset Location of the haplotype matrix (chromosomes x SNP).
#' @param out_group Group under which \code{models} and \code{ModelsProb} are written.
#' @param range Optional c(start, end) region to tile; defaults to the span of positions.
#' @param window Window width in bp (default 1000).
#' @param blockSize,distance,maxSteps,prob0 Mixture-model parameters.
#' @return An object of class \code{"recombProbs"}: a list with \code{models} and
#'   \code{windows} descriptors, the \code{filename}/\code{dataset} of ModelsProb, and
#'   \code{params}.
#' @export
recomb_probabilities <- function(filename, positions,
                                 hap_group = "haplotypes", hap_dataset = "H",
                                 out_group = "recombClust", range = NULL,
                                 window = 1000, blockSize = 2, distance = 1e4,
                                 maxSteps = 100, prob0 = 0.5) {
  desc <- run_region(filename, hap_group, hap_dataset, positions,
                        out_group, "models", blockSize, distance, maxSteps, prob0)
  BigDataStatMeth::hdf5_close_all()

  rs <- if (is.null(range)) min(positions) else range[1]
  re <- if (is.null(range)) max(positions) else range[2]
  w <- windows_hdf5(filename, out_group, "models", positions, desc$pairs,
                       rs, re, window, out_group, "ModelsProb")

  structure(list(models = desc, windows = w, filename = filename,
                 dataset = paste0(out_group, "/ModelsProb"),
                 params = list(window = window, blockSize = blockSize,
                               distance = distance)),
            class = "recombProbs")
}

#' Print a recombProbs result
#' @param x A \code{"recombProbs"} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.recombProbs <- function(x, ...) {
  cat("recombProbs - windowed recombination probabilities (on disk)\n")
  cat(sprintf("  file    : %s\n", x$filename))
  cat(sprintf("  dataset : %s  (%d chromosomes x %d windows)\n",
              x$dataset, x$windows$nrow, x$windows$ncol))
  cat(sprintf("  models  : %d (blockSize=%d, distance=%g)  window=%g bp\n",
              x$models$ncol, x$params$blockSize, x$params$distance, x$params$window))
  invisible(x)
}
