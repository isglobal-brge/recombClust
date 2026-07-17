#' Detect recombination regions (BigLD) from windowed probabilities
#'
#' Region-detection step of the Path B pipeline (v2.0). From a windowed
#' recombination-probability matrix \code{ModelsProb} (chromosomes x windows,
#' on disk) it (1) computes the windows x windows correlation with
#' BigDataStatMeth's out-of-core C++ correlation (dropping windows with too much
#' missing data), and (2) partitions the windows into LD blocks with the C++
#' BigLD/CLQD port. All HDF5 I/O and the correlation run on the BigDataStatMeth
#' C++ API; only the graph primitives inside BigLD use igraph (see
#' \code{recombclust_bigld.hpp}, the single documented exception).
#'
#' The partition reproduces \code{gpart::CLQD} / \code{gpart::BigLD} exactly.
#' Note the correlation threshold acts on \code{|r|} (as in gpart's LD="r2"
#' branch, which thresholds \code{abs(cor)}, not \code{cor^2}).
#'
#' @param x A \code{"recombProbs"} object (from \code{\link{recomb_probabilities}}).
#' @param CLQcut Correlation threshold |r| for two windows to be linked (0-1).
#' @param clstgap Physical-distance gap (bp) that splits a clique.
#' @param minLength Minimum number of windows in a reported block.
#' @param hrstType CLQD heuristic: "fast" (default), "near-nonhrst" or "nonhrst".
#' @param hrstParam Heuristic core threshold.
#' @param CLQmode "maximal" (default) or "density".
#' @param chrN Chromosome label for the returned blocks.
#' @param group HDF5 group holding \code{ModelsProb} (default "recombClust").
#' @return An object of class \code{"recombRegions"}: a data.frame of LD blocks
#'   (chr, start, end, start.idx, end.idx) with the call parameters attached.
#' @export
define_Regions <- function(x, CLQcut = 0.5, clstgap = 40000, minLength = 1L,
                           hrstType = c("fast", "near-nonhrst", "nonhrst"),
                           hrstParam = 200L, CLQmode = c("maximal", "density"),
                           chrN = "chr1", group = "recombClust") {
  hrstType <- match.arg(hrstType)
  CLQmode <- match.arg(CLQmode)
  if (!inherits(x, "recombProbs"))
    stop("`x` must be a 'recombProbs' object (see recomb_probabilities()).")

  filename <- x$filename
  win <- x$windows$windows                 # ModelsProb windows: (start, end) per column
  if (is.null(win) || ncol(win) != 2L)
    stop("could not find ModelsProb window ranges in `x$windows$windows`.")

  ## (1) windows x windows |correlation| via BigDataStatMeth (out-of-core).
  desc <- cor_hdf5(filename, group, "ModelsProb", group, "Corr",
                      absval = TRUE)
  BigDataStatMeth::hdf5_close_all()

  kept <- desc$kept                        # 1-based windows surviving the NA filter
  ws <- win[kept, 1]
  we <- win[kept, 2]

  ## (2) BigLD/CLQD partition into LD blocks (C++ port; reads Corr via BDSM).
  blocks <- bigld_hdf5(filename, group, "Corr", ws, we,
                          CLQcut = CLQcut, clstgap = clstgap, hrstType = hrstType,
                          hrstParam = hrstParam, CLQmode = CLQmode, chrN = chrN)

  if (nrow(blocks) > 0L && minLength > 1L) {
    len <- blocks$end.idx - blocks$start.idx + 1L
    blocks <- blocks[len >= minLength, , drop = FALSE]
  }

  structure(blocks,
            class = c("recombRegions", "data.frame"),
            params = list(CLQcut = CLQcut, clstgap = clstgap, minLength = minLength,
                          hrstType = hrstType, CLQmode = CLQmode, nwin = desc$nwin,
                          filename = filename, corr_group = group, corr_dataset = "Corr",
                          win_start = ws, win_end = we))
}

#' Print recombRegions
#' @param x A \code{"recombRegions"} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.recombRegions <- function(x, ...) {
  p <- attr(x, "params")
  cat("recombRegions - LD blocks from BigLD (on-disk correlation)\n")
  cat(sprintf("  windows : %d (after NA filter)  CLQcut=%.2f  clstgap=%g\n",
              p$nwin, p$CLQcut, p$clstgap))
  cat(sprintf("  blocks  : %d  (hrstType=%s, CLQmode=%s)\n",
              nrow(x), p$hrstType, p$CLQmode))
  print(utils::head(as.data.frame(x)))
  invisible(x)
}

#' Plot detected recombination regions
#'
#' Draws the windows x windows absolute-correlation matrix of the recombination
#' profiles as a heat map with the detected LD blocks outlined (the recombClust
#' analogue of an LD heat map): strong off-diagonal blocks are groups of windows
#' that share a recombination history, and each outlined square is a reported
#' region. If the correlation matrix is no longer on disk, a genomic block diagram
#' (blocks along the base-pair axis) is drawn instead.
#'
#' @param x A \code{"recombRegions"} object from \code{\link{define_Regions}}.
#' @param ... Passed to \code{image}/\code{plot}.
#' @return \code{x}, invisibly.
#' @export
plot.recombRegions <- function(x, ...) {
  p <- attr(x, "params")
  C <- tryCatch({
    hm <- BigDataStatMeth::hdf5_matrix(p$filename, paste0(p$corr_group, "/", p$corr_dataset))
    m <- as.matrix(hm); hm$close(); BigDataStatMeth::hdf5_close_all(); m
  }, error = function(e) NULL)

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))

  if (!is.null(C)) {
    n <- nrow(C)
    pal <- grDevices::colorRampPalette(
      c("#FFFFFF", "#FEE0A0", "#F59E42", "#D55E00", "#7F2704"))(100)
    graphics::par(mar = c(4.2, 4.2, 2.5, 1))
    graphics::image(seq_len(n), seq_len(n), C[, rev(seq_len(n))], col = pal,
                    zlim = c(0, 1), xlab = "window", ylab = "window",
                    main = "recombClust regions: |correlation| with LD blocks", ...)
    if (nrow(x) > 0)
      for (i in seq_len(nrow(x))) {
        s <- x$start.idx[i]; e <- x$end.idx[i]
        graphics::rect(s - 0.5, n - e + 0.5, e + 0.5, n - s + 0.5,
                       border = "#111111", lwd = 2)
      }
  } else {
    ## fallback: blocks along the genomic (bp) axis
    xr <- range(c(p$win_start, p$win_end, x$start, x$end))
    graphics::par(mar = c(4.2, 2, 2.5, 1))
    plot(NA, xlim = xr, ylim = c(0, 1), yaxt = "n", ylab = "",
         xlab = "genomic position (bp)",
         main = "recombClust regions (LD blocks)", ...)
    if (nrow(x) > 0)
      graphics::rect(x$start, 0.35, x$end, 0.65,
                     col = grDevices::adjustcolor("#D55E00", 0.6), border = "#7F2704")
  }
  invisible(x)
}
