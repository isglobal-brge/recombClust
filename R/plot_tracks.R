#' Genomic-track view of the recombination signal
#'
#' Draws the recombClust result along the genome as stacked panels sharing a
#' base-pair axis: the per-window **recombination proportion** (overall, or one line
#' per chromosome cluster), the **difference between clusters** as a mountain plot,
#' and the detected regions as a strip. This complements the coordinate-free
#' correlation heat map of \code{\link{plot.recombRegions}} with a view that lines
#' up with genomic position.
#'
#' The recombination proportion of a window is the average recombination
#' responsibility of the chromosomes in that window (the \code{ModelsProb} column
#' mean); with \code{class} it is computed within each cluster, so their divergence
#' localises where the recombination histories differ. The panels are drawn with
#' base graphics, so no Bioconductor plotting package is required; for a fully
#' annotated genomic figure (ideogram, gene models) the same per-window values can
#' be fed to \pkg{Gviz} \code{DataTrack}s using real coordinates.
#'
#' @param probs A \code{"recombProbs"} object (from \code{\link{recomb_probabilities}}).
#' @param class Optional chromosome cluster labels (e.g. \code{recombClust(...)$class})
#'   in the same order as the chromosomes of \code{ModelsProb}; adds the per-cluster
#'   proportion and their difference.
#' @param regions Optional \code{"recombRegions"} object to mark as a region strip.
#' @param col Optional colours, one per cluster.
#' @param ... Ignored.
#' @return \code{probs}, invisibly (the figure is drawn as a side effect).
#' @export
plotRecombTracks <- function(probs, class = NULL, regions = NULL, col = NULL, ...) {
  if (!inherits(probs, "recombProbs"))
    stop("`probs` must be a 'recombProbs' object (see recomb_probabilities()).")

  MP <- as.matrix(BigDataStatMeth::hdf5_matrix(probs$filename, probs$dataset))
  BigDataStatMeth::hdf5_close_all()
  win <- probs$windows$windows
  mid <- (win[, 1] + win[, 2]) / 2
  xr <- range(win)

  if (is.null(class)) {
    props <- matrix(colMeans(MP), ncol = 1); cls <- "all"
  } else {
    if (length(class) != nrow(MP))
      stop("length(class) must equal the number of chromosomes in ModelsProb.")
    cls <- sort(unique(class))
    props <- vapply(cls, function(k) colMeans(MP[class == k, , drop = FALSE]),
                    numeric(ncol(MP)))
  }
  pal <- if (is.null(col)) c("#0072B2", "#D55E00", "#009E73", "#CC79A7") else col
  do_diff <- length(cls) >= 2
  has_reg <- !is.null(regions) && nrow(regions) > 0

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  npanel <- 1 + do_diff
  graphics::layout(matrix(seq_len(npanel + 1), ncol = 1),
                   heights = c(rep(1, npanel), 0.4))

  ## Panel 1: recombination proportion per window (one line per cluster)
  graphics::par(mar = c(0.4, 5, 2, 1))
  plot(NA, xlim = xr, ylim = c(0, max(props, na.rm = TRUE) * 1.05),
       xaxt = "n", xlab = "", ylab = "Recomb.\nproportion",
       main = "recombClust: recombination signal along the region")
  for (j in seq_along(cls))
    graphics::lines(mid, props[, j], col = pal[j], lwd = 2)
  graphics::legend("topleft", legend = if (identical(cls, "all")) "all" else
                   paste("cluster", cls), col = pal[seq_along(cls)], lwd = 2, bty = "n")

  ## Panel 2: difference between the first two clusters, as a mountain
  if (do_diff) {
    d <- props[, 1] - props[, 2]
    graphics::par(mar = c(0.4, 5, 0.6, 1))
    plot(NA, xlim = xr, ylim = range(0, d) * 1.1, xaxt = "n", xlab = "",
         ylab = "Difference\n(clusters 1-2)")
    graphics::abline(h = 0, col = "grey60")
    graphics::polygon(c(mid[1], mid, mid[length(mid)]), c(0, pmax(d, 0), 0),
                      col = "#F59E42", border = NA)
    graphics::polygon(c(mid[1], mid, mid[length(mid)]), c(0, pmin(d, 0), 0),
                      col = "#56B4E9", border = NA)
    graphics::lines(mid, d, col = "#7F2704", lwd = 1.5)
  }

  ## Bottom strip: detected regions + genomic axis
  graphics::par(mar = c(3.2, 5, 0.6, 1))
  plot(NA, xlim = xr, ylim = c(0, 1), yaxt = "n", xlab = "genomic position (bp)",
       ylab = "regions")
  if (has_reg)
    graphics::rect(regions$start, 0.3, regions$end, 0.7,
                   col = grDevices::adjustcolor("#D55E00", 0.6), border = "#7F2704")
  invisible(probs)
}
