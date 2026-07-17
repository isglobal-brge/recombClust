#' Print a recombClust result
#'
#' @param x A \code{"recombClust"} object from \code{\link{recombClust}}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.recombClust <- function(x, ...) {
  cat("recombClust - chromosome clustering by recombination history\n")
  cat(sprintf("  chromosomes : %d\n", length(x$class)))
  cat(sprintf("  models      : %d (blockSize=%d, distance=%g)\n",
              x$models$ncol, x$params$blockSize, x$params$distance))
  cat(sprintf("  PCs         : %d  |  clusters: %d\n",
              x$params$PCs, x$params$clusters))
  sizes <- table(x$class)
  cat(sprintf("  cluster sizes: %s\n",
              paste(sprintf("%s=%d", names(sizes), as.integer(sizes)),
                    collapse = ", ")))
  cat(sprintf("  sdev[1:%d]  : %s\n", min(3, length(x$sdev)),
              paste(round(utils::head(x$sdev, 3), 4), collapse = ", ")))
  invisible(x)
}

#' Summarize a recombClust result
#'
#' @param object A \code{"recombClust"} object.
#' @param ... Ignored.
#' @return A list with the cluster sizes and the variance explained by each PC.
#' @export
summary.recombClust <- function(object, ...) {
  var <- object$sdev^2
  list(cluster_sizes = table(object$class),
       prop_variance = utils::head(var / sum(var), object$params$PCs))
}

## Colour-blind-friendly qualitative palette (Okabe-Ito), recycled if needed.
.rc_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                 "#E69F00", "#56B4E9", "#F0E442", "#999999")

#' Plot a recombClust result (PCA of recombination profiles)
#'
#' Scatter plot of the chromosomes on two principal components of their
#' recombination-model responsibilities, coloured by the assigned cluster. This is
#' the signature recombClust visualisation: chromosomes that share a recombination
#' history fall in the same cloud, and the clusters recover the latent
#' sub-populations (e.g. inversion orientations).
#'
#' @param x A \code{"recombClust"} object from \code{\link{recombClust}}.
#' @param comps Integer pair: which principal components to draw (default c(1, 2)).
#' @param col Optional vector of colours, one per cluster.
#' @param pch Plotting symbol (default 19).
#' @param legend Draw a cluster legend (default TRUE).
#' @param ... Passed to \code{plot}.
#' @return \code{x}, invisibly.
#' @export
plot.recombClust <- function(x, comps = c(1, 2), col = NULL, pch = 19,
                             legend = TRUE, ...) {
  sc <- x$scores
  cl <- as.integer(as.factor(x$class))
  k <- length(unique(cl))
  if (is.null(col)) col <- .rc_palette[((seq_len(k) - 1) %% length(.rc_palette)) + 1]
  ve <- x$sdev^2 / sum(x$sdev^2) * 100
  lab <- function(i) sprintf("PC%d (%.1f%%)", i, ve[i])
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::par(mar = c(4.2, 4.2, 2.5, 1))
  plot(sc[, comps[1]], sc[, comps[2]], col = grDevices::adjustcolor(col[cl], 0.75),
       pch = pch, xlab = lab(comps[1]), ylab = lab(comps[2]),
       main = "recombClust: chromosomes by recombination history", ...)
  if (legend)
    graphics::legend("topright", legend = paste("cluster", sort(unique(x$class))),
                     col = col, pch = pch, bty = "n")
  invisible(x)
}
