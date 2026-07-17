// recombclust_pairs.cpp — Rcpp entry point for SNP-block pairing.

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "recombclust_pairs.h"

using namespace Rcpp;

//' Enumerate SNP-block pairs for the LD mixture model
//'
//' Forms blocks of \code{blockSize} consecutive SNPs and returns the pairs of
//' blocks that are non-overlapping and closer than \code{distance} (block j's
//' start within \code{distance} of block i's start). Replaces the
//' GenomicRanges/findOverlaps logic of the reference pipeline.
//'
//' @param positions Numeric vector of SNP positions (bp), sorted ascending.
//' @param blockSize Integer, number of SNPs per block (default 2).
//' @param distance Numeric, maximum bp between paired blocks' starts (default 1e4).
//' @return Integer matrix with two columns (\code{left}, \code{right}), one row per
//'   pair, holding 1-based block indices (i < j).
//' @keywords internal
//' @export
// [[Rcpp::export]]
IntegerMatrix block_pairs(NumericVector positions, int blockSize = 2,
                             double distance = 1e4) {
    std::vector<double> pos(positions.begin(), positions.end());
    std::vector<std::pair<int, int>> pairs =
        recombclust::block_pairs(pos, blockSize, distance);

    IntegerMatrix out(pairs.size(), 2);
    for (std::size_t k = 0; k < pairs.size(); ++k) {
        out(k, 0) = pairs[k].first + 1;    // 1-based block index (left)
        out(k, 1) = pairs[k].second + 1;   // 1-based block index (right)
    }
    colnames(out) = CharacterVector::create("left", "right");
    return out;
}
