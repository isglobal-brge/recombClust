// recombclust_bigld.cpp — thin Rcpp wrappers for the in-memory BigLD/CLQD port
// (validation entry points). HDF5-backed variants live in recombclust_regions.cpp
// (they route I/O + correlation through BigDataStatMeth). These take a correlation
// matrix already in memory and delegate the graph work to igraph via the header.

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "recombclust_bigld.h"

using namespace Rcpp;

//' CLQD clique decomposition (in-memory; validation entry point)
//'
//' Port of \code{gpart::CLQD}. Takes a pre-threshold correlation matrix and
//' returns the per-column bin assignment (0 = unassigned). The partition matches
//' \code{gpart::CLQD} exactly (bin labels may be permuted).
//'
//' @param OCM Correlation matrix (windows x windows), pre-threshold.
//' @param SNPbps Numeric bp position per column (unique).
//' @param CLQcut Correlation threshold |r| for adjacency.
//' @param clstgap Physical-distance gap (bp) that splits a clique.
//' @param hrstType One of "near-nonhrst", "fast", "nonhrst".
//' @param hrstParam Heuristic core threshold.
//' @param CLQmode "density" or "maximal".
//' @return Integer vector of bin ids (0 = NA/unassigned).
//' @keywords internal
//' @export
// [[Rcpp::export]]
IntegerVector clqd(Eigen::MatrixXd OCM, NumericVector SNPbps,
                      double CLQcut = 0.5, double clstgap = 40000,
                      std::string hrstType = "fast", int hrstParam = 200,
                      std::string CLQmode = "maximal") {
    try {
        if (OCM.rows() != OCM.cols())
            throw std::runtime_error("OCM must be square");
        if (SNPbps.size() != OCM.cols())
            throw std::runtime_error("length(SNPbps) must equal ncol(OCM)");
        std::vector<double> bps(SNPbps.begin(), SNPbps.end());
        std::vector<int> bin = recombclust::CLQD(OCM, CLQcut, clstgap, hrstType,
                                                 hrstParam, CLQmode, bps);
        IntegerVector out(bin.size());
        for (std::size_t i = 0; i < bin.size(); ++i)
            out[i] = bin[i] == 0 ? NA_INTEGER : bin[i];
        return out;
    } catch (std::exception& ex) {
        Rf_error("clqd: %s", ex.what());
    } catch (...) {
        Rf_error("clqd: unknown C++ exception");
    }
    return IntegerVector();
}

//' constructLDblock (validation entry point)
//'
//' Port of \code{gpart:::constructLDblock}. Takes a list of integer index sets
//' and returns LD blocks as a matrix of (start, end) index ranges.
//'
//' @param clstlist List of integer vectors (index sets, values as-is).
//' @return Integer matrix (blocks x 2): start and end index of each block.
//' @keywords internal
//' @export
// [[Rcpp::export]]
IntegerMatrix constructldblock(List clstlist) {
    try {
        std::vector<std::vector<int>> cl;
        for (int i = 0; i < clstlist.size(); ++i) {
            IntegerVector v = clstlist[i];
            cl.emplace_back(v.begin(), v.end());
        }
        std::vector<std::array<int, 2>> blocks = recombclust::constructLDblock(cl);
        IntegerMatrix out(static_cast<int>(blocks.size()), 2);
        for (std::size_t i = 0; i < blocks.size(); ++i) {
            out(i, 0) = blocks[i][0];
            out(i, 1) = blocks[i][1];
        }
        return out;
    } catch (std::exception& ex) {
        Rf_error("constructldblock: %s", ex.what());
    } catch (...) {
        Rf_error("constructldblock: unknown C++ exception");
    }
    return IntegerMatrix();
}

//' BigLD region blocks (in-memory; validation entry point)
//'
//' CLQD -> constructLDblock, mirroring the recombClust BigLD wrapper. Returns a
//' data.frame of LD blocks with 1-based window indices and their bp span.
//'
//' @param OCM Correlation matrix (windows x windows), pre-threshold.
//' @param win_start,win_end bp span of each window (SNPbps for CLQD = win_start).
//' @param CLQcut,clstgap,hrstType,hrstParam,CLQmode CLQD parameters.
//' @param chrN Chromosome label.
//' @return data.frame(chr, start, end, start.idx, end.idx), sorted by start.idx.
//' @keywords internal
//' @export
// [[Rcpp::export]]
DataFrame bigld(Eigen::MatrixXd OCM, NumericVector win_start, NumericVector win_end,
                   double CLQcut = 0.5, double clstgap = 40000,
                   std::string hrstType = "fast", int hrstParam = 200,
                   std::string CLQmode = "maximal", std::string chrN = "chr1") {
    try {
        if (OCM.rows() != OCM.cols()) throw std::runtime_error("OCM must be square");
        if (win_start.size() != OCM.cols() || win_end.size() != OCM.cols())
            throw std::runtime_error("win_start/win_end length must equal ncol(OCM)");
        std::vector<double> ws(win_start.begin(), win_start.end());
        std::vector<std::array<int, 2>> blocks =
            recombclust::BigLD(OCM, CLQcut, clstgap, hrstType, hrstParam, CLQmode, ws);

        const int nb = static_cast<int>(blocks.size());
        CharacterVector chr(nb);
        NumericVector start(nb), end(nb);
        IntegerVector sidx(nb), eidx(nb);
        for (int i = 0; i < nb; ++i) {
            chr[i] = chrN;
            start[i] = win_start[blocks[i][0]];
            end[i] = win_end[blocks[i][1]];
            sidx[i] = blocks[i][0] + 1;    // 1-based for R
            eidx[i] = blocks[i][1] + 1;
        }
        return DataFrame::create(Named("chr") = chr, Named("start") = start,
                                 Named("end") = end, Named("start.idx") = sidx,
                                 Named("end.idx") = eidx,
                                 Named("stringsAsFactors") = false);
    } catch (std::exception& ex) {
        Rf_error("bigld: %s", ex.what());
    } catch (...) {
        Rf_error("bigld: unknown C++ exception");
    }
    return DataFrame::create();
}
