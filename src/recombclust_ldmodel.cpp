// recombclust_ldmodel.cpp — Rcpp entry points for the pure-Eigen mixture-model core.
//
// Thin [[Rcpp::export]] wrappers over recombclust_core.hpp. These carry no HDF5 /
// BigDataStatMeth dependency; they exist to (a) expose the numerical core to R and
// (b) enable the N1/N2 validation against the pure-R v1.0.0 reference.

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "recombclust_core.h"

using namespace Rcpp;

//' Reestimate recomb (no-LD) block frequencies
//' @param Resp Numeric responsibilities for the recomb population.
//' @param Block Character vector of "left+right" haplotype tokens.
//' @param nSNP Number of SNPs per block.
//' @return Named numeric vector of joint recomb frequencies.
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericVector recombFreq(NumericVector Resp, StringVector Block, int nSNP) {
    return recombclust::CRecombFreq(Resp, Block, nSNP);
}

//' Reestimate linkage (full-LD) block frequencies
//' @inheritParams recombFreq
//' @return Named numeric vector of joint linkage frequencies.
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericVector linkageFreq(NumericVector Resp, StringVector Block, int nSNP) {
    return recombclust::CLinkageFreq(Resp, Block, nSNP);
}

//' One EM step of the LD mixture model
//' @param par List with props1, props2, prob0, inds, blocksize.
//' @return Updated list (r1, r2, props1, props2, prob0, inds, blocksize).
//' @keywords internal
//' @export
// [[Rcpp::export]]
List updateModel(List par) {
    return recombclust::CupdateModel(par);
}

//' In-memory LD mixture model for a pair of SNP-blocks
//' @param dat Character matrix, 2 columns (left/right allele strings), one row per chromosome.
//' @param maxSteps Maximum EM iterations (default 100).
//' @param prob0 Initial mixture probability (default 0.5).
//' @param blocksize SNPs per block (default 2).
//' @return list(logNoLD, prob, r1).
//' @keywords internal
//' @export
// [[Rcpp::export]]
List LDmixtureModel(RObject dat, Nullable<int> maxSteps = R_NilValue,
                           Nullable<double> prob0 = R_NilValue,
                           Nullable<int> blocksize = R_NilValue) {
    int imaxSteps = maxSteps.isNotNull() ? as<int>(maxSteps) : 100;
    double dprob0 = prob0.isNotNull() ? as<double>(prob0) : 0.5;
    int iblocksize = blocksize.isNotNull() ? as<int>(blocksize) : 2;
    return recombclust::LDmixtureModel_mem(dat, imaxSteps, dprob0, iblocksize);
}

//' Windowed reduction: mean over selected model columns
//' @param mat Chromosomes x models responsibility matrix.
//' @param sel 1-based column indices of the models within a window.
//' @return Per-chromosome mean recomb probability for the window.
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericVector getProbs(RObject mat, RObject sel) {
    return recombclust::getProbs(mat, sel);
}
