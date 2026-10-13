// recombclust_core.h — R-facing adapters over the pure C++ mixture-model core.
//
// Phase B of the EM-core restructuring: the arithmetic lives in
// recombclust_core_pure.h (no R types, thread-safe); this header keeps the
// historical function names and signatures so the exported wrappers
// (recombclust_ldmodel.cpp) and the validation harness (tests 01/02 against the
// pure-R v1.0.0 oracle) compile and pass UNCHANGED. Conversion R <-> pure
// happens only here, once per call — never inside the EM loop.
//
// The removed string/map helpers of the old implementation (getNamesbyBlocks,
// VectortoOrderedMap, getMatrixfromMap, ...) are in git history (tag
// v2.0.1-validada); their arithmetic survives, reproduced bit-for-bit by the
// pure core (see the contract note in recombclust_core_pure.h).

#ifndef RECOMBCLUST_CORE_HPP
#define RECOMBCLUST_CORE_HPP

#include <RcppEigen.h>
#include <string>
#include <vector>
#include "recombclust_core_pure.h"

namespace recombclust {

using namespace Rcpp;

// Concatenate two string vectors position by position with `par` in the middle.
// (Kept: the LDmixtureModel wrapper uses it to form the per-chromosome tokens.)
inline CharacterVector concatenate(StringVector x, StringVector y, std::string par) {
    CharacterVector res(x.size());
    for (int i = 0; i < x.size(); i++)
        res[i] = ((Rcpp::as<std::string>(x[i]) + par) + Rcpp::as<std::string>(y[i]));
    return res;
}

// --- boundary conversions ---------------------------------------------------

// Encode per-chromosome "left+right" tokens to lexicographic slots.
inline std::vector<int> encode_tokens(const StringVector& Block, int nSNP) {
    std::vector<int> codes(Block.size());
    for (int i = 0; i < Block.size(); i++)
        codes[i] = pure::encode_token(as<std::string>(Block[i]), nSNP);
    return codes;
}

// Rebuild the named NumericVector a FreqVec stands for.
inline NumericVector freqvec_to_R(const pure::FreqVec& f, int nSNP) {
    NumericVector v(f.val.size());
    StringVector names(f.val.size());
    for (std::size_t i = 0; i < f.val.size(); i++) {
        v[i] = f.val[i];
        names[i] = pure::decode_token(f.slot[i], nSNP);
    }
    v.attr("names") = names;
    return v;
}

// Read a named frequency NumericVector back into a FreqVec.
inline pure::FreqVec freqvec_from_R(const NumericVector& v, int nSNP) {
    pure::FreqVec f;
    StringVector names = v.names();
    f.val.resize(v.size());
    f.slot.resize(v.size());
    for (int i = 0; i < v.size(); i++) {
        f.val[i] = v[i];
        f.slot[i] = pure::encode_token(as<std::string>(names[i]), nSNP);
    }
    return f;
}

// --- historical API, delegating to the pure core ----------------------------

/**
 * @brief Re-estimate the recomb (no-LD) joint block-pair frequencies.
 * @param Resp  Per-chromosome responsibilities for the recomb population.
 * @param Block Character vector of "left+right" haplotype tokens (one per chromosome).
 * @param nSNP  Number of SNPs per block.
 * @return Named numeric vector of joint recomb frequencies, keyed by "left+right".
 */
inline NumericVector CRecombFreq(NumericVector Resp, StringVector Block, int nSNP) {
    std::vector<double> resp(Resp.begin(), Resp.end());
    return freqvec_to_R(pure::recomb_freq(resp, encode_tokens(Block, nSNP), nSNP), nSNP);
}

/**
 * @brief Re-estimate the linkage (full-LD) joint block-pair frequencies.
 * @param Resp  Per-chromosome responsibilities for the linkage population.
 * @param Block Character vector of "left+right" haplotype tokens (one per chromosome).
 * @param nSNP  Number of SNPs per block.
 * @return Named numeric vector of joint linkage frequencies, keyed by "left+right".
 */
inline NumericVector CLinkageFreq(NumericVector Resp, StringVector Block, int nSNP) {
    std::vector<double> resp(Resp.begin(), Resp.end());
    return freqvec_to_R(pure::linkage_freq(resp, encode_tokens(Block, nSNP), nSNP), nSNP);
}

/**
 * @brief One EM step (E + M) of the two-population LD mixture model.
 * @param par List with props1, props2, prob0, inds (tokens), blocksize.
 * @return Updated list (r1, r2, props1, props2, prob0, inds, blocksize).
 */
inline List CupdateModel(List par) {
    StringVector inds = as<StringVector>(par["inds"]);
    int blocksize = as<int>(par["blocksize"]);
    const int fullcombs = pow(2, 2 * blocksize);

    std::vector<int> codes = encode_tokens(inds, blocksize);
    pure::EMParams p;
    p.props1 = freqvec_from_R(as<NumericVector>(par["props1"]), blocksize);
    p.props2 = freqvec_from_R(as<NumericVector>(par["props2"]), blocksize);
    p.prob0 = as<double>(par["prob0"]);

    pure::EMParams np = pure::update_model(p, codes, blocksize);

    // The historical list also carried the unnormalised npr1/npr2 (never read
    // back by the EM, but part of this function's contract).
    const std::vector<double> p1 = p.props1.by_slot(fullcombs);
    const std::vector<double> p2 = p.props2.by_slot(fullcombs);
    NumericVector npr1(codes.size()), npr2(codes.size());
    for (std::size_t i = 0; i < codes.size(); i++) {
        npr1[i] = p.prob0 * p1[codes[i]];
        npr2[i] = (1 - p.prob0) * p2[codes[i]];
    }

    return List::create(Named("r1") = npr1, Named("r2") = npr2,
                        Named("props1") = freqvec_to_R(np.props1, blocksize),
                        Named("props2") = freqvec_to_R(np.props2, blocksize),
                        Named("prob0") = np.prob0,
                        Named("inds") = inds, Named("blocksize") = blocksize);
}

/**
 * @brief Fit the LD mixture model to one pair of SNP-blocks (in memory).
 * @param dat        Character matrix, 2 columns (left/right allele strings), one row per chromosome.
 * @param imaxSteps  Maximum EM iterations.
 * @param dprob0     Initial mixture probability.
 * @param iblocksize SNPs per block.
 * @return list(logNoLD, prob, r1).
 */
inline List LDmixtureModel_mem(RObject dat, int imaxSteps, double dprob0, int iblocksize) {
    CharacterMatrix cdat = transpose(as<CharacterMatrix>(dat));
    StringVector inds = concatenate(cdat(0, _), cdat(1, _), "+");

    pure::MixtureFit fit =
        pure::ld_mixture_model(encode_tokens(inds, iblocksize), imaxSteps, dprob0, iblocksize);

    NumericVector R1(fit.r1.begin(), fit.r1.end());
    R1.attr("names") = inds;
    return List::create(Named("logNoLD") = fit.logNoLD, Named("prob") = fit.prob,
                        Named("r1") = R1);
}

/**
 * @brief Windowed reduction: per-chromosome mean recomb probability over a window.
 * @param mat Chromosomes x models responsibility matrix.
 * @param sel 1-based column indices of the models within the window.
 * @return Per-chromosome mean recomb probability (length = number of chromosomes).
 */
inline NumericVector getProbs(RObject mat, RObject sel) {
    Eigen::MatrixXd dmat = as<Eigen::MatrixXd>(mat);
    Eigen::VectorXd dsel = as<Eigen::VectorXd>(sel);
    Eigen::MatrixXd dmatsel = Eigen::MatrixXd::Zero(dmat.rows(), dsel.size());
    for (int i = 0; i < dsel.size(); i++) dmatsel.col(i) = dmat.col(dsel[i] - 1);
    return wrap(dmatsel.rowwise().mean());
}

}  // namespace recombclust
#endif  // RECOMBCLUST_CORE_HPP
