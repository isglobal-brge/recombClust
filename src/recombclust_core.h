// recombclust_core.hpp — pure Eigen/Rcpp numerical core of the LD mixture model.
//
// Ported (Phase 2) from the current recombClust package's C++ core, consolidated
// as header-only `inline` functions so the dev loop (Rcpp::sourceCpp) compiles in a
// single translation unit. NO HDF5 / BigDataStatMeth dependency here — this is the
// mixture-model mathematics only, validated 1:1 against the pure-R v1.0.0 reference
// (recombFreq.R / linkageFreq.R / updateModel.R). See docs/DESIGN.md sections 3.1, 8.
//
// Numerical behaviour is preserved verbatim from the source (rounding to 1e-15, the
// 1e-5 floor + renormalisation in the linkage model, greedy allele matching, etc.).
#ifndef RECOMBCLUST_CORE_HPP
#define RECOMBCLUST_CORE_HPP

#include <RcppEigen.h>
#include <map>
#include <unordered_map>
#include <string>
#include <functional>
#include <algorithm>
#include <cmath>

namespace recombclust {

using namespace Rcpp;

// ---------------------------------------------------------------------------
// String / combination helpers (from CUtils)
// ---------------------------------------------------------------------------

// Concatenate two string vectors position by position with `par` in the middle.
inline CharacterVector concatenate(StringVector x, StringVector y, std::string par) {
    CharacterVector res(x.size());
    for (int i = 0; i < x.size(); i++)
        res[i] = ((Rcpp::as<std::string>(x[i]) + par) + Rcpp::as<std::string>(y[i]));
    return res;
}

// Duplicate array content by block of x vector.
inline StringVector getNamesbyBlocks(StringVector x) {
    int vlen = x.length();
    int vlength = pow(2, x.length());
    StringVector vnames;
    for (int i = 0; i < vlength; i++) vnames.push_back(x[i % vlen]);
    return vnames;
}

// Duplicate n times each position of array x.
inline StringVector getNamesbyWord(StringVector x) {
    int vlen = x.length();
    int vlength = pow(2, x.length());
    StringVector vnames;
    for (int i = 0; i < vlength; i++) vnames.push_back(x[(int)(i / vlen)]);
    return vnames;
}

// Generate all binary strings of length n (recursive).
inline void getAllBinaryStrings(int n, int arr[], int i, StringVector* k) {
    if (i == n) {
        std::string word;
        for (int j = 0; j < n; j++) word += std::to_string(arr[j]);
        k->push_back(word);
        return;
    }
    arr[i] = 0;
    getAllBinaryStrings(n, arr, i + 1, k);
    arr[i] = 1;
    getAllBinaryStrings(n, arr, i + 1, k);
}

// Split allele token by inserting `symbol` at position nSNP.
inline void splitAllBinaryStrings(StringVector* combs, int nSNP, char symbol) {
    int clenght = combs->length();
    for (int i = 0; i < clenght; i++) {
        std::string comb = (std::string)combs->at(i);
        comb = comb.substr(0, nSNP) + symbol + comb.substr(nSNP, nSNP);
        combs->at(i) = comb;
    }
}

// Sum of responsibilities by substring [start, start+length).
inline std::map<std::string, double> getFreqSum(NumericVector R, StringVector Block,
                                                int start, int length) {
    std::map<std::string, double> mapfreq;
    for (int i = 0; i < R.length(); i++)
        mapfreq[as<std::string>(Block(i)).substr(start, length)] += R[i];
    return mapfreq;
}

// Sum of responsibilities over all 2^ilength combinations of a sub-block.
inline std::map<std::string, double> getFreqAllCombs(NumericVector R, StringVector Block,
                                                     int start, int ilength) {
    std::map<std::string, double> mapfreq;
    int arr[ilength];
    StringVector combs;
    getAllBinaryStrings(ilength, arr, 0, &combs);
    for (int i = 0; i < pow(2, ilength); i++) mapfreq[as<std::string>(combs(i))] = 0;
    for (int i = 0; i < R.length(); i++)
        mapfreq[as<std::string>(Block(i)).substr(start, ilength)] += R[i];
    return mapfreq;
}

// ---------------------------------------------------------------------------
// Map helpers (from CMapUtils)
// ---------------------------------------------------------------------------

inline std::map<std::string, double> VectortoOrderedMap(NumericVector v) {
    StringVector namesv = v.names();
    std::map<std::string, double> mapv;
    for (int i = 0; i < namesv.size(); ++i) mapv[(std::string)namesv[i]] = v[i];
    return mapv;
}

inline NumericVector getNumericVectorfromStringVector(std::map<std::string, double> mapv,
                                                      StringVector strvalues) {
    int vsize = strvalues.length();
    NumericVector v(vsize);
    for (int i = 0; i < vsize; i++) v[i] = mapv[(std::string)strvalues[i]];
    return v;
}

// Flip <key,value> and sort by value descending -> multimap<value, key>.
inline std::multimap<double, std::string, std::greater<double>>
sortMultimapbyValue(std::map<std::string, double> sNPFreq) {
    std::multimap<double, std::string, std::greater<double>> mmsortFreq;
    for (std::pair<std::string, double> element : sNPFreq)
        mmsortFreq.insert(std::pair<double, std::string>(element.second, element.first));
    return mmsortFreq;
}

// ---------------------------------------------------------------------------
// Linear algebra helpers (from CQRSolve, CFlatMat, CgetMatrixFromMap)
// ---------------------------------------------------------------------------

inline Eigen::VectorXd CQR_Solve(Eigen::MatrixXd A, Eigen::VectorXd b) {
    return A.colPivHouseholderQr().solve(b);
}

// Column-major flatten of a matrix into a vector.
inline NumericVector flatmatcm(NumericMatrix x) {
    double xcol = x.ncol(), xrow = x.nrow();
    NumericVector fc(xrow * xcol);
    for (double i = 0; i < xrow; i++)
        for (double j = 0; j < xcol; j++) fc[i + j * xrow] = x(i, j);
    return fc;
}

// Build the AA / bb linear system and names vector from a frequency map.
inline void getMatrixfromMap(NumericMatrix* AA, StringVector* names, NumericVector* bb,
                             std::map<std::string, double> datamap, int ipos) {
    double imaxf = 999999;
    int i = ipos;
    std::multimap<double, std::string, std::greater<double>> sdatamap =
        sortMultimapbyValue(datamap);
    bb->at(ipos) = 1;
    for (auto x : sdatamap) {
        if (i == ipos) imaxf = x.first;
        AA->at(i, ipos) = ((double)x.first * (-1)) / imaxf;
        AA->at(ipos, i) = 1;
        names->at(i) = x.second;
        i++;
    }
}

// ---------------------------------------------------------------------------
// Recomb (no-LD) model — independent block frequencies (from CrecombFreq)
// ---------------------------------------------------------------------------

/**
 * @brief Re-estimate the recomb (no-LD) joint block-pair frequencies.
 *
 * Under the recomb model the two SNP-blocks are independent, so the joint
 * frequency factorises into the outer product of the left and right marginals.
 * The marginals are obtained by solving a small QR system built from the
 * responsibility-weighted block counts.
 *
 * @param Resp  Per-chromosome responsibilities for the recomb population.
 * @param Block Character vector of "left+right" haplotype tokens (one per chromosome).
 * @param nSNP  Number of SNPs per block.
 * @return Named numeric vector of joint recomb frequencies, keyed by "left+right".
 */
inline NumericVector CRecombFreq(NumericVector Resp, StringVector Block, int nSNP) {
    std::map<std::string, double> freqLeft, freqRight;
    int ncombs = pow(2, nSNP);

    freqLeft = getFreqAllCombs(Resp, Block, 0, nSNP);
    freqRight = getFreqAllCombs(Resp, Block, nSNP + 1, nSNP);

    NumericMatrix AA(2 * ncombs, 2 * ncombs);
    NumericMatrix bb(1, 2 * ncombs);
    StringVector names(2 * ncombs);

    AA.fill_diag(1);
    getMatrixfromMap(&AA, &names, &bb, freqLeft, 0);
    getMatrixfromMap(&AA, &names, &bb, freqRight, ncombs);

    Eigen::VectorXd props = CQR_Solve(as<Eigen::MatrixXd>(AA), as<Eigen::VectorXd>(bb));
    Eigen::MatrixXd ans = props.head(ncombs) * props.tail(ncombs).adjoint();

    StringVector namesf = concatenate(getNamesbyBlocks(names[Range(0, ncombs - 1)]),
                                      getNamesbyWord(names[Range(ncombs, names.length() - 1)]),
                                      "+");
    NumericVector ansf = flatmatcm(wrap(ans));
    ansf.attr("names") = namesf;
    for (NumericVector::iterator i = ansf.begin(); i != ansf.end(); ++i)
        *i = round(*i * 1000000000000000.0) / 1000000000000000.0;
    return ansf;
}

// ---------------------------------------------------------------------------
// Linkage (full-LD) model — greedy allele matching + QR (from CLinkageFreq)
// ---------------------------------------------------------------------------

/**
 * @brief Re-estimate the linkage (full-LD) joint block-pair frequencies.
 *
 * Under the linkage model the right block is a deterministic function of the left,
 * so haplotypes are matched greedily by decreasing frequency and the resulting
 * frequencies are obtained from a small QR system, then floored at 1e-5 and
 * renormalised for numerical stability.
 *
 * @param Resp  Per-chromosome responsibilities for the linkage population.
 * @param Block Character vector of "left+right" haplotype tokens (one per chromosome).
 * @param nSNP  Number of SNPs per block.
 * @return Named numeric vector of joint linkage frequencies, keyed by "left+right".
 */
inline NumericVector CLinkageFreq(NumericVector Resp, StringVector Block, int nSNP) {
    int ncombs = pow(2, nSNP);
    int sizeblock = Block(0).size();

    std::map<std::string, double> freqs = getFreqSum(Resp, Block, 0, sizeblock);
    std::multimap<double, std::string, std::greater<double>> sfreqs =
        sortMultimapbyValue(freqs);

    int arr[nSNP];
    StringVector leftlevs, rightlevs, leftLevsvec, rightLevsvec;
    std::map<std::string, double> selAlleles;

    getAllBinaryStrings(nSNP, arr, 0, &leftlevs);
    rightLevsvec = leftlevs;
    rightlevs = rightLevsvec;
    leftLevsvec = leftlevs;

    std::multimap<double, std::string>::iterator itm;

    while (sfreqs.size() > 0 && leftLevsvec.size() > 0) {
        itm = sfreqs.begin();
        std::string allele = itm->second;
        if ((std::find(leftLevsvec.begin(), leftLevsvec.end(), allele.substr(0, nSNP)) !=
             leftLevsvec.end()) &&
            (std::find(rightLevsvec.begin(), rightLevsvec.end(), allele.substr(nSNP + 1, nSNP)) !=
             rightLevsvec.end())) {
            selAlleles.insert(selAlleles.end(), std::pair<std::string, double>(allele, 0));
            auto itl = std::find(leftLevsvec.begin(), leftLevsvec.end(), allele.substr(0, nSNP));
            if (itl != leftLevsvec.end()) leftLevsvec.erase(itl);
            auto itr =
                std::find(rightLevsvec.begin(), rightLevsvec.end(), allele.substr(nSNP + 1, nSNP));
            if (itr != rightLevsvec.end()) rightLevsvec.erase(itr);
        }
        sfreqs.erase(itm);
    }

    if (leftLevsvec.length() > 0) {
        StringVector leftright = concatenate(leftLevsvec, rightLevsvec, "+");
        for (int i = 0; i < leftright.length(); i++)
            selAlleles.insert(selAlleles.end(),
                              std::pair<std::string, double>((std::string)leftright[i], 0));
    }

    for (std::pair<std::string, double> x : selAlleles) {
        std::map<std::string, double>::iterator ito = freqs.find(x.first);
        if (ito != freqs.end()) selAlleles[x.first] = ito->second;
    }

    NumericMatrix AA(ncombs, ncombs);
    NumericMatrix bb(1, ncombs);
    int fullcombs = pow(2, 2 * nSNP);
    StringVector names(ncombs);

    AA.fill_diag(1);
    getMatrixfromMap(&AA, &names, &bb, selAlleles, 0);

    Eigen::VectorXd props = Eigen::VectorXd::Zero(fullcombs);
    props.segment(0, ncombs) = CQR_Solve(as<Eigen::MatrixXd>(AA), as<Eigen::VectorXd>(bb));
    props = props.unaryExpr([](double v) { return v > 1e-5 ? v : 1e-5; });
    props = props / props.array().sum();

    StringVector namesf;
    StringVector combs;
    getAllBinaryStrings(2 * nSNP, arr, 0, &combs);
    splitAllBinaryStrings(&combs, nSNP, '+');

    for (int i = 0; i < fullcombs; i++) {
        if (i < names.size()) {
            namesf.push_back(names[i]);
            combs.erase(std::remove(combs.begin(), combs.end(), names[i]), combs.end());
        } else {
            namesf.push_back(combs[i - names.size()]);
        }
    }

    Eigen::Map<Eigen::MatrixXd> ans(props.data(), props.rows(), props.cols());
    NumericVector ansf = flatmatcm(wrap(ans));
    ansf.attr("names") = namesf;
    for (NumericVector::iterator i = ansf.begin(); i != ansf.end(); ++i)
        *i = round(*i * 1000000000000000.0) / 1000000000000000.0;
    return ansf;
}

// ---------------------------------------------------------------------------
// EM step (from CUpdateModel)
// ---------------------------------------------------------------------------

/**
 * @brief One EM step (E + M) of the two-population LD mixture model.
 *
 * E-step: computes per-chromosome responsibilities from the current mixing
 * proportion and the two frequency tables. M-step: updates the mixing proportion
 * (mean recomb responsibility) and re-estimates both frequency tables via
 * CRecombFreq / CLinkageFreq.
 *
 * @param par List with props1, props2, prob0, inds (tokens), blocksize.
 * @return Updated list (r1, r2, props1, props2, prob0, inds, blocksize).
 */
inline List CupdateModel(List par) {
    StringVector inds = as<StringVector>(par["inds"]);
    int blocksize = as<int>(par["blocksize"]);

    NumericVector npr1 = as<double>(par["prob0"]) *
        getNumericVectorfromStringVector(VectortoOrderedMap(as<NumericVector>(par["props1"])), inds);
    NumericVector npr2 = (1 - as<double>(par["prob0"])) *
        getNumericVectorfromStringVector(VectortoOrderedMap(as<NumericVector>(par["props2"])), inds);

    NumericVector R1 = npr1 / (npr1 + npr2);
    NumericVector R2 = npr2 / (npr1 + npr2);

    double npprob0 = mean(R1);
    NumericVector npprops1 = CRecombFreq(R1, inds, blocksize);
    NumericVector npprops2 = CLinkageFreq(R2, inds, blocksize);

    return List::create(Named("r1") = npr1, Named("r2") = npr2, Named("props1") = npprops1,
                        Named("props2") = npprops2, Named("prob0") = npprob0,
                        Named("inds") = inds, Named("blocksize") = blocksize);
}

// ---------------------------------------------------------------------------
// In-memory EM driver (from CLDmixtureModel, HDF5 write path stripped)
// ---------------------------------------------------------------------------
// dat: character matrix, 2 columns (left-block allele string, right-block allele
//      string), one row per chromosome. Returns list(logNoLD, prob, r1, annot).
/**
 * @brief Fit the LD mixture model to one pair of SNP-blocks (in memory).
 *
 * Initialises the two frequency tables, then iterates CupdateModel until the L2
 * change in the frequency tables plus |delta prob0| falls below 1e-5 (or maxSteps,
 * or a degenerate mixing proportion of 0/1). Returns the recomb log-likelihood, the
 * final mixing proportion, and the per-chromosome recomb responsibilities r1.
 *
 * @param dat        Character matrix, 2 columns (left/right allele strings), one row per chromosome.
 * @param imaxSteps  Maximum EM iterations.
 * @param dprob0     Initial mixture probability.
 * @param iblocksize SNPs per block.
 * @return list(logNoLD, prob, r1).
 */
inline List LDmixtureModel_mem(RObject dat, int imaxSteps, double dprob0, int iblocksize) {
    CharacterMatrix cdat = transpose(as<CharacterMatrix>(dat));

    StringVector inds = concatenate(cdat(0, _), cdat(1, _), "+");

    NumericVector r1(inds.length());
    r1.fill(1);

    NumericVector propsRecomb = CRecombFreq(r1, inds, iblocksize);
    NumericVector propsLink = CLinkageFreq(r1, inds, iblocksize);

    r1 = getNumericVectorfromStringVector(VectortoOrderedMap(propsRecomb), inds);
    NumericVector r2 = getNumericVectorfromStringVector(VectortoOrderedMap(propsLink), inds);

    r1.attr("names") = inds;
    r2.attr("names") = inds;

    double LoglikeRecomb = sum(log(r1));

    r1 = dprob0 * r1;
    r2 = (1 - dprob0) * r2;

    List params = List::create(Named("r1") = r1, Named("r2") = r2, Named("props1") = propsRecomb,
                               Named("props2") = propsLink, Named("prob0") = dprob0,
                               Named("inds") = inds, Named("blocksize") = iblocksize);

    double tol = 1;
    int steps = 1;
    double MINTOL = .00001;

    while (tol > MINTOL && steps <= imaxSteps) {
        List newparams = CupdateModel(params);

        NumericVector parprops1 = as<NumericVector>(params["props1"]);
        NumericVector nparprops1 = as<NumericVector>(newparams["props1"]);
        NumericVector parprops2 = as<NumericVector>(params["props2"]);
        NumericVector nparprops2 = as<NumericVector>(newparams["props2"]);

        NumericVector difprops1 = (parprops1 - nparprops1);
        NumericVector difprops2 = (parprops2 - nparprops2);

        tol = sqrt(std::inner_product(difprops1.begin(), difprops1.end(), difprops1.begin(), 0.0) +
                   std::inner_product(difprops2.begin(), difprops2.end(), difprops2.begin(), 0.0) +
                   std::abs(as<double>(params["prob0"]) - as<double>(newparams["prob0"])));

        params = clone(newparams);
        steps = steps + 1;

        if (as<double>(newparams["prob0"]) == 1 || as<double>(newparams["prob0"]) == 0) break;
    }

    r1 = as<double>(params["prob0"]) *
        getNumericVectorfromStringVector(VectortoOrderedMap(as<NumericVector>(params["props1"])), inds);
    r2 = (1 - as<double>(params["prob0"])) *
        getNumericVectorfromStringVector(VectortoOrderedMap(as<NumericVector>(params["props2"])), inds);

    NumericVector R1 = r1 / (r1 + r2);
    R1.attr("names") = inds;

    return List::create(Named("logNoLD") = LoglikeRecomb, Named("prob") = as<double>(params["prob0"]),
                        Named("r1") = R1);
}

// ---------------------------------------------------------------------------
// Windowed reduction (from CgetProbs): mean over selected model columns.
// mat: chromosomes x models. sel: 1-based column indices of the models that fall
// in a window. Returns the per-chromosome mean recomb probability for that window.
// ---------------------------------------------------------------------------
/**
 * @brief Windowed reduction: per-chromosome mean recomb probability over a window.
 *
 * Selects the model columns that fall in a genomic window and averages them
 * row-wise, giving one recomb probability per chromosome for that window.
 *
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
