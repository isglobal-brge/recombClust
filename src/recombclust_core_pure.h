// recombclust_core_pure.h — R-free rewrite of the LD mixture model core.
//
// Phase B of the EM-core restructuring (docs/REFACTOR_BRIEF.md): the per-pair
// hot path with NO R types, so it can later run on worker threads and does not
// pay R allocation/GC/protection costs (measured at 27.6% of the EM stage —
// docs/REFACTOR_NOTES.md R1). The R-facing wrappers in recombclust_core.h
// delegate here and rebuild the R objects (names included) at the boundary.
//
// BIT-FOR-BIT CONTRACT with the previous Rcpp implementation (tag
// v2.0.1-validada). Everything below preserves the exact floating-point
// operation sequence and the exact orderings that shape the arithmetic:
//
//  * Haplotype tokens ("01+10") are encoded once per pair as their
//    lexicographic rank (binary value of the 2L bits, '+' dropped). For
//    fixed-length binary tokens, rank order == lexicographic string order,
//    so every std::map<std::string,double> iteration of the old code maps to
//    ascending-slot iteration over dense arrays (values and FP accumulation
//    order unchanged).
//  * sortMultimapbyValue (multimap<double,string,greater>) == stable_sort by
//    descending value with ties kept in ascending-slot (insertion) order.
//  * The EM tolerance subtracts the props vectors POSITIONALLY, in each
//    iteration's own frequency-descending name order (which changes between
//    iterations). That behaviour is preserved verbatim: props carry their
//    positional (namesf) order in `val` plus the slot of each position.
//  * mean() replicates Rcpp sugar's two-pass long-double algorithm;
//    sum(log(x)) replicates sugar's plain left-to-right double accumulation.
//  * The small Eigen expressions (colPivHouseholderQr, outer product,
//    unaryExpr floor, .array().sum() redux) are kept as the very same Eigen
//    calls, so their internal reduction orders are unchanged.
//  * Rounding keeps `round(x * 1e15) / 1e15` on every emitted frequency.
//
// The only intentional behaviour change: tokens must be strictly binary
// ('0'/'1' plus the '+' separator). The old map-based code silently produced
// garbage-in/garbage-out results for non-binary tokens; encode_token() now
// throws instead. The pipeline guarantees binary haplotypes (imputation runs
// before), so no valid caller is affected.

#ifndef RECOMBCLUST_CORE_PURE_H
#define RECOMBCLUST_CORE_PURE_H

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace recombclust {
namespace pure {

// ---------------------------------------------------------------------------
// Token codec. slot = lexicographic rank of the token among all 2^(2L)
// binary "left+right" combinations (leftmost char = most significant bit).
// ---------------------------------------------------------------------------

inline int encode_token(const std::string& tok, int nSNP) {
    if ((int)tok.size() != 2 * nSNP + 1 || tok[nSNP] != '+')
        throw std::runtime_error("recombclust core: malformed haplotype token '" + tok + "'");
    int slot = 0;
    for (int j = 0; j < (int)tok.size(); ++j) {
        if (j == nSNP) continue;
        const char ch = tok[j];
        if (ch != '0' && ch != '1')
            throw std::runtime_error("recombclust core: non-binary haplotype token '" + tok + "'");
        slot = (slot << 1) | (ch == '1');
    }
    return slot;
}

// Rank of one L-bit side (no separator), e.g. a left or right block string.
inline int encode_side(const std::string& side) {
    int slot = 0;
    for (char ch : side) {
        if (ch != '0' && ch != '1')
            throw std::runtime_error("recombclust core: non-binary block token '" + side + "'");
        slot = (slot << 1) | (ch == '1');
    }
    return slot;
}

inline std::string decode_side(int slot, int nSNP) {
    std::string s(nSNP, '0');
    for (int j = 0; j < nSNP; ++j)
        if (slot & (1 << (nSNP - 1 - j))) s[j] = '1';
    return s;
}

inline std::string decode_token(int slot, int nSNP) {
    return decode_side(slot >> nSNP, nSNP) + "+" + decode_side(slot & ((1 << nSNP) - 1), nSNP);
}

// ---------------------------------------------------------------------------
// Frequency table in positional (namesf) order: val[i] belongs to the
// combination slot[i]. Mirrors the old named NumericVector: positional ops
// (the EM tolerance) use val as-is; map-order lookups go through by_slot().
// ---------------------------------------------------------------------------

struct FreqVec {
    std::vector<double> val;   // frequencies in this iteration's namesf order
    std::vector<int> slot;     // slot (lexicographic rank) of each position

    // VectortoOrderedMap equivalent: dense values in ascending-slot order.
    // (A std::map keyed by token holds exactly these values; building it does
    // no floating-point work, so a permutation copy is bit-identical.)
    std::vector<double> by_slot(int fullcombs) const {
        std::vector<double> v(fullcombs, 0.0);
        for (std::size_t i = 0; i < val.size(); ++i) v[slot[i]] = val[i];
        return v;
    }
};

// Rcpp sugar mean() (REALSXP): R's two-pass long-double algorithm.
inline double mean_rcpp(const std::vector<double>& x) {
    const std::size_t n = x.size();
    long double s = std::accumulate(x.begin(), x.end(), 0.0L);
    s /= n;
    if (std::isfinite((double)s)) {
        long double t = 0.0;
        for (std::size_t i = 0; i < n; ++i) t += x[i] - s;
        s += t / n;
    }
    return (double)s;
}

// sortMultimapbyValue equivalent: (value, slot) pairs sorted by descending
// value, ties in ascending-slot order. Input must be in ascending-slot order
// (as a std::map iteration would be); stable_sort then reproduces the
// multimap<greater> tie behaviour (equal keys keep insertion order).
inline void sort_desc_stable(std::vector<std::pair<double, int>>& fv) {
    std::stable_sort(fv.begin(), fv.end(),
                     [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                         return a.first > b.first;
                     });
}

// getMatrixfromMap equivalent: fill the [ipos, ipos+n) block of the linear
// system from a frequency list already sorted by descending value.
// Statement order preserved (AA(ipos,ipos) is first set to -1, then
// overwritten with 1 by the row assignment, exactly as before).
inline void fill_system(Eigen::MatrixXd& AA, std::vector<int>& order_slots,
                        Eigen::VectorXd& bb,
                        const std::vector<std::pair<double, int>>& sorted_fv, int ipos) {
    double imaxf = 999999;
    int i = ipos;
    bb[ipos] = 1;
    for (const auto& x : sorted_fv) {
        if (i == ipos) imaxf = x.first;
        AA(i, ipos) = (x.first * (-1)) / imaxf;
        AA(ipos, i) = 1;
        order_slots[i] = x.second;
        i++;
    }
}

// ---------------------------------------------------------------------------
// Recomb (no-LD) model — CRecombFreq without R types.
// Resp/codes: per-chromosome responsibilities and full-token slots.
// ---------------------------------------------------------------------------

inline FreqVec recomb_freq(const std::vector<double>& Resp, const std::vector<int>& codes,
                           int nSNP) {
    const int ncombs = pow(2, nSNP);

    // getFreqAllCombs on the left / right sub-blocks: all combos pre-set to 0,
    // then accumulated in chromosome order (same FP order as the old maps).
    std::vector<double> freqLeft(ncombs, 0.0), freqRight(ncombs, 0.0);
    for (std::size_t i = 0; i < Resp.size(); ++i) {
        freqLeft[codes[i] >> nSNP] += Resp[i];
        freqRight[codes[i] & (ncombs - 1)] += Resp[i];
    }

    Eigen::MatrixXd AA = Eigen::MatrixXd::Zero(2 * ncombs, 2 * ncombs);
    AA.diagonal().setOnes();                       // NumericMatrix zero-init + fill_diag(1)
    Eigen::VectorXd bb = Eigen::VectorXd::Zero(2 * ncombs);
    std::vector<int> order_slots(2 * ncombs);

    std::vector<std::pair<double, int>> fv(ncombs);
    for (int s = 0; s < ncombs; ++s) fv[s] = {freqLeft[s], s};
    sort_desc_stable(fv);
    fill_system(AA, order_slots, bb, fv, 0);
    for (int s = 0; s < ncombs; ++s) fv[s] = {freqRight[s], s};
    sort_desc_stable(fv);
    fill_system(AA, order_slots, bb, fv, ncombs);

    Eigen::VectorXd props = AA.colPivHouseholderQr().solve(bb);
    Eigen::MatrixXd ans = props.head(ncombs) * props.tail(ncombs).adjoint();

    // flatmatcm + namesf: element k = ans(k % ncombs, k / ncombs), named
    // left[k % ncombs] "+" right[k / ncombs] (both in frequency-descending order).
    FreqVec out;
    out.val.resize(ncombs * ncombs);
    out.slot.resize(ncombs * ncombs);
    for (int k = 0; k < ncombs * ncombs; ++k) {
        const int r = k % ncombs, c = k / ncombs;
        out.val[k] = round(ans(r, c) * 1000000000000000.0) / 1000000000000000.0;
        out.slot[k] = (order_slots[r] << nSNP) | order_slots[ncombs + c];
    }
    return out;
}

// ---------------------------------------------------------------------------
// Linkage (full-LD) model — CLinkageFreq without R types.
// ---------------------------------------------------------------------------

inline FreqVec linkage_freq(const std::vector<double>& Resp, const std::vector<int>& codes,
                            int nSNP) {
    const int ncombs = pow(2, nSNP);
    const int fullcombs = pow(2, 2 * nSNP);

    // getFreqSum over the whole token: only observed tokens exist in the map.
    std::vector<double> freqs(fullcombs, 0.0);
    std::vector<char> present(fullcombs, 0);
    for (std::size_t i = 0; i < Resp.size(); ++i) {
        freqs[codes[i]] += Resp[i];
        present[codes[i]] = 1;
    }

    // sfreqs: observed tokens in ascending-slot order, then stable descending sort.
    std::vector<std::pair<double, int>> sfreqs;
    sfreqs.reserve(fullcombs);
    for (int s = 0; s < fullcombs; ++s)
        if (present[s]) sfreqs.push_back({freqs[s], s});
    sort_desc_stable(sfreqs);

    // Greedy allele matching by decreasing frequency (erase-from-front walk of
    // the old multimap). selAlleles is a sorted map keyed by slot.
    std::vector<char> leftFree(ncombs, 1), rightFree(ncombs, 1);
    int nLeftFree = ncombs;
    std::vector<char> selected(fullcombs, 0);
    for (const auto& x : sfreqs) {
        if (nLeftFree == 0) break;
        const int lft = x.second >> nSNP, rgt = x.second & (ncombs - 1);
        if (leftFree[lft] && rightFree[rgt]) {
            selected[x.second] = 1;
            leftFree[lft] = 0;
            rightFree[rgt] = 0;
            nLeftFree--;
        }
    }
    if (nLeftFree > 0) {
        // Pair the remaining left/right levels positionally (ascending order),
        // as concatenate(leftLevsvec, rightLevsvec, "+") did.
        std::vector<int> lrem, rrem;
        for (int s = 0; s < ncombs; ++s) {
            if (leftFree[s]) lrem.push_back(s);
            if (rightFree[s]) rrem.push_back(s);
        }
        for (std::size_t i = 0; i < lrem.size(); ++i)
            selected[(lrem[i] << nSNP) | rrem[i]] = 1;
    }

    // Values: freqs where the token was observed, else 0 (the old code looked
    // each selected token up in the observed-frequency map).
    std::vector<std::pair<double, int>> selfv;
    selfv.reserve(ncombs);
    for (int s = 0; s < fullcombs; ++s)
        if (selected[s]) selfv.push_back({present[s] ? freqs[s] : 0.0, s});
    sort_desc_stable(selfv);

    Eigen::MatrixXd AA = Eigen::MatrixXd::Zero(ncombs, ncombs);
    AA.diagonal().setOnes();
    Eigen::VectorXd bb = Eigen::VectorXd::Zero(ncombs);
    std::vector<int> order_slots(ncombs);
    fill_system(AA, order_slots, bb, selfv, 0);

    Eigen::VectorXd props = Eigen::VectorXd::Zero(fullcombs);
    props.segment(0, ncombs) = AA.colPivHouseholderQr().solve(bb);
    props = props.unaryExpr([](double v) { return v > 1e-5 ? v : 1e-5; });
    props = props / props.array().sum();

    // namesf: the ncombs selected tokens in frequency-descending order, then
    // every remaining combination in ascending (lexicographic) order.
    FreqVec out;
    out.val.resize(fullcombs);
    out.slot.resize(fullcombs);
    std::vector<char> named(fullcombs, 0);
    for (int i = 0; i < ncombs; ++i) {
        out.slot[i] = order_slots[i];
        named[order_slots[i]] = 1;
    }
    int k = ncombs;
    for (int s = 0; s < fullcombs; ++s)
        if (!named[s]) out.slot[k++] = s;
    for (int i = 0; i < fullcombs; ++i)
        out.val[i] = round(props[i] * 1000000000000000.0) / 1000000000000000.0;
    return out;
}

// ---------------------------------------------------------------------------
// EM state and step — CupdateModel without R types. The params list becomes a
// POD struct; the dead r1/r2 slots of the old list (recomputed from props by
// every consumer, never read back) are dropped.
// ---------------------------------------------------------------------------

struct EMParams {
    FreqVec props1;   // recomb frequencies (positional order + slots)
    FreqVec props2;   // linkage frequencies
    double prob0;
};

struct MixtureFit {
    double logNoLD;
    double prob;
    std::vector<double> r1;   // per-chromosome recomb responsibility
};

inline EMParams update_model(const EMParams& par, const std::vector<int>& codes,
                             int blocksize) {
    const int fullcombs = pow(2, 2 * blocksize);
    const std::size_t n = codes.size();

    const std::vector<double> p1 = par.props1.by_slot(fullcombs);
    const std::vector<double> p2 = par.props2.by_slot(fullcombs);

    // npr1 = prob0 * lookup; npr2 = (1-prob0) * lookup; R = npr/(npr1+npr2).
    // Each vector is materialised in its own pass, exactly like the Rcpp sugar
    // expressions were: keeping the multiply and the add in separate statements
    // also keeps the compiler from contracting them into an FMA, which would
    // change the last bit (observed on Apple clang / ARM64).
    std::vector<double> npr1(n), npr2(n), R1(n), R2(n);
    for (std::size_t i = 0; i < n; ++i) npr1[i] = par.prob0 * p1[codes[i]];
    for (std::size_t i = 0; i < n; ++i) npr2[i] = (1 - par.prob0) * p2[codes[i]];
    for (std::size_t i = 0; i < n; ++i) R1[i] = npr1[i] / (npr1[i] + npr2[i]);
    for (std::size_t i = 0; i < n; ++i) R2[i] = npr2[i] / (npr1[i] + npr2[i]);

    EMParams out;
    out.prob0 = mean_rcpp(R1);
    out.props1 = recomb_freq(R1, codes, blocksize);
    out.props2 = linkage_freq(R2, codes, blocksize);
    return out;
}

// ---------------------------------------------------------------------------
// In-memory EM driver — LDmixtureModel_mem without R types.
// codes: per-chromosome full-token slots (one per chromosome).
// ---------------------------------------------------------------------------

inline MixtureFit ld_mixture_model(const std::vector<int>& codes, int imaxSteps,
                                   double dprob0, int iblocksize) {
    const int fullcombs = pow(2, 2 * iblocksize);
    const std::size_t n = codes.size();

    std::vector<double> ones(n, 1.0);
    EMParams params;
    params.props1 = recomb_freq(ones, codes, iblocksize);
    params.props2 = linkage_freq(ones, codes, iblocksize);
    params.prob0 = dprob0;

    // Initial responsibilities and the recomb log-likelihood: sum(log(r1)) with
    // sugar's plain left-to-right double accumulation.
    const std::vector<double> p1 = params.props1.by_slot(fullcombs);
    double LoglikeRecomb = 0;
    for (std::size_t i = 0; i < n; ++i) LoglikeRecomb += std::log(p1[codes[i]]);

    double tol = 1;
    int steps = 1;
    const double MINTOL = .00001;

    while (tol > MINTOL && steps <= imaxSteps) {
        EMParams newparams = update_model(params, codes, iblocksize);

        // Positional subtraction in each iteration's own namesf order — the
        // orders differ between iterations; that is the validated behaviour.
        const std::vector<double>& a1 = params.props1.val;
        const std::vector<double>& b1 = newparams.props1.val;
        const std::vector<double>& a2 = params.props2.val;
        const std::vector<double>& b2 = newparams.props2.val;
        std::vector<double> difprops1(a1.size()), difprops2(a2.size());
        for (std::size_t i = 0; i < a1.size(); ++i) difprops1[i] = a1[i] - b1[i];
        for (std::size_t i = 0; i < a2.size(); ++i) difprops2[i] = a2[i] - b2[i];

        tol = sqrt(std::inner_product(difprops1.begin(), difprops1.end(), difprops1.begin(), 0.0) +
                   std::inner_product(difprops2.begin(), difprops2.end(), difprops2.begin(), 0.0) +
                   std::abs(params.prob0 - newparams.prob0));

        params = newparams;
        steps = steps + 1;

        if (newparams.prob0 == 1 || newparams.prob0 == 0) break;
    }

    const std::vector<double> f1 = params.props1.by_slot(fullcombs);
    const std::vector<double> f2 = params.props2.by_slot(fullcombs);

    // Final responsibilities, materialised pass by pass like the sugar code
    // (see the FMA-contraction note in update_model).
    std::vector<double> r1(n), r2(n);
    for (std::size_t i = 0; i < n; ++i) r1[i] = params.prob0 * f1[codes[i]];
    for (std::size_t i = 0; i < n; ++i) r2[i] = (1 - params.prob0) * f2[codes[i]];

    MixtureFit fit;
    fit.logNoLD = LoglikeRecomb;
    fit.prob = params.prob0;
    fit.r1.resize(n);
    for (std::size_t i = 0; i < n; ++i) fit.r1[i] = r1[i] / (r1[i] + r2[i]);
    return fit;
}

}  // namespace pure
}  // namespace recombclust

#endif  // RECOMBCLUST_CORE_PURE_H
