// omics_impute.cpp — hdf5omics missing-genotype imputation for 0/1 phased
// haplotype matrices (namespace hdf5omics; docs/DESIGN.md 6.3). Does NOT depend on
// recombClust. Uses the BigDataStatMeth C++ API (LinkingTo), never its R S3 layer.
//
// Two methods (literature: LinkImpute / LD-kNNi, Money et al. G3 2015; Beagle-style
// reference-panel imputation is not applicable to panel-free phased 0/1 data):
//   "ldknni" (default) — LD-aware k-nearest-neighbour imputation. For each SNP with
//     missing entries it takes the L SNPs most correlated with it (LD = |cor|), then
//     imputes each missing allele from the k nearest haplotypes (taxicab/Hamming
//     distance over those L SNPs) by inverse-distance-weighted majority vote.
//   "mode" — per-SNP major allele (fast baseline; biases MAF, see the literature).
//
// Imputation is intrinsically a preload (PATH1) operation — kNN needs the whole
// region matrix in memory — so the SNP x SNP correlation is computed with
// BigDataStatMeth's IN-MEMORY correlation kernel (RcppbdCorr_matrix_single) on the
// already-loaded matrix, rather than the HDF5 (out-of-core) variant, to avoid an
// nSNP x nSNP temporary dataset. The out-of-core correlation path stays where the
// data is large and not otherwise loaded (BigLD's windows x windows).

// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <RcppEigen.h>
#include <BigDataStatMeth.hpp>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

// Read a full HDF5 dataset into an Eigen matrix (R-view; missing stored as NaN).
Eigen::MatrixXd read_rview(const std::string& fn, const std::string& grp,
                           const std::string& ds) {
    std::unique_ptr<BigDataStatMeth::hdf5Dataset> d(
        new BigDataStatMeth::hdf5Dataset(fn, grp, ds, false));
    d->openDataset();
    if (d->getDatasetptr() == nullptr)
        throw std::runtime_error("cannot open " + grp + "/" + ds);
    const std::size_t nr = d->nrows_r(), nc = d->ncols_r();
    Eigen::MatrixXd M(nr, nc);
    std::vector<hsize_t> off = {0, 0}, cnt = {(hsize_t)nc, (hsize_t)nr},
                         st = {1, 1}, bl = {1, 1};
    d->readDatasetBlock(off, cnt, st, bl, M.data());
    return M;
}

// Per-SNP major allele (0/1) over non-missing entries; -1 if a SNP is all-missing.
std::vector<int> mode_alleles(const Eigen::MatrixXd& H) {
    const int nS = (int)H.cols(), nC = (int)H.rows();
    std::vector<int> mode(nS, -1);
    for (int p = 0; p < nS; ++p) {
        double ones = 0.0, n = 0.0;
        for (int h = 0; h < nC; ++h) {
            double v = H(h, p);
            if (!std::isnan(v)) { ones += v; ++n; }
        }
        if (n > 0) mode[p] = (ones / n >= 0.5) ? 1 : 0;
    }
    return mode;
}

}  // namespace

//' Impute missing entries of a 0/1 haplotype matrix (LD-kNNi or mode)
//'
//' Fills missing entries (NaN) of an HDF5 haplotype matrix (chromosomes x SNP,
//' 0/1) and writes the completed matrix to a new dataset. \code{method = "ldknni"}
//' (default) uses LD-aware k-nearest-neighbour imputation (LinkImpute; Money et al.
//' 2015): for each SNP with missing data it takes the \code{L} most-correlated SNPs
//' and imputes each missing allele from the \code{k} nearest haplotypes (Hamming
//' distance over those SNPs) by inverse-distance-weighted majority vote.
//' \code{method = "mode"} uses the per-SNP major allele (fast baseline).
//'
//' The SNP x SNP correlation (LD) is computed with BigDataStatMeth's in-memory
//' correlation kernel; all HDF5 I/O goes through the BigDataStatMeth C++ API.
//'
//' @param filename HDF5 file.
//' @param in_group,in_dataset Source haplotype matrix (chromosomes x SNP, 0/1, NaN = missing).
//' @param out_group,out_dataset Destination for the completed matrix.
//' @param method "ldknni" (default) or "mode".
//' @param k Number of nearest neighbours (LD-kNNi; default 5).
//' @param L Number of most-correlated SNPs used for the distance (LD-kNNi; default 20).
//' @return List: filename, group, dataset, nrow, ncol, \code{n_missing} (imputed
//'   entries) and \code{method}.
//' @keywords internal
//' @export
// [[Rcpp::export]]
List hdf5omics_impute_hdf5(std::string filename, std::string in_group, std::string in_dataset,
                           std::string out_group, std::string out_dataset,
                           std::string method = "ldknni", int k = 5, int L = 20) {
    try {
        H5::Exception::dontPrint();
        if (method != "ldknni" && method != "mode")
            throw std::runtime_error("method must be 'ldknni' or 'mode'");

        Eigen::MatrixXd H = read_rview(filename, in_group, in_dataset);   // nChrom x nSNP
        const int nC = (int)H.rows(), nS = (int)H.cols();
        Eigen::MatrixXd Himp = H;

        // Missing entries per SNP.
        std::vector<std::vector<int>> missRows(nS);      // haplotypes missing at each SNP
        long n_missing = 0;
        for (int p = 0; p < nS; ++p)
            for (int h = 0; h < nC; ++h)
                if (std::isnan(H(h, p))) { missRows[p].push_back(h); ++n_missing; }

        std::vector<int> mode = mode_alleles(H);

        if (method == "mode") {
            for (int p = 0; p < nS; ++p)
                for (int h : missRows[p])
                    Himp(h, p) = (mode[p] < 0) ? 0.0 : (double)mode[p];
        } else {
            // LD: SNP x SNP correlation via BigDataStatMeth's in-memory kernel.
            BigDataStatMeth::corr_result cr =
                BigDataStatMeth::RcppbdCorr_matrix_single(H, "pearson", true, false, false,
                                                          R_NilValue);
            if (!cr.bcomputed) throw std::runtime_error("correlation computation failed");
            Eigen::MatrixXd C = cr.correlation_matrix.cwiseAbs();        // nSNP x nSNP, |cor|

            const double cbase = 1.0;      // LinkImpute distance constant (avoid zeros)
            for (int p = 0; p < nS; ++p) {
                if (missRows[p].empty()) continue;

                // L SNPs most correlated with p (finite |cor|, excluding p itself).
                std::vector<int> cand;
                cand.reserve(nS);
                for (int q = 0; q < nS; ++q)
                    if (q != p && std::isfinite(C(p, q)) && C(p, q) > 0.0) cand.push_back(q);
                const int Lp = std::min((int)cand.size(), L);
                // stable_sort: ties in |cor| keep ascending SNP index (matches R order()).
                std::stable_sort(cand.begin(), cand.end(),
                                 [&](int a, int b) { return C(p, a) > C(p, b); });
                std::vector<int> topL(cand.begin(), cand.begin() + Lp);

                for (int h : missRows[p]) {
                    // Distance from h to every haplotype known at p, over topL SNPs.
                    std::vector<std::pair<double, double>> cand_ds;   // (distance, allele)
                    cand_ds.reserve(nC);
                    for (int s = 0; s < nC; ++s) {
                        if (s == h) continue;
                        double vs = H(s, p);
                        if (std::isnan(vs)) continue;                 // must know allele at p
                        double sumAbs = 0.0; int nAvail = 0;
                        for (int q : topL) {
                            double a = H(h, q), b = H(s, q);
                            if (!std::isnan(a) && !std::isnan(b)) { sumAbs += std::abs(a - b); ++nAvail; }
                        }
                        double d = cbase + (nAvail > 0 ? sumAbs / nAvail : 1.0);
                        cand_ds.emplace_back(d, vs);
                    }
                    if (cand_ds.empty()) {                            // no neighbour: mode
                        Himp(h, p) = (mode[p] < 0) ? 0.0 : (double)mode[p];
                        continue;
                    }
                    const int kk = std::min((int)cand_ds.size(), k);
                    // stable_sort: ties in distance keep ascending haplotype index
                    // (cand_ds is built in ascending s), matching the R reference.
                    std::stable_sort(cand_ds.begin(), cand_ds.end(),
                                     [](const std::pair<double, double>& a,
                                        const std::pair<double, double>& b) { return a.first < b.first; });
                    double num = 0.0, den = 0.0;
                    for (int i = 0; i < kk; ++i) {
                        double w = 1.0 / cand_ds[i].first;
                        num += w * cand_ds[i].second;
                        den += w;
                    }
                    Himp(h, p) = (num / den >= 0.5) ? 1.0 : 0.0;      // weighted majority vote
                }
            }
        }

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> out(
            new BigDataStatMeth::hdf5Dataset(filename, out_group, out_dataset, true));
        out->createDataset(nC, nS, "real");
        out->openDataset();
        out->writeDataset(Rcpp::wrap(Himp));

        return List::create(Named("filename") = filename, Named("group") = out_group,
                            Named("dataset") = out_dataset, Named("nrow") = nC,
                            Named("ncol") = nS, Named("n_missing") = (double)n_missing,
                            Named("method") = method);
    } catch (std::exception& ex) {
        Rf_error("hdf5omics_impute_hdf5: %s", ex.what());
    } catch (...) {
        Rf_error("hdf5omics_impute_hdf5: unknown C++ exception");
    }
    return R_NilValue;
}
