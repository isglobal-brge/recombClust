// recombclust_region.cpp — the region orchestrator (Phase 3, Family-B / C++).
//
// Given a haplotype matrix in HDF5 (chromosomes x SNP, 0/1), the SNP positions and
// the model parameters, this runs the whole per-region core in one C++ call:
//   block_pairs -> per-pair EM (LDmixtureModel_mem) -> stream r1 to the models
//   matrix (write_col, column-block, never held whole in memory).
// Replaces runLDmixtureModel.R's GenomicRanges + bplapply loop.
//
// Out-of-core note (docs/DESIGN.md): the OUTPUT models matrix (chromosomes x nPairs,
// nPairs large) is the big object and is STREAMED. The INPUT haplotype matrix
// (chromosomes x nSNP, comparatively small) is read whole once via BDSM's
// readDatasetBlock. If a region's haplotypes no longer fit, the input read would
// also move to per-block; the output streaming already scales.

// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <RcppEigen.h>
#include <BigDataStatMeth.hpp>
#include "recombclust_core.h"
#include "recombclust_pairs.h"
#include "recombclust_io.h"

using namespace Rcpp;

namespace {

/**
 * @brief Read a full HDF5 dataset into an Eigen matrix in R-view order.
 *
 * Thin inline use of BigDataStatMeth's readDatasetBlock primitive (no wrapper kept
 * in the shared header — whole-matrix reads are BDSM's job). File stores transposed,
 * so the block count is {ncol_r, nrow_r}.
 */
Eigen::MatrixXd read_matrix(const std::string& fn, const std::string& grp,
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

/// @brief Append the '0'/'1' character for a 0/1-coded haplotype value.
inline void append_allele(std::string& s, double v) {
    s.push_back(static_cast<char>('0' + static_cast<int>(std::llround(v))));
}

}  // namespace

//' Run the recombClust region core (block pairing + EM + streamed models matrix)
//'
//' @param filename HDF5 file holding the haplotype matrix (and receiving the models).
//' @param hap_group,hap_dataset Location of the haplotype matrix (chromosomes x SNP, 0/1).
//' @param positions Numeric vector of SNP positions (bp), sorted ascending, length = nSNP.
//' @param out_group,out_dataset Where to write the models matrix (chromosomes x models).
//' @param blockSize SNPs per block (default 2).
//' @param distance Max bp between paired blocks' starts (default 1e4).
//' @param maxSteps,prob0 EM parameters.
//' @return List: filename, group, dataset, nrow (chromosomes), ncol (models), and
//'   \code{pairs} (integer matrix of 1-based block-index pairs, one row per model).
//' @keywords internal
//' @export
// [[Rcpp::export]]
List run_region(std::string filename, std::string hap_group, std::string hap_dataset,
                   NumericVector positions, std::string out_group, std::string out_dataset,
                   int blockSize = 2, double distance = 1e4,
                   int maxSteps = 100, double prob0 = 0.5) {
    try {
        H5::Exception::dontPrint();

        // Input haplotypes (read whole; see out-of-core note above).
        Eigen::MatrixXd H = read_matrix(filename, hap_group, hap_dataset);
        const int nChrom = static_cast<int>(H.rows());
        const int nSNP = static_cast<int>(H.cols());
        if (static_cast<int>(positions.size()) != nSNP)
            throw std::runtime_error("positions length != number of SNP columns");

        // Block pairs (0-based block indices).
        std::vector<double> pos(positions.begin(), positions.end());
        std::vector<std::pair<int, int>> pairs =
            recombclust::block_pairs(pos, blockSize, distance);
        const int nPairs = static_cast<int>(pairs.size());
        if (nPairs == 0) throw std::runtime_error("no SNP-block pairs for this region");

        // Pre-size the models matrix and open it once for streaming.
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> models(
            new BigDataStatMeth::hdf5Dataset(filename, out_group, out_dataset, true));
        models->createDataset(nChrom, nPairs, "real");
        models->openDataset();
        if (models->getDatasetptr() == nullptr)
            throw std::runtime_error("cannot open " + out_group + "/" + out_dataset);

        IntegerMatrix pairsOut(nPairs, 2);

        for (int k = 0; k < nPairs; ++k) {
            const int i = pairs[k].first, j = pairs[k].second;
            pairsOut(k, 0) = i + 1;
            pairsOut(k, 1) = j + 1;

            // Build the left/right block haplotype strings per chromosome.
            CharacterMatrix dat(nChrom, 2);
            for (int c = 0; c < nChrom; ++c) {
                std::string left, right;
                left.reserve(blockSize);
                right.reserve(blockSize);
                for (int b = 0; b < blockSize; ++b) append_allele(left, H(c, i + b));
                for (int b = 0; b < blockSize; ++b) append_allele(right, H(c, j + b));
                dat(c, 0) = left;
                dat(c, 1) = right;
            }

            List m = recombclust::LDmixtureModel_mem(dat, maxSteps, prob0, blockSize);
            Eigen::VectorXd r1 = as<Eigen::VectorXd>(as<NumericVector>(m["r1"]));
            recombclust::write_col(models.get(), k, r1);
        }

        colnames(pairsOut) = CharacterVector::create("left", "right");
        return List::create(Named("filename") = filename, Named("group") = out_group,
                            Named("dataset") = out_dataset, Named("nrow") = nChrom,
                            Named("ncol") = nPairs, Named("pairs") = pairsOut);
    } catch (std::exception& ex) {
        Rf_error("run_region: %s", ex.what());
    } catch (...) {
        Rf_error("run_region: unknown C++ exception");
    }
    return R_NilValue;
}
