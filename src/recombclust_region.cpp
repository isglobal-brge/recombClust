// recombclust_region.cpp — the region orchestrator (Phase 3, Family-B / C++).
//
// Given a haplotype matrix in HDF5 (chromosomes x SNP, 0/1), the SNP positions and
// the model parameters, this runs the whole per-region core in one C++ call:
//   block_pairs -> per-pair EM (pure C++ core, OpenMP over pairs) -> batched
//   streaming of r1 columns to the models matrix (never held whole in memory).
// Replaces runLDmixtureModel.R's GenomicRanges + bplapply loop. Workers touch
// neither the R API nor HDF5; the master thread writes each finished batch.
//
// Out-of-core note (docs/DESIGN.md): the OUTPUT models matrix (chromosomes x nPairs,
// nPairs large) is the big object and is STREAMED. The INPUT haplotype matrix
// (chromosomes x nSNP, comparatively small) is read whole once via BDSM's
// readDatasetBlock. If a region's haplotypes no longer fit, the input read would
// also move to per-block; the output streaming already scales.

// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <RcppEigen.h>
#include <BigDataStatMeth.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "recombclust_core.h"
#include "recombclust_pairs.h"

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

/// @brief 0/1 bit for a 0/1-coded haplotype value (same rounding the old
/// string tokens used); anything else is a corrupt haplotype matrix.
inline int allele_bit(double v) {
    const int b = static_cast<int>(std::llround(v));
    if (b != 0 && b != 1)
        throw std::runtime_error("haplotype matrix is not 0/1-coded (value " +
                                 std::to_string(v) + ")");
    return b;
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
//' @param compression gzip level (0-9) for the on-disk models matrix; -1 keeps
//'   BigDataStatMeth's default (6).
//' @param threads Number of OpenMP threads for the per-pair EM fits; 0 (default)
//'   uses all available cores (respects \code{OMP_NUM_THREADS}). The result is
//'   identical for any thread count: pairs are independent, deterministic fits.
//' @return List: filename, group, dataset, nrow (chromosomes), ncol (models), and
//'   \code{pairs} (integer matrix of 1-based block-index pairs, one row per model).
//' @keywords internal
//' @export
// [[Rcpp::export]]
List run_region(std::string filename, std::string hap_group, std::string hap_dataset,
                   NumericVector positions, std::string out_group, std::string out_dataset,
                   int blockSize = 2, double distance = 1e4,
                   int maxSteps = 100, double prob0 = 0.5, int compression = -1,
                   int threads = 0) {
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
        models->createDataset(nChrom, nPairs, "real", compression);
        models->openDataset();
        if (models->getDatasetptr() == nullptr)
            throw std::runtime_error("cannot open " + out_group + "/" + out_dataset);

        IntegerMatrix pairsOut(nPairs, 2);
        for (int k = 0; k < nPairs; ++k) {
            pairsOut(k, 0) = pairs[k].first + 1;
            pairsOut(k, 1) = pairs[k].second + 1;
        }

        int nthreads = threads;
#ifdef _OPENMP
        if (nthreads <= 0 || nthreads > omp_get_max_threads())
            nthreads = omp_get_max_threads();
#else
        nthreads = 1;
#endif

        // Batched parallel EM: worker threads fit the pure-C++ mixture model
        // (no R API, no HDF5 from workers — both are thread-unsafe) into a
        // plain buffer; the master thread then writes the whole batch as one
        // block (R-view columns k0..k0+nb-1 = file rows, buffer row-major).
        // Pairs are independent deterministic fits, so the result is identical
        // for any thread count and any batch size.
        const int BATCH = 256;
        std::vector<double> buf;
        std::string worker_error;

        for (int k0 = 0; k0 < nPairs; k0 += BATCH) {
            const int nb = std::min(BATCH, nPairs - k0);
            buf.assign((std::size_t)nb * nChrom, 0.0);

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
#endif
            for (int t = 0; t < nb; ++t) {
                try {
                    const int k = k0 + t;
                    const int i = pairs[k].first, j = pairs[k].second;
                    // Encode each chromosome's left+right haplotype pair once
                    // as its lexicographic slot — no strings in the hot path.
                    std::vector<int> codes(nChrom);
                    for (int c = 0; c < nChrom; ++c) {
                        int code = 0;
                        for (int b = 0; b < blockSize; ++b)
                            code = (code << 1) | allele_bit(H(c, i + b));
                        for (int b = 0; b < blockSize; ++b)
                            code = (code << 1) | allele_bit(H(c, j + b));
                        codes[c] = code;
                    }
                    recombclust::pure::MixtureFit m =
                        recombclust::pure::ld_mixture_model(codes, maxSteps, prob0, blockSize);
                    std::copy(m.r1.begin(), m.r1.end(),
                              buf.begin() + (std::size_t)t * nChrom);
                } catch (const std::exception& ex) {
                    // Exceptions must not escape an OpenMP region.
#ifdef _OPENMP
                    #pragma omp critical(recombclustWorkerError)
#endif
                    if (worker_error.empty()) worker_error = ex.what();
                }
            }
            if (!worker_error.empty())
                throw std::runtime_error(std::string("EM worker: ") + worker_error);

            std::vector<hsize_t> off = {(hsize_t)k0, 0},
                                 cnt = {(hsize_t)nb, (hsize_t)nChrom},
                                 st = {1, 1}, bl = {1, 1};
            models->writeDatasetBlock(buf, off, cnt, st, bl);
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
