// recombclust_windows.cpp — Path B: windowed recomb-probability reduction (ModelsProb).
//
// Tiles the region into windows of width `window` and, for each window, averages the
// per-chromosome recomb responsibilities (r1) over the models whose genomic range
// CONTAINS the window (window ⊆ [pos[left], pos[right]] — the type="within" overlap of
// the reference pipeline). Reads the selected model columns from HDF5 out-of-core and
// writes one ModelsProb column per non-empty window (chromosomes x windows). Replaces
// getRecombProb_hdf5.R / CgetProbsHdf5. All I/O via the BigDataStatMeth C++ API.

// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <RcppEigen.h>
#include <BigDataStatMeth.hpp>

using namespace Rcpp;

namespace {

// Read R-view column j of an open dataset into `out` (length nrow).
inline void read_col(BigDataStatMeth::hdf5Dataset* d, hsize_t j, hsize_t nrow,
                     std::vector<double>& out) {
    std::vector<hsize_t> off = {j, 0}, cnt = {1, nrow}, st = {1, 1}, bl = {1, 1};
    d->readDatasetBlock(off, cnt, st, bl, out.data());
}

}  // namespace

//' Windowed recomb-probability matrix (Path B / ModelsProb)
//'
//' @param filename HDF5 file holding the models matrix.
//' @param models_group,models_dataset Models matrix (chromosomes x models).
//' @param positions SNP positions (bp), one per SNP/block-start index.
//' @param pairs Integer matrix (models x 2) of 1-based block indices (left, right),
//'   as returned by \code{run_region}.
//' @param range_start,range_end Genomic region tiled into windows.
//' @param window Window width (bp).
//' @param out_group,out_dataset Destination for ModelsProb.
//' @return List: filename, group, dataset, nrow (chromosomes), ncol (non-empty
//'   windows) and \code{windows} (matrix of kept window start/end).
//' @keywords internal
//' @export
// [[Rcpp::export]]
List windows_hdf5(std::string filename, std::string models_group,
                     std::string models_dataset, NumericVector positions,
                     IntegerMatrix pairs, double range_start, double range_end,
                     double window, std::string out_group, std::string out_dataset) {
    try {
        H5::Exception::dontPrint();

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> models(
            new BigDataStatMeth::hdf5Dataset(filename, models_group, models_dataset, false));
        models->openDataset();
        if (models->getDatasetptr() == nullptr)
            throw std::runtime_error("cannot open " + models_group + "/" + models_dataset);
        const hsize_t nChrom = models->nrows_r();
        const int nM = pairs.nrow();

        // Model genomic ranges [pos[left], pos[right]] (block indices are 1-based).
        std::vector<double> Lpos(nM), Rpos(nM);
        for (int k = 0; k < nM; ++k) {
            Lpos[k] = positions[pairs(k, 0) - 1];
            Rpos[k] = positions[pairs(k, 1) - 1];
        }

        // Windows tiling [range_start, range_end]; keep only non-empty ones.
        std::vector<double> keptStart, keptEnd;
        std::vector<std::vector<int>> keptSel;
        for (double ws = range_start; ws <= range_end; ws += window) {
            const double we = ws + window - 1.0;
            std::vector<int> sel;
            for (int k = 0; k < nM; ++k)
                if (Lpos[k] <= ws && Rpos[k] >= we) sel.push_back(k);   // window within model
            if (!sel.empty()) {
                keptStart.push_back(ws);
                keptEnd.push_back(we);
                keptSel.push_back(std::move(sel));
            }
        }
        const int nW = (int)keptSel.size();
        if (nW == 0) throw std::runtime_error("no non-empty window in the region");

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> out(
            new BigDataStatMeth::hdf5Dataset(filename, out_group, out_dataset, true));
        out->createDataset((int)nChrom, nW, "real");
        out->openDataset();

        std::vector<double> col(nChrom);
        std::vector<hsize_t> st = {1, 1}, bl = {1, 1};
        for (int w = 0; w < nW; ++w) {
            Eigen::VectorXd acc = Eigen::VectorXd::Zero(nChrom);
            for (int k : keptSel[w]) {
                read_col(models.get(), (hsize_t)k, nChrom, col);
                acc += Eigen::Map<Eigen::VectorXd>(col.data(), nChrom);
            }
            acc /= (double)keptSel[w].size();               // per-window mean r1
            std::vector<hsize_t> woff = {(hsize_t)w, 0}, wcnt = {nChrom, 1};
            out->writeDatasetBlock(Rcpp::wrap(acc), woff, wcnt, st, bl, false);
        }

        NumericMatrix windows(nW, 2);
        for (int w = 0; w < nW; ++w) { windows(w, 0) = keptStart[w]; windows(w, 1) = keptEnd[w]; }
        colnames(windows) = CharacterVector::create("start", "end");

        return List::create(Named("filename") = filename, Named("group") = out_group,
                            Named("dataset") = out_dataset, Named("nrow") = (int)nChrom,
                            Named("ncol") = nW, Named("windows") = windows);
    } catch (std::exception& ex) {
        Rf_error("windows_hdf5: %s", ex.what());
    } catch (...) {
        Rf_error("windows_hdf5: unknown C++ exception");
    }
    return R_NilValue;
}
