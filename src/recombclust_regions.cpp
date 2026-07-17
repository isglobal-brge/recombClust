// recombclust_regions.cpp — HDF5-backed region detection (Path B, v2.0).
//
// Two exports that route all HDF5 I/O and the correlation through the
// BigDataStatMeth C++ API (DESIGN 4.6); only the graph primitives inside BigLD
// are delegated to igraph (see recombclust_bigld.hpp — the documented exception).
//
//   cor_hdf5     ModelsProb (chromosomes x windows) -> windows x windows
//                   correlation via BigDataStatMeth::RcppbdCorr_hdf5_single
//                   (real Pearson, out-of-core), after dropping windows with
//                   <=5% non-missing (the define_Regions column filter). Fixes
//                   the legacy getCorrelationMatrix bug (it returned covariance).
//   bigld_hdf5   reads that correlation matrix (via BDSM) and runs the BigLD
//                   port -> data.frame of LD blocks.

// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <RcppEigen.h>
#include <BigDataStatMeth.hpp>
#include "recombclust_bigld.h"

using namespace Rcpp;

namespace {

// Read a full HDF5 dataset into an Eigen matrix (R-view). File stores transposed,
// hence readDatasetBlock count {ncol_r, nrow_r}.
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

// Write an Eigen matrix (R-view) to HDF5, overwriting any existing dataset.
void write_rview(const std::string& fn, const std::string& grp, const std::string& ds,
                 const Eigen::MatrixXd& M) {
    std::unique_ptr<BigDataStatMeth::hdf5Dataset> d(
        new BigDataStatMeth::hdf5Dataset(fn, grp, ds, true));
    d->createDataset((int)M.rows(), (int)M.cols(), "real");
    d->openDataset();
    d->writeDataset(Rcpp::wrap(M));
}

}  // namespace

//' Windowed correlation matrix of ModelsProb (out-of-core, BigDataStatMeth)
//'
//' Drops windows with 5% or fewer non-missing chromosomes (the
//' \code{define_Regions} column filter), then computes the windows x windows
//' Pearson correlation of the surviving windows with
//' \code{BigDataStatMeth::RcppbdCorr_hdf5_single} (pairwise-complete
//' observations, out-of-core). Optionally stores absolute values. All HDF5 I/O
//' and the correlation go through the BigDataStatMeth C++ API.
//'
//' @param filename HDF5 file holding ModelsProb; the result is written here.
//' @param group,dataset Location of ModelsProb (chromosomes x windows).
//' @param out_group,out_dataset Destination for the correlation matrix.
//' @param absval If TRUE (default) store \code{abs(cor)}.
//' @param minNonNA Minimum fraction of non-missing entries to keep a window.
//' @return List with the correlation descriptor and \code{kept} (1-based indices
//'   of the surviving windows, for mapping back to window positions).
//' @keywords internal
//' @export
// [[Rcpp::export]]
List cor_hdf5(std::string filename, std::string group, std::string dataset,
                 std::string out_group, std::string out_dataset,
                 bool absval = true, double minNonNA = 0.05) {
    try {
        H5::Exception::dontPrint();

        // Per-window (R-column) non-missing fraction via block reads.
        int nChrom, nWin;
        std::vector<int> kept;   // 0-based window indices that survive the filter
        {
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> mp(
                new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false));
            mp->openDataset();
            if (mp->getDatasetptr() == nullptr)
                throw std::runtime_error("cannot open " + group + "/" + dataset);
            nChrom = (int)mp->nrows_r();
            nWin = (int)mp->ncols_r();
            std::vector<double> col(nChrom);
            std::vector<hsize_t> st = {1, 1}, bl = {1, 1};
            for (int j = 0; j < nWin; ++j) {
                std::vector<hsize_t> off = {(hsize_t)j, 0}, cnt = {1, (hsize_t)nChrom};
                mp->readDatasetBlock(off, cnt, st, bl, col.data());
                int nonNA = 0;
                for (int i = 0; i < nChrom; ++i) if (!std::isnan(col[i])) ++nonNA;
                if ((double)nonNA / (double)nChrom > minNonNA) kept.push_back(j);
            }
        }
        if (kept.empty()) throw std::runtime_error("no window passes the non-missing filter");

        // Correlation of ALL windows via BigDataStatMeth (out-of-core), written
        // to the output dataset; then subset to kept windows and optionally abs.
        BigDataStatMeth::RcppbdCorr_hdf5_single(
            filename, group, dataset, /*method*/ "pearson",
            /*use_complete_obs*/ true, /*compute_pvalues*/ false,
            /*block_size*/ 1024, /*bforce*/ true,
            Rcpp::wrap(out_group), Rcpp::wrap(out_dataset), R_NilValue,
            /*trans_x*/ true, R_NilValue);

        Eigen::MatrixXd Cfull = read_rview(filename, out_group, out_dataset);  // nWin x nWin
        const int nk = (int)kept.size();
        bool rewrite = (nk < nWin) || absval;
        if (rewrite) {
            Eigen::MatrixXd C = Cfull(kept, kept);         // Eigen 3.4 index slicing
            if (absval) C = C.cwiseAbs();
            write_rview(filename, out_group, out_dataset, C);
        }

        IntegerVector keptR(nk);
        for (int a = 0; a < nk; ++a) keptR[a] = kept[a] + 1;   // 1-based for R
        return List::create(Named("filename") = filename, Named("group") = out_group,
                            Named("dataset") = out_dataset, Named("nwin") = nk,
                            Named("kept") = keptR);
    } catch (std::exception& ex) {
        Rf_error("cor_hdf5: %s", ex.what());
    } catch (...) {
        Rf_error("cor_hdf5: unknown C++ exception");
    }
    return R_NilValue;
}

//' BigLD region blocks from an on-disk correlation matrix (BigDataStatMeth I/O)
//'
//' Reads the windows x windows correlation matrix through the BigDataStatMeth
//' C++ API and runs the BigLD/CLQD port. \code{win_start}/\code{win_end} give
//' the bp span of each window (already restricted to the windows present in the
//' correlation matrix); \code{win_start} is used as the CLQD SNP position.
//'
//' @param filename HDF5 file holding the correlation matrix.
//' @param group,dataset Location of the correlation matrix.
//' @param win_start,win_end bp span per window (length == ncol of the matrix).
//' @param CLQcut,clstgap,hrstType,hrstParam,CLQmode BigLD/CLQD parameters.
//' @param chrN Chromosome label.
//' @return data.frame(chr, start, end, start.idx, end.idx), sorted by start.idx.
//' @keywords internal
//' @export
// [[Rcpp::export]]
DataFrame bigld_hdf5(std::string filename, std::string group, std::string dataset,
                        NumericVector win_start, NumericVector win_end,
                        double CLQcut = 0.5, double clstgap = 40000,
                        std::string hrstType = "fast", int hrstParam = 200,
                        std::string CLQmode = "maximal", std::string chrN = "chr1") {
    try {
        H5::Exception::dontPrint();
        Eigen::MatrixXd OCM = read_rview(filename, group, dataset);
        if (OCM.rows() != OCM.cols())
            throw std::runtime_error("correlation matrix must be square");
        if (win_start.size() != OCM.cols() || win_end.size() != OCM.cols())
            throw std::runtime_error("win_start/win_end length must equal ncol(matrix)");

        std::vector<double> ws(win_start.begin(), win_start.end());
        std::vector<std::array<int, 2>> blocks =
            recombclust::BigLD(OCM, CLQcut, clstgap, hrstType, hrstParam, CLQmode, ws);

        const int nb = (int)blocks.size();
        CharacterVector chr(nb);
        NumericVector start(nb), end(nb);
        IntegerVector sidx(nb), eidx(nb);
        for (int i = 0; i < nb; ++i) {
            chr[i] = chrN;
            start[i] = win_start[blocks[i][0]];
            end[i] = win_end[blocks[i][1]];
            sidx[i] = blocks[i][0] + 1;
            eidx[i] = blocks[i][1] + 1;
        }
        return DataFrame::create(Named("chr") = chr, Named("start") = start,
                                 Named("end") = end, Named("start.idx") = sidx,
                                 Named("end.idx") = eidx,
                                 Named("stringsAsFactors") = false);
    } catch (std::exception& ex) {
        Rf_error("bigld_hdf5: %s", ex.what());
    } catch (...) {
        Rf_error("bigld_hdf5: unknown C++ exception");
    }
    return DataFrame::create();
}
