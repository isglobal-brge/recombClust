// recombclust_pca.cpp — PCA of the models matrix via BigDataStatMeth's C++ SVD.
//
// Calls BigDataStatMeth::RcppbdSVD_hdf5 directly (LinkingTo), NOT the R prcomp S3
// method (DESIGN 4.6: compute on the C++ API). The SVD (centered, unscaled) of the
// models matrix X (chromosomes x models) is X = U D V'; PCA scores = U D and
// sdev = d / sqrt(n-1), matching stats::prcomp. Only the k-means downstream stays in R.

// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <RcppEigen.h>
#include <BigDataStatMeth.hpp>

using namespace Rcpp;

namespace {

// Read a full HDF5 dataset into an Eigen matrix (R-view). Thin use of BDSM's
// readDatasetBlock; the file stores transposed, hence count {ncol_r, nrow_r}.
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

}  // namespace

//' PCA of an HDF5 matrix via BigDataStatMeth's C++ SVD
//'
//' Computes the SVD of \code{group/dataset} (centered, unscaled) with
//' BigDataStatMeth's C++ \code{RcppbdSVD_hdf5} (out-of-core; writes
//' \code{SVD/<dataset>/{u,d,v}}), then forms the PCA scores \code{U D} and standard
//' deviations \code{d / sqrt(n-1)} (\code{stats::prcomp} contract). Observations are
//' the matrix rows (R-view).
//'
//' @param filename HDF5 file.
//' @param group,dataset Location of the matrix (rows = observations).
//' @param ncomp Number of principal components.
//' @return List with \code{scores} (observations x ncomp) and \code{sdev} (length ncomp).
//' @keywords internal
//' @export
// [[Rcpp::export]]
List pca_hdf5(std::string filename, std::string group, std::string dataset, int ncomp) {
    try {
        H5::Exception::dontPrint();

        // Centered, unscaled SVD on the BigDataStatMeth C++ API.
        BigDataStatMeth::RcppbdSVD_hdf5(filename, group, dataset,
                                        /*k*/ 2, /*q*/ 1, /*nev*/ ncomp,
                                        /*bcenter*/ true, /*bscale*/ false,
                                        /*dthreshold*/ 0.0, /*bforce*/ true,
                                        /*asRowMajor*/ false, R_NilValue, R_NilValue);

        const std::string svd_group = "SVD/" + dataset;
        Eigen::MatrixXd U = read_rview(filename, svd_group, "u");   // n x r
        Eigen::MatrixXd D = read_rview(filename, svd_group, "d");   // 1 x r

        const int n = (int)U.rows();
        int r = std::min((int)U.cols(), (int)D.size());
        if (ncomp > 0) r = std::min(r, ncomp);

        Eigen::VectorXd d(r);
        for (int i = 0; i < r; ++i) d[i] = D(0, i);

        Eigen::MatrixXd scores = U.leftCols(r) * d.asDiagonal();     // n x r
        Eigen::VectorXd sdev = d / std::sqrt((double)(n - 1));

        return List::create(Named("scores") = wrap(scores), Named("sdev") = wrap(sdev));
    } catch (std::exception& ex) {
        Rf_error("pca_hdf5: %s", ex.what());
    } catch (...) {
        Rf_error("pca_hdf5: unknown C++ exception");
    }
    return R_NilValue;
}
