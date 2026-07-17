// omics_hdf5.cpp — hdf5omics out-of-core HDF5 helpers (namespace hdf5omics).
//
// Part of the embedded, extraction-ready hdf5omics layer (docs/DESIGN.md 6.3). Does
// NOT depend on recombClust. Uses the BigDataStatMeth C++ API directly (LinkingTo),
// never its R S3 layer (DESIGN 4.6). Files are prefixed omics_* (a movable unit)
// rather than placed in a src/ subdirectory, which R's package build does not compile.

// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <RcppEigen.h>
#include <BigDataStatMeth.hpp>

using namespace Rcpp;

namespace {

/**
 * @brief Copy given columns of an open HDF5 dataset into another, one at a time.
 * @param in Open source dataset.
 * @param out Open, pre-sized destination dataset (nrow x keep0.size()).
 * @param keep0 0-based source column indices, in output order.
 */
void copy_columns(BigDataStatMeth::hdf5Dataset* in, BigDataStatMeth::hdf5Dataset* out,
                  const std::vector<int>& keep0) {
    const hsize_t nrow = in->nrows_r();
    std::vector<double> col(nrow);
    std::vector<hsize_t> st = {1, 1}, bl = {1, 1};
    for (std::size_t c = 0; c < keep0.size(); ++c) {
        const hsize_t j = (hsize_t)keep0[c];
        std::vector<hsize_t> roff = {j, 0}, rcnt = {1, nrow};   // R-view col j == file row j
        in->readDatasetBlock(roff, rcnt, st, bl, col.data());
        std::vector<hsize_t> woff = {(hsize_t)c, 0}, wcnt = {nrow, 1};
        out->writeDatasetBlock(Rcpp::wrap(Eigen::Map<Eigen::VectorXd>(col.data(), nrow)),
                               woff, wcnt, st, bl, false);
    }
}

}  // namespace

//' Subset columns of an HDF5 matrix into a new dataset (out-of-core)
//'
//' Copies the requested columns of an HDF5-backed matrix to a new dataset one column
//' at a time, so the full matrix is never materialised in memory.
//'
//' @param filename HDF5 file.
//' @param in_group,in_dataset Source matrix (rows x SNPs).
//' @param out_group,out_dataset Destination dataset.
//' @param keep 1-based column indices to keep, in output order.
//' @return List descriptor of the new dataset (filename, group, dataset, nrow, ncol).
//' @keywords internal
//' @export
// [[Rcpp::export]]
List hdf5omics_subset_cols(std::string filename, std::string in_group, std::string in_dataset,
                           std::string out_group, std::string out_dataset,
                           IntegerVector keep) {
    try {
        H5::Exception::dontPrint();
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> in(
            new BigDataStatMeth::hdf5Dataset(filename, in_group, in_dataset, false));
        in->openDataset();
        if (in->getDatasetptr() == nullptr)
            throw std::runtime_error("cannot open " + in_group + "/" + in_dataset);
        const hsize_t nrow = in->nrows_r();

        std::vector<int> keep0(keep.size());
        for (int i = 0; i < keep.size(); ++i) keep0[i] = keep[i] - 1;

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> out(
            new BigDataStatMeth::hdf5Dataset(filename, out_group, out_dataset, true));
        out->createDataset((int)nrow, (int)keep0.size(), "real");
        out->openDataset();
        copy_columns(in.get(), out.get(), keep0);

        return List::create(Named("filename") = filename, Named("group") = out_group,
                            Named("dataset") = out_dataset, Named("nrow") = (int)nrow,
                            Named("ncol") = (int)keep0.size());
    } catch (std::exception& ex) {
        Rf_error("hdf5omics_subset_cols: %s", ex.what());
    } catch (...) {
        Rf_error("hdf5omics_subset_cols: unknown C++ exception");
    }
    return R_NilValue;
}

//' Filter SNP columns by minor allele frequency (0/1 haplotypes, out-of-core, C++)
//'
//' Computes each SNP's frequency of "1" with BigDataStatMeth's C++ out-of-core column
//' means (\code{get_HDF5_colMeans}), forms the phased-haplotype MAF
//' \code{min(freq, 1 - freq)} and copies the columns with \code{MAF >= threshold}
//' to a new dataset. All computation is on the BigDataStatMeth C++ API (no R S3 layer).
//' Unlike \code{BigDataStatMeth::filter_maf}, which assumes 0/1/2 dosage coding.
//'
//' @param filename HDF5 file.
//' @param in_group,in_dataset Source haplotype matrix (chromosomes x SNP, 0/1).
//' @param out_group,out_dataset Destination for the filtered matrix.
//' @param threshold Minimum minor allele frequency to keep a SNP.
//' @return List: filename, group, dataset, nrow, ncol, \code{keep} (1-based kept
//'   column indices) and \code{maf} (MAF of the kept SNPs).
//' @keywords internal
//' @export
// [[Rcpp::export]]
List hdf5omics_filter_maf_hdf5(std::string filename, std::string in_group,
                               std::string in_dataset, std::string out_group,
                               std::string out_dataset, double threshold) {
    try {
        H5::Exception::dontPrint();
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> in(
            new BigDataStatMeth::hdf5Dataset(filename, in_group, in_dataset, false));
        in->openDataset();
        if (in->getDatasetptr() == nullptr)
            throw std::runtime_error("cannot open " + in_group + "/" + in_dataset);
        const hsize_t nrow = in->nrows_r();

        // Per-SNP frequency of "1" via BigDataStatMeth's C++ out-of-core column means.
        Eigen::VectorXd freq1 =
            BigDataStatMeth::get_HDF5_colMeans(in.get(), false, R_NilValue, R_NilValue);

        std::vector<int> keep0;
        std::vector<double> mafKept;
        for (int c = 0; c < freq1.size(); ++c) {
            const double maf = std::min(freq1[c], 1.0 - freq1[c]);
            if (maf >= threshold) { keep0.push_back(c); mafKept.push_back(maf); }
        }
        if (keep0.empty())
            throw std::runtime_error("no SNP passes the MAF threshold");

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> out(
            new BigDataStatMeth::hdf5Dataset(filename, out_group, out_dataset, true));
        out->createDataset((int)nrow, (int)keep0.size(), "real");
        out->openDataset();
        copy_columns(in.get(), out.get(), keep0);

        IntegerVector keepR(keep0.size());
        for (std::size_t i = 0; i < keep0.size(); ++i) keepR[i] = keep0[i] + 1;   // 1-based

        return List::create(Named("filename") = filename, Named("group") = out_group,
                            Named("dataset") = out_dataset, Named("nrow") = (int)nrow,
                            Named("ncol") = (int)keep0.size(), Named("keep") = keepR,
                            Named("maf") = wrap(mafKept));
    } catch (std::exception& ex) {
        Rf_error("hdf5omics_filter_maf_hdf5: %s", ex.what());
    } catch (...) {
        Rf_error("hdf5omics_filter_maf_hdf5: unknown C++ exception");
    }
    return R_NilValue;
}
