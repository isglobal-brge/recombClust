// recombclust_stream.cpp — stream the LD mixture models to HDF5 block-wise.
//
// Phase 3 (write path). Runs the per-pair EM and writes each model's r1 as a
// column block into a pre-sized HDF5 dataset that is opened once, so the full
// models matrix is never held in memory (PATH2 / out-of-core). Uses BigDataStatMeth's
// block write (writeDatasetBlock) via the recombclust_io.hpp R-view idiom.
//
// This entry point takes a list of already-formed SNP-block pairs (each a 2-column
// character matrix). Block pairing from SNP positions moves to C++ in a later step;
// here we validate the streaming write path against the in-memory result.

// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <RcppEigen.h>
#include <BigDataStatMeth.hpp>
#include "recombclust_core.h"
#include "recombclust_io.h"

using namespace Rcpp;

//' Stream LD mixture models to an HDF5 dataset, block-wise
//'
//' @param filename HDF5 file (created here; in the real pipeline it exists from ingestion).
//' @param group Output group.
//' @param dataset Output dataset name (chromosomes x models).
//' @param pairs List of 2-column character matrices, one per SNP-block pair.
//' @param maxSteps,prob0,blocksize EM parameters.
//' @return Named list with filename, group, dataset, nrow, ncol.
//' @keywords internal
//' @export
// [[Rcpp::export]]
List stream_models(std::string filename, std::string group, std::string dataset,
                      List pairs, int maxSteps = 100, double prob0 = 0.5,
                      int blocksize = 2) {
    try {
        H5::Exception::dontPrint();

        int nPairs = pairs.size();
        if (nPairs == 0) throw std::runtime_error("empty pairs list");

        // Number of chromosomes from the first pair (rows of the char matrix).
        CharacterMatrix first = as<CharacterMatrix>(pairs[0]);
        int nChrom = first.nrow();

        // Create the file (real pipeline: already exists from ingestion).
        {
            std::unique_ptr<BigDataStatMeth::hdf5File> f(
                new BigDataStatMeth::hdf5File(filename, true));
            f->createFile();
        }

        // Pre-size the models dataset and open it once for streaming.
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> models(
            new BigDataStatMeth::hdf5Dataset(filename, group, dataset, true));
        models->createDataset(nChrom, nPairs, "real");
        models->openDataset();
        if (models->getDatasetptr() == nullptr)
            throw std::runtime_error("cannot open " + group + "/" + dataset);

        for (int j = 0; j < nPairs; j++) {
            List m = recombclust::LDmixtureModel_mem(pairs[j], maxSteps, prob0, blocksize);
            Eigen::VectorXd r1 = as<Eigen::VectorXd>(as<NumericVector>(m["r1"]));
            recombclust::write_col(models.get(), j, r1);
        }

        return List::create(Named("filename") = filename, Named("group") = group,
                            Named("dataset") = dataset, Named("nrow") = nChrom,
                            Named("ncol") = nPairs);
    } catch (std::exception& ex) {
        Rf_error("stream_models: %s", ex.what());
    } catch (...) {
        Rf_error("stream_models: unknown C++ exception");
    }
    return R_NilValue;
}

//' Stream an arbitrary numeric matrix to HDF5 column by column (orientation test)
//'
//' Writes M one R-view column at a time via write_col (block write on an open
//' dataset). Used to prove the R-view <-> HDF5 transposition is handled correctly.
//' @param filename,group,dataset Output location.
//' @param M Numeric matrix.
//' @return Named list with the descriptor.
//' @keywords internal
//' @export
// [[Rcpp::export]]
List stream_matrix(std::string filename, std::string group, std::string dataset,
                      Rcpp::NumericMatrix M) {
    try {
        H5::Exception::dontPrint();
        Eigen::Map<Eigen::MatrixXd> em(REAL(M), M.nrow(), M.ncol());
        {
            std::unique_ptr<BigDataStatMeth::hdf5File> f(
                new BigDataStatMeth::hdf5File(filename, true));
            f->createFile();
        }
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> d(
            new BigDataStatMeth::hdf5Dataset(filename, group, dataset, true));
        d->createDataset(em.rows(), em.cols(), "real");
        d->openDataset();
        for (int j = 0; j < em.cols(); j++)
            recombclust::write_col(d.get(), j, em.col(j));
        return List::create(Named("filename") = filename, Named("group") = group,
                            Named("dataset") = dataset, Named("nrow") = (int)em.rows(),
                            Named("ncol") = (int)em.cols());
    } catch (std::exception& ex) {
        Rf_error("stream_matrix: %s", ex.what());
    } catch (...) {
        Rf_error("stream_matrix: unknown C++ exception");
    }
    return R_NilValue;
}
