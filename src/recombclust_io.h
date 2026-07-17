// recombclust_io.hpp — the single recombClust-specific HDF5 I/O idiom.
//
// We deliberately keep ONLY what BigDataStatMeth does not already provide cleanly:
// a streaming column-block write that encapsulates the R-view <-> HDF5 (row/col
// major) transposition. Whole-matrix read/write is NOT here on purpose — BDSM
// already offers it (R: `hdf5_create_matrix` / `hdf5_matrix`; C++: `createDataset`
// + `writeDataset` + `readDatasetBlock`). Our job is to leverage BDSM, not to
// re-wrap it. This idiom is a candidate to migrate into BigDataStatMeth later.
#ifndef RECOMBCLUST_IO_HPP
#define RECOMBCLUST_IO_HPP

#include <BigDataStatMeth.hpp>

namespace recombclust {

/**
 * @brief Write an R-view column into an already-open HDF5 dataset as a block.
 *
 * Streaming write idiom: the caller creates + opens the dataset once and writes
 * columns one at a time as they are produced, so the full matrix is never
 * materialised in memory (PATH2 / out-of-core). This is the only I/O helper
 * recombClust adds; whole-matrix read/write is left to BigDataStatMeth.
 *
 * Encapsulates the R-view <-> HDF5 transposition: HDF5 stores transposed, and for a
 * NumericVector @c writeDatasetBlock uses count = {1, vCount[0]}, so the vector
 * length goes in @c cnt[0] and @c offset[0]=j selects the column (a "file row").
 * This matches the file-order convention BDSM's readers use, so streaming and BDSM's
 * own whole-matrix write yield the identical on-disk layout (05_check_orientation.R).
 *
 * @param d   Pointer to an already-opened BigDataStatMeth::hdf5Dataset.
 * @param j   0-based R-view column index to write.
 * @param col Column values (length = number of R-view rows).
 */
inline void write_col(BigDataStatMeth::hdf5Dataset* d, int j, const Eigen::VectorXd& col) {
    std::vector<hsize_t> off = {(hsize_t)j, 0}, cnt = {(hsize_t)col.size(), 1},
                         st = {1, 1}, bl = {1, 1};
    d->writeDatasetBlock(Rcpp::wrap(col), off, cnt, st, bl, false);
}

}  // namespace recombclust
#endif  // RECOMBCLUST_IO_HPP
