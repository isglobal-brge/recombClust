// recombclust_pairs.hpp — SNP-block pairing (pure geometry, no HDF5).
//
// Ported to C++ from runLDmixtureModel.R (recombClust 1.0.0), replacing the
// GenomicRanges/findOverlaps logic. Given sorted SNP positions, forms blocks of
// `blockSize` consecutive SNPs and pairs blocks that are close but non-overlapping.
// See docs/DESIGN.md section 9.3 (block pairing -> C++).
#ifndef RECOMBCLUST_PAIRS_HPP
#define RECOMBCLUST_PAIRS_HPP

#include <RcppEigen.h>
#include <vector>
#include <utility>

namespace recombclust {

/**
 * @brief Enumerate SNP-block pairs to test with the LD mixture model.
 *
 * Blocks are runs of @p blockSize consecutive SNPs: block @c i spans SNP indices
 * @c [i, i+blockSize-1] and, genomically, @c [pos[i], pos[i+blockSize-1]]. There are
 * @c nSNP-blockSize+1 blocks. A pair @c (i,j) with @c i<j is kept iff the two blocks
 * are non-overlapping and block @c j starts within @p distance of block @c i:
 *   - non-overlap:  @c pos[j] > pos[i+blockSize-1]
 *   - within range: @c pos[j] <= pos[i] + distance - 1
 * This reproduces exactly the `findOverlaps(resize(GRblocks, distance), GRblocks)`
 * minus `findOverlaps(GRblocks)` set used by the reference R pipeline (validated in
 * tests/harness/06_check_pairs.R). Positions are assumed sorted ascending.
 *
 * @param pos        SNP positions (bp), sorted ascending, length nSNP.
 * @param blockSize  Number of SNPs per block (>= 1).
 * @param distance   Maximum bp from block i's start to block j's start to pair them.
 * @return Vector of 0-based block-index pairs (left block i, right block j), i<j.
 */
inline std::vector<std::pair<int, int>> block_pairs(const std::vector<double>& pos,
                                                    int blockSize, double distance) {
    std::vector<std::pair<int, int>> pairs;
    const int nSNP = static_cast<int>(pos.size());
    const int nBlocks = nSNP - blockSize + 1;
    if (nBlocks <= 1) return pairs;

    for (int i = 0; i < nBlocks; ++i) {
        const double iEnd = pos[i + blockSize - 1];   // block i's last SNP position
        const double wEnd = pos[i] + distance - 1.0;  // resize(distance) window end
        for (int j = i + 1; j < nBlocks; ++j) {
            if (pos[j] > wEnd) break;                 // sorted: no further j in window
            if (pos[j] > iEnd) pairs.emplace_back(i, j);  // non-overlapping -> keep
        }
    }
    return pairs;
}

}  // namespace recombclust
#endif  // RECOMBCLUST_PAIRS_HPP
