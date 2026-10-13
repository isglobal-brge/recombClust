# recombClust 2.1.0

* The EM model core is reimplemented in pure C++ (block-pair models no longer
  round-trip through R), and the block-pair model stage is parallelised with OpenMP:
  at the largest benchmarked problem (2,000 chromosomes; 159,080 block-pair models)
  the serial run takes 1,122 s and the 8-thread run 187 s, with identical clustering.
* Corrected default execution parameters: the number of blocks is selected
  automatically from the problem's dimensions (the shallowest blocking with at most
  1e7 elements per block), and intermediate HDF5 datasets default to gzip level 1.
  Explicit arguments override both, as before.

# recombClust 2.0.0

* Computational core reimplemented in C++ on HDF5 files through the
  [BigDataStatMeth](https://github.com/isglobal-brge/BigDataStatMeth) API: the
  haplotypes, the fitted models and the intermediate matrices stay on disk, so the
  analysis runs out-of-core and scales to large cohorts.
* New main entry point `recombClust()` (returns an S3 `"recombClust"` object with
  `print`, `summary` and `plot` methods).
* Region detection with `define_Regions()` (BigLD/CLQD), plus `plot.recombRegions()`
  LD-block heat map and `plotRecombTracks()` genomic-track view.
* Embedded `hdf5omics` layer for data handling: `hdf5omics_import()` (VCF/GDS/BED to
  an HDF5 haplotype matrix), `hdf5omics_filter_maf()` (MAF filtering) and
  `hdf5omics_impute()` (reference-panel-free LD-kNNi imputation of 0/1 haplotypes).
* New practical vignette guiding the user through a complete analysis.

# recombClust 1.0.0

* Methods described in Ruiz-Arenas et al. (2020), *Genome Research*
  (doi:10.1101/gr.258301.119).

# recombClust 0.0.9

* Methods described in the Ruiz-Arenas et al. (2019) bioRxiv preprint
  (doi:10.1101/792747).
