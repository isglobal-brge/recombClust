# recombClust

## Description

recombClust is an R package that groups chromosomes by their **recombination history**. Two chromosomes can carry the same alleles yet differ in history: one may descend from an ancestral haplotype that has not recombined across a region, another may be a mosaic produced by recombination. recombClust reads that difference from phased genotypes and answers two questions for a genomic region: *are there groups of chromosomes with different recombination histories?* and *which part of the region drives the difference?* A typical application is characterising the two orientations of a chromosomal inversion.

Recombination history is modelled with a mixture that, given a pair of SNP-blocks, separates chromosomes into two populations, one with high Linkage Disequilibrium (LD) and low recombination (linkage) and another with low LD and high recombination. The classification of many SNP-block pairs in a region groups chromosomes into clusters with different recombination history. The package takes phased genotype data as input.

From release **v2.0.0** the computational core is reimplemented in C++ on HDF5 files through the [BigDataStatMeth](https://github.com/isglobal-brge/BigDataStatMeth) API, so the analysis runs **out-of-core** and scales to large cohorts without loading the data into memory.

The methods described in the paper **"Identifying chromosomal subpopulations based on their recombination histories advances the study of the genetic basis of phenotypic traits" by Carlos Ruiz-Arenas et al. (2020)**[[1]](#1) available at **[Genome Research](https://genome.cshlp.org/content/early/2020/11/17/gr.258301.119)** can be found in release [v 1.0.0](https://github.com/isglobal-brge/recombClust/releases/tag/v1.0.0).

Release [v 0.0.9](https://github.com/isglobal-brge/recombClust/releases/tag/v0.0.9) contains the methods described in preprint **"Historical recombination variability contributes to deciphering the genetic basis of phenotypic traits" by Carlos Ruiz-Arenas et al. (2019)**[[2]](#2) publicly available at [bioRxiv](https://www.biorxiv.org/content/10.1101/792747v1)

### Pre-requisites

recombClust builds on [BigDataStatMeth](https://github.com/isglobal-brge/BigDataStatMeth) (its C++/HDF5 backend) and on a few Bioconductor packages for reading genomic files. Install them first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Rhdf5lib", "SeqArray", "gdsfmt"))
install.packages("igraph")

# BigDataStatMeth is available from GitHub
remotes::install_github("isglobal-brge/BigDataStatMeth")
```

### Installation

recombClust can be installed from the GitHub repository:

```r
remotes::install_github("isglobal-brge/recombClust")
```

### Usage

```r
library(recombClust)

# 1. Import phased genotypes (VCF/GDS/BED) into an HDF5 haplotype matrix
h <- hdf5omics_import("region.vcf.gz", out_file = "region.h5")

# 2. Classify chromosomes by their recombination history
res <- recombClust(h$file, h$positions, clusters = 2)
res$class          # cluster of each chromosome
plot(res)          # PCA of the recombination profiles

# 3. Locate the region driving the signal
probs   <- recomb_probabilities(h$file, h$positions)
regions <- define_Regions(probs)
regions            # detected LD blocks (coordinates)
plot(regions)      # LD-block heat map
```

The vignette (`vignette("recombClust")`) is a step-by-step guide with a worked example: it explains the options of each function, what they change, and how to read and interpret the results. Missing genotypes can be imputed beforehand with `hdf5omics_impute()`.

### Coming from a previous version

The analysis is unchanged from earlier versions; a few entry points were renamed or merged into cleaner ones:

| Task | Earlier versions | Now |
|---|---|---|
| Load phased genotypes | `getData()` | `hdf5omics_import()` |
| Classify chromosomes | `runRecombClust()` | `recombClust()` |
| Windowed recombination probabilities | `getRecombProb()`, `computeRecombProbabilities()` | `recomb_probabilities()` |
| Detect regions | `define_Regions()` | `define_Regions()` (unchanged) |

The core model functions keep their names (`recombFreq()`, `linkageFreq()`, `updateModel()`, `LDmixtureModel()`, `getProbs()`). The hand-written HDF5 helpers of the older versions are no longer part of the package: that layer is now provided by [BigDataStatMeth](https://github.com/isglobal-brge/BigDataStatMeth), which is what lets the analysis run out-of-core.

### References
<a id="1">[1]</a> 
Ruiz-Arenas C, Cáceres A, López M, Pelegrí-Sisó D, González J, González JR. Identifying chromosomal subpopulations based on their recombination histories advances the study of the genetic basis of phenotypic traits. Genome Res. 2020 Nov 17. doi: [10.1101/gr.258301.119](https://doi.org/10.1101/gr.258301.119). Epub ahead of print. PMID: 33203765.
\  

<a id="2">[2]</a> 
Historical recombination variability contributes to deciphering the genetic basis of phenotypic traits
Carlos Ruiz-Arenas, Alejandro Cáceres, Marcos López, Dolors Pelegrí-Sisó, Josefa González, Juan R. González
bioRxiv 792747; doi: [https://doi.org/10.1101/792747](https://doi.org/10.1101/792747)
