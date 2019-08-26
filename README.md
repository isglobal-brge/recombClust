# recombClust

## Description

recombClust is a R package that groups chromosomes by their recombination history. Recombination history is based on a mixture model that, given a pair of SNP-blocks, separates chromosomes in two populations, one with high Linkage Disequilibrium (LD)  and low recombination (linkage) and another with low LD and high recombination. The method use the classification of several SNP-block pairs in a region to group chromosomes in clusters with different recombination history. This package takes as input genotype phased data.

### Pre-requisites

recombClust depends on some Bioconductor packages. These can be installed using BiocManager package:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicRanges", "BiocParallel", "S4Vectors", "gtools", "stats", "VariantAnnotation"))
```

### Installation

recombClust can be installed from github repository using devtools package:

```r
devtools::install_github("isglobal-brge/recombClust")
```

### Usage

In the folder vignettes, overview.Rmd contains an example on how to apply recombClust to a dataset. 
