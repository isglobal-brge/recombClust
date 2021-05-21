# recombClust

## Description

recombClust is a R package that groups chromosomes by their recombination history. Recombination history is based on a mixture model that, given a pair of SNP-blocks, separates chromosomes in two populations, one with high Linkage Disequilibrium (LD)  and low recombination (linkage) and another with low LD and high recombination. The method use the classification of several SNP-block pairs in a region to group chromosomes in clusters with different recombination history. This package takes as input genotype phased data.

The methods described in the paper **"Identifying chromosomal subpopulations based on their recombination histories advances the study of the genetic basis of phenotypic traits" by Carlos Ruiz-Arenas et al. (2020)**[[1]](#1) available at **[Genome Research](https://genome.cshlp.org/content/early/2020/11/17/gr.258301.119)** can be found in release [v 1.0.0](https://github.com/isglobal-brge/recombClust/releases/tag/v1.0.0).

Release [v 0.0.9](https://github.com/isglobal-brge/recombClust/releases/tag/v0.0.9) contains the methods described in preprint **"Historical recombination variability contributes to deciphering the genetic basis of phenotypic traits" by Carlos Ruiz-Arenas et al. (2019)**[[2]](#2) publicly available at [bioRxiv](https://www.biorxiv.org/content/10.1101/792747v1)



### Pre-requisites

recombClust depends on some Bioconductor packages. These can be installed using BiocManager package:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicRanges", "BiocParallel", "S4Vectors", "gtools", "stats", "VariantAnnotation", "SeqVarTools"))
```

### Installation

recombClust can be installed from github repository using devtools package:

```r
devtools::install_github("isglobal-brge/recombClust")
```

### Usage

In the folder vignettes, overview.Rmd contains an example on how to apply recombClust to a dataset. 


### References
<a id="1">[1]</a> 
Ruiz-Arenas C, Cáceres A, López M, Pelegrí-Sisó D, González J, González JR. Identifying chromosomal subpopulations based on their recombination histories advances the study of the genetic basis of phenotypic traits. Genome Res. 2020 Nov 17. doi: [10.1101/gr.258301.119](https://doi.org/10.1101/gr.258301.119). Epub ahead of print. PMID: 33203765.
\  

<a id="2">[2]</a> 
Historical recombination variability contributes to deciphering the genetic basis of phenotypic traits
Carlos Ruiz-Arenas, Alejandro Cáceres, Marcos López, Dolors Pelegrí-Sisó, Josefa González, Juan R. González
bioRxiv 792747; doi: [https://doi.org/10.1101/792747](https://doi.org/10.1101/792747)