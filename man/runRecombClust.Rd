% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runRecombClust.R
\name{runRecombClust}
\alias{runRecombClust}
\title{Main function to run all the recombClust pipeline}
\usage{
runRecombClust(haplos, annot, clusters = 2, PCs = 3, ...)
}
\arguments{
\item{clusters}{Numeric with the clusters used in k-means}

\item{PCs}{Numeric with the number of PCA components used to make the clustering.}

\item{...}{Further arguments passed to runLDmixtureModel}

\item{filename}{File with genotypic data (SNPs and samples) the data file type must be .vcf (Variant Call Format), 
.gds (Genomic Data Structure) or .bed (Browser Extensible Data).}

\item{start}{start position to be treated}

\item{end}{end position to be treated}
}
\value{
A list with four elements:
\itemize{
 \item{class: Cluster classification of the chromosomes}
 \item{pc: Responsibilities PCA}
 \item{mat: Responsibilities matrix}
 \item{models: List of models}
}
}
\description{
This function takes as input phased chromosomes and returns the responsibilities
PCA and the cluster classification of the samples.
}
