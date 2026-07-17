#' recombClust: Group chromosomes by recombination history
#'
#' recombClust infers recombination history from a mixture model that, given a
#' pair of SNP-blocks, separates chromosomes into two populations: one with high
#' Linkage Disequilibrium and low recombination (linkage) and another with low LD
#' and high recombination (recomb). This version reimplements the computational
#' core on top of the BigDataStatMeth C++/HDF5 API so the pipeline runs
#' out-of-core.
#'
#' @keywords internal
#' @useDynLib recombClust, .registration = TRUE
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom BigDataStatMeth hdf5_close_all hdf5_create_matrix
## igraph is called from C++ (BigLD graph primitives) via Rcpp::Environment; the
## import makes the dependency explicit and ensures the namespace is loaded.
#' @importFrom igraph graph_from_adjacency_matrix coreness components max_cliques cliques
"_PACKAGE"
