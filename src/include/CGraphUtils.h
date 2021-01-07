#ifndef CGraphUtils
#define CGraphUtils

   #include <RcppEigen.h>


   Rcpp::RObject get_local_cores(Rcpp::NumericVector cores, int hrst);

   Rcpp::RObject get_graph_matrix_data( Rcpp::RObject OCM, int hrstParam );

#endif