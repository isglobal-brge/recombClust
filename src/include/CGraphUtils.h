#ifndef CGraphUtils
#define CGraphUtils

   #include <RcppEigen.h>


   //..// Rcpp::RObject get_local_cores(Rcpp::NumericVector cores, int hrst);
   //..// std::map<int, int> get_local_cores(Rcpp::NumericVector cores, int hrst);
   std::map<int, std::vector<int>> get_local_cores(Rcpp::NumericVector cores, int hrst);
   
   Rcpp::RObject get_graph_matrix_data( Rcpp::RObject OCM, int hrstParam );
   Rcpp::NumericVector get_localcores_map_as_vector(std::map<int, std::vector<int>> localCores);
   

#endif