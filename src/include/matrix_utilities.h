#ifndef matrix_utilities
#define matrix_utilities

  #include <RcppEigen.h>

  Rcpp::NumericVector flatmatrm(Rcpp::NumericMatrix x);
  Rcpp::NumericVector flatmatcm(Rcpp::NumericMatrix x);
  Eigen::MatrixXd RcppNormalize_Data ( Eigen::MatrixXd  X, bool bc, bool bs );
  Rcpp::NumericMatrix RcppNormalize_Data_r ( Rcpp::NumericMatrix  x );
  
  

#endif


