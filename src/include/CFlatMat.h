#ifndef CFlatMat
#define CFlatMat

   #include<RcppEigen.h>

   Rcpp::NumericVector flatmatrm(Rcpp::NumericMatrix x);
   Rcpp::NumericVector flatmatcm(Rcpp::NumericMatrix x);

#endif