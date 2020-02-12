#ifndef VectorUtils
#define VectorUtils

   #include <RcppEigen.h>

   using namespace Rcpp;

   NumericVector generate_seq (int start, int end, int inc);
   double vecmin(NumericVector x);
   double vecmax(NumericVector x) ;
   
#endif
