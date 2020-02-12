#ifndef ClinkageFreq
#define ClinkageFreq

   // [[Rcpp::plugins(cpp11)]]
   
   #include<map>
   #include<unordered_map>
   #include<string>
   
   #include<RcppEigen.h>
   #include"CUtils.h"
   #include"CQRSolve.h"
   // #include"CFlatMat.h"
   #include "matrix_utilities.h"
   #include"CgetMatrixFromMap.h"
   #include"CMapUtils.h"

   Rcpp::NumericVector CLinkageFreq( Rcpp::NumericVector Resp, Rcpp::StringVector Block, int nSNP );

#endif