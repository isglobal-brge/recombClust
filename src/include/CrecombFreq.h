#ifndef CrecombFreq
#define CrecombFreq

// [[Rcpp::plugins(cpp11)]]

   #include<map>
   #include<unordered_map>
   #include<string>

   #include<RcppEigen.h>
   #include"CUtils.h"
   #include"CQRSolve.h"
   #include"CFlatMat.h"
   #include"CgetMatrixFromMap.h"
   

   // dgemv_ is a symbol in the LAPACK-BLAS - Level 3
   //    DGEQP3  DGEQP3 computes a QR factorization with column pivoting of a matrix A:  A*P = Q*R  using Level 3 BLAS.
   extern "C" {
      extern void dgeqp3_(int*, int*, double*, int*, int*, double*, double*, int*, int*);
   }

   Rcpp::NumericVector CRecombFreq( Rcpp::NumericVector Resp, Rcpp::StringVector Block, int nSNP );
   // void getMatrixfromMap( Rcpp::NumericMatrix *AA, Rcpp::StringVector *names,
   //                       Rcpp::NumericVector *bb, std::map<std::string, int> datamap, int ipos);
   //
#endif