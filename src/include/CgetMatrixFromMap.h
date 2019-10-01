#ifndef CgetMatrixFromMap
#define CgetMatrixFromMap

   #include<map>
   #include<unordered_map>
   #include<string>
      
   #include<RcppEigen.h>
   #include"CUtils.h"
   #include"CMapUtils.h"
      
   void getMatrixfromMap( Rcpp::NumericMatrix *AA, Rcpp::StringVector *names, 
                          Rcpp::NumericVector *bb, std::map<std::string, double> datamap, int ipos);
   
#endif