#ifndef CBigLD
#define CBigLD

   #include<RcppEigen.h>
   #include "CUtils_str.h"

   bool in_array(const std::string &value, const std::vector<std::string> &array);

   
   Rcpp::RObject BigLD( Rcpp::RObject corMat, Rcpp::Nullable<double> CLQcut, Rcpp::Nullable<int> clstgap, std::string CLQmode, 
                        std::string hrstType, Rcpp::Nullable<int> hrstParam, std::string chrN);

   

#endif