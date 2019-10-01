#ifndef CUpdateModel
#define CUpdateModel

   #include <RcppEigen.h>
   #include <map>
   #include <unordered_map>

   #include "CMapUtils.h"
   #include "CrecombFreq.h"
   #include "CLinkageFreq.h"

   // [[Rcpp::plugins(cpp11)]]
   
   Rcpp::List CupdateModel( Rcpp::List par);



#endif