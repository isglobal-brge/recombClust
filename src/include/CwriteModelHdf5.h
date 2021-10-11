#ifndef CwriteModelHdf5
#define CwriteModelHdf5

   #include <RcppEigen.h>
   #include "CHdf5Functions.h"
   #include "CHdf5Checks.h"
   
   using namespace Rcpp;

   bool writeResultModel( std::string filename, std::string group,
                          double LoglikeRecomb, double prob0, 
                          NumericVector R1, 
                          std::string grchr, double grstart, double grend );

#endif

