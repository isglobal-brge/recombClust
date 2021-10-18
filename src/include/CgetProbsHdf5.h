#ifndef CgetProbsHdf5
#define CgetProbsHdf5

    #include <RcppEigen.h>
    #include <Rcpp.h>
    #include "CHdf5Functions.h"
    #include "CHdf5Checks.h"

    Rcpp::RObject getProbs_hdf5( std::string filename, 
                                 std::string group, 
                                 Rcpp::RObject selection, 
                                 Rcpp::Nullable<std::string> outgroup,
                                 Rcpp::Nullable<std::string> outdataset);

#endif