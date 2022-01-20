#ifndef CAnnotationData
#define CAnnotationData

    #include <RcppEigen.h>
    #include "CHdf5Functions.h"
    #include "H5Cpp.h"
    #include "CHdf5Checks.h"

    // R functions
    Rcpp::RObject getAnnotationDataHdf5(std::string resfilename, std::string resgroup, Rcpp::StringVector datasets);

#endif