#ifndef hdf5_getDim
#define hdf5_getDim

    #include<RcppEigen.h>
    #include "H5Cpp.h"
#include "CHdf5Functions.h"
    
    // [[Rcpp::depends(RcppEigen)]]
    using namespace H5;
    using namespace Rcpp;

    
    // C++ funcions
    Rcpp::IntegerVector get_HDF5_dataset_size(DataSet dataset);
    
    // R functions
    Rcpp::RObject get_dimHdf5(std::string filename, std::string element);

#endif