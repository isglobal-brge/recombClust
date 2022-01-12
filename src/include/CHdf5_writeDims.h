#ifndef CHdf5_writeDims
#define CHdf5_writeDims

    #include<RcppEigen.h>
    #include "H5Cpp.h"
    #include "CHdf5Functions.h"

    // [[Rcpp::depends(RcppEigen)]]
    using namespace H5;
    using namespace Rcpp;


    // C++ funcions
    int write_hdf5_matrix_dimnames(H5File* file, std::string groupname, std::string datasetname, 
                                    StringVector rownames, StringVector colnames );
    
    // R functions
    void write_dimNames(std::string filename, std::string group, std::string dataset,
                            Rcpp::Nullable<Rcpp::StringVector> rownames, 
                            Rcpp::Nullable<Rcpp::StringVector> colnames,
                            Rcpp::Nullable<bool> force  );

#endif