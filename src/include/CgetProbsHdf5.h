#ifndef CgetProbsHdf5
#define CgetProbsHdf5

    #include <RcppEigen.h>
    #include <Rcpp.h>
    #include "CHdf5Functions.h"
    #include "CHdf5Checks.h"

    struct BLOCKS {   // Declare PERSON struct type
        int start;   // Declare member types
        int end;
        int len;
        int startSelection;
        int endESelection;
    };
    
    // C++ functions
    std::vector<BLOCKS> getChunkSelection( std::vector<int>val);

    // R functions
    void getProbs_hdf5( std::string filename, 
                        std::string group, 
                        std::string dataset, 
                        Rcpp::RObject selection, 
                        int nCols,
                        Rcpp::Nullable<std::string> outgroup,
                        Rcpp::Nullable<std::string> outdataset );

#endif