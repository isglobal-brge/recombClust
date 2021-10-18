#ifndef CHdf5Checks
#define CHdf5Checks

  #include <RcppEigen.h>
  #include "CHdf5Functions.h"
  #include "H5Cpp.h"

  bool exist_FileGroupDataset(std::string filename, std::string group, std::string dataset);
  int exist_File(std::string filename);
  
  double prepare_outDataset(H5File* file, std::string outDataset, bool bforce);
  double prepare_outGroup(H5File* file, std::string outGroup, bool bforce);
  
    
  

#endif