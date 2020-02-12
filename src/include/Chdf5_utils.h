#ifndef Chdf5_utils
#define Chdf5_utils


   #include <RcppEigen.h>
   #include <iostream>
   #include <string>
   #include "H5Cpp.h"
   #include <sys/stat.h>
   
   using namespace H5;
   using namespace Rcpp;
   
   const int	 RANK1 = 1;
   const int	 RANK2 = 2;
   const int	DIM1 = 1;
   const int	DIM2 = 2;
   const int	MAXSTRING = 32;
   
   
   
   inline bool ResFileExist (const std::string& name);
   std::vector<std::string> get_Unique_Values(StringVector vect);
   std::vector<int> get_Position_Elements( std::vector<std::string> v, StringVector vtocount );
   std::vector<int> count_Elements_Value( std::vector<std::string> v, StringVector vtocount );
   
   
   extern "C" int create_HiCBrick_dataset(H5std_string filename, const std::string hiCDatasetName, RObject hiCDatasetValues);
   extern "C" int create_HiCBrick_dataset_bintable(H5std_string filename, const std::string hiCDatasetName, RObject hiCDatasetValues);
   extern "C" int create_HiCBrick_dataset_CharNum(H5std_string filename, const std::string hiCDatasetName, RObject hiCDatasetValues);
   extern "C" int create_HDF5_HiCBrick_group(H5std_string filename, const H5std_string hiCGroup);
   extern "C" int create_HDF5_HiCBrick_group_attribute(H5std_string HiCfilename, H5std_string HiCObject, 
                                                      H5std_string HiCattrName, RObject HiCattr_data);


#endif