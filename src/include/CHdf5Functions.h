#ifndef CHdf5Functions
#define CHdf5Functions

   #include <RcppEigen.h>
   #include <fstream>
    #include <sys/stat.h>
   #include <boost/algorithm/string.hpp>
   #include "H5Cpp.h"
    
   
   // [[Rcpp::depends(RcppEigen)]]

   using namespace H5;
   using namespace Rcpp;

   // Constant definitions
   const int	RANK1 = 1;
   const int	RANK2 = 2;
   const int	MAXSTRING = 20;
   const hsize_t MAXSTRBLOCK = 100000;


   // FUNCTION DECLARATION

   // Utils - generics with files
   bool pathExists(hid_t id, const std::string& path);
   bool ResFileExist_filestream(std::string name);
   bool ResFileExist(const std::string& name);


   // HDF5 specific functions (files and groups)

   extern "C" int Create_hdf5_file(std::string filename);
   extern "C" int create_HDF5_group_ptr( H5File* file, const H5std_string mGroup);
   extern "C" int create_HDF5_groups_ptr( H5File* file, const H5std_string mGroup);
   extern "C" int create_HDF5_dataset_ptr(H5File* file, const std::string CDatasetName,
                                         const size_t rows, const size_t cols, std::string strdatatype);
   extern "C" bool exists_HDF5_element_ptr(H5File* file, const H5std_string element);
   IntegerVector get_HDF5_dataset_size(DataSet dataset);
   extern "C" bool remove_HDF5_element_ptr(H5File* file, const H5std_string element);



   // Write data
   extern "C" int write_HDF5_matrix_ptr(H5File* file, const std::string CDatasetName, RObject DatasetValues);

   extern "C" int write_HDF5_matrix_subset( H5File* file, DataSet* dataset,
                                              IntegerVector ivoffset, IntegerVector ivcount,
                                              IntegerVector ivstride, IntegerVector ivblock,
                                              RObject DatasetValues);
   extern "C" int write_hdf5_string_vector(H5File* file, std::string datasetname, StringVector DatasetValues);


   // Get data
   StringVector get_hdf5_matrix_dimnames(H5File* file, std::string groupname, std::string datasetname, int idim );

   extern "C" int read_HDF5_matrix_subset( H5File* file, DataSet* dataset,
                                          IntegerVector ivoffset, IntegerVector ivcount,
                                          IntegerVector ivstride, IntegerVector ivblock,
                                          double* rdatablock);
   extern "C" int extend_HDF5_matrix_subset_ptr(H5File* file, DataSet* dataset, const size_t rows, const size_t cols);
   extern "C" int create_HDF5_unlimited_matrix_dataset_ptr(H5File* file, const std::string CDatasetName, 
                                                          const size_t rows, const size_t cols, 
                                                          std::string strdatatype);

#endif