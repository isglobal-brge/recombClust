#include "include/CwriteModelHdf5.h"

using namespace Rcpp;

// Only write results does not test nothing because all tests have already been preformed
bool writeResultModel( std::string filename, std::string group,
                       double LoglikeRecomb, double prob0, 
                       NumericVector R1, 
                       std::string grchr, double grstart, double grend )
{
   // Test if file NOT EXISTS :
   //       ==> Create file
   //       ==> Create group
   //       ==> Create datasets with data  (with extension option)
   //       
   // if file EXISTS  
   //       ==> Extend dataset
   //       ==> Add data
   //       ==> ALL DATA IN THE SAME MATRIX ???!!!! (Only one read)
   //             OR CREATE 2 MATRIX ONE WITH MODEL RESULTS AND ANOTHER WITH R1, LOGLIKERECOMB, GSTART AND GEND?? (minimum 2 reads and 2 writes??)
   
   H5File* file;
   DataSet* unlimDataset = nullptr;
   
   // hdf5 parameters
   try{
      
      IntegerVector count = IntegerVector::create(0, 0);
      IntegerVector offset = IntegerVector::create(0, 0);
      IntegerVector stride = IntegerVector::create(1, 1);
      IntegerVector block = IntegerVector::create(1, 1);
      IntegerVector dims_out;
      
      std::string stroutDatasetName = group + "/" + grchr;
      
      H5::Exception::dontPrint();  
      
      // Open or create new file and prepare group and dataset
      file = new H5File( filename, H5F_ACC_RDWR ); 
      
      // extend dataset
      unlimDataset = new DataSet(file->openDataSet(stroutDatasetName));
      dims_out = get_HDF5_dataset_size(*unlimDataset);
      extend_HDF5_matrix_subset_ptr(file, unlimDataset, 1, 0);
      
      Eigen::Map<Eigen::VectorXd> vData(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(R1));
      
      // 4 : LoglikeRecomb + prob0 + gstart + gend
      Eigen::VectorXd otherData(4);
      otherData << LoglikeRecomb, prob0, grstart, grend;
      
      Eigen::VectorXd completeData(vData.size() + otherData.size());
      completeData << vData, otherData;
      
      count[0] = completeData.size();
      count[1] = 1;
      offset[1] = dims_out[1];
      
      write_HDF5_matrix_subset(file, unlimDataset, offset, count, stride, block, Rcpp::wrap(completeData) );  
      
      
      
   } catch( FileIException& error ) { // catch failure caused by the H5File operations
      file->close();
      ::Rf_error( "c++ exception blockmult_hdf5 (File IException)" );
      return wrap(-1);
   } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
      file->close();
      ::Rf_error( "c++ exception blockmult_hdf5 (DataSet IException)" );
      return wrap(-1);
   } catch(std::exception &ex) {
      Rcpp::Rcout<< ex.what();
      return wrap(-1);
   }
   
   
   
   
}