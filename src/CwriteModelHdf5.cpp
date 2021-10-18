#include "include/CwriteModelHdf5.h"

using namespace Rcpp;

// Only write results does not test nothing because all tests have already been preformed
bool writeResultModel( std::string filename, std::string group,
                       double LoglikeRecomb, double prob0, 
                       NumericVector R1, 
                       std::string grchr, double grstart, double grend )
{

    H5File* file;
    DataSet* unlimDataset = nullptr;
   
    // hdf5 parameters
    try{
      
        IntegerVector count = IntegerVector::create(0, 0);
        IntegerVector offset = IntegerVector::create(0, 0);
        IntegerVector stride = IntegerVector::create(1, 1);
        IntegerVector block = IntegerVector::create(1, 1);
        IntegerVector dims_out;
        bool datasetexists;
        
        std::string stroutDatasetName = group + "/" + grchr;
        
        H5::Exception::dontPrint();  
        
        // Open file
        file = new H5File( filename, H5F_ACC_RDWR ); 
        
        datasetexists = exists_HDF5_element_ptr( file, stroutDatasetName );
        
        count[0] = 1;
        count[1] = R1.size() + 4 ; // + 4 : LoglikeRecomb + prob0 + gstart + gend
        
        // if dataset doesn't exist --> Create unlimited dataset
        if( datasetexists == 0 ) {
            create_HDF5_unlimited_matrix_dataset_ptr(file, stroutDatasetName, count[0], count[1], "numeric");
            dims_out[0] = 0; dims_out[1] = 0;
        } 
        
        // extend dataset (if datasets exists previous execution)
        unlimDataset = new DataSet(file->openDataSet(stroutDatasetName));
        
        // do actions 
        if(datasetexists == 0) {
            CharacterVector cvrowsmames =  R1.names(),
                            cvcolnames;
            cvrowsmames.push_back("LoglikeRecomb");
            cvrowsmames.push_back("prob0");
            cvrowsmames.push_back("gstart");
            cvrowsmames.push_back("gend");
            write_hdf5_matrix_dimnames(file, group, grchr, cvrowsmames, cvcolnames );
        } else {
            dims_out = get_HDF5_dataset_size(*unlimDataset);
        }
        
        if( datasetexists != 0 ) {
            extend_HDF5_matrix_subset_ptr(file, unlimDataset, 1, 0);
        }
        
        Eigen::Map<Eigen::VectorXd> vData(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(R1));
        
        // 4 : LoglikeRecomb + prob0 + gstart + gend
        Eigen::VectorXd otherData(4);
        otherData << LoglikeRecomb, prob0, grstart, grend;
        
        Eigen::VectorXd completeData(vData.size() + otherData.size());
        completeData << vData, otherData;
        
        count[0] = 1;
        count[1] = completeData.size();
        offset[0] = dims_out[0];
        
        write_HDF5_matrix_subset(file, unlimDataset, offset, count, stride, block, Rcpp::wrap(completeData) );  
        
        unlimDataset->close();
        file->close();
      
    } catch( FileIException& error ) { // catch failure caused by the H5File operations
        unlimDataset->close();
        file->close();
        ::Rf_error( "c++ exception blockmult_hdf5 (File IException)" );
        return wrap(-1);
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        unlimDataset->close();
        file->close();
        ::Rf_error( "c++ exception blockmult_hdf5 (DataSet IException)" );
        return wrap(-1);
    } catch(std::exception &ex) {
        unlimDataset->close();
        file->close();
        Rcpp::Rcout<< ex.what();
        return wrap(-1);
    }
    
    return(false);

   
}