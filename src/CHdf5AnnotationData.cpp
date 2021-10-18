#include "include/CHdf5AnnotationData.h"


//' Get annotation data from hdf5
//'
//' Get model annotation data from hdf5 data file
//'
//' @param resfilename string, path and file name where models data is store.
//' @param resgroup Annotation obtained for the different models in LDMixture
//' @param datasets Range to look in
//' @return Annotation data from models
//' @export
// [[Rcpp::export]]
Rcpp::RObject getAnnotationDataHdf5(std::string resfilename, std::string resgroup, Rcpp::StringVector datasets)
{
    H5File* file;
    DataSet* pdataset = nullptr;
    Eigen::MatrixXd allAnnotations(0,0);
    
    // hdf5 parameters
    try {
        
        
        IntegerVector dims_out;
        bool datasetexists;
        

        if(exist_File(resfilename) == 1){
            file = new H5File( resfilename, H5F_ACC_RDWR ); 
        } else {
            return(wrap(-1));
        }
        
        // Open different datasets from datasets
        // Seek all datasets to perform calculus
        for( int i=0; i < datasets.size(); i++ ) 
        {
            std::string strdataset = resgroup +"/" + datasets(i);
            
            Rcpp::Rcout<<"\nDataset a llegir : "<<strdataset<<"\n";

            if( exists_HDF5_element_ptr(file, strdataset ) == 0 ) {
                file->close();
                Rcpp::Rcout<<"Group or dataset does not exists, please create the input dataset before proceed";
                return wrap(false);
            }
            
            // Open dataset and get dimension
            pdataset = new DataSet(file->openDataSet(strdataset));
            IntegerVector dims_out = get_HDF5_dataset_size(*pdataset);
            
            // Define positions and sizes
            IntegerVector stride = IntegerVector::create(1, 1);
            IntegerVector block = IntegerVector::create(1, 1);
            
            // 2 : start - end
            IntegerVector count = IntegerVector::create(dims_out[0], 2);
            IntegerVector offset = IntegerVector::create(0, dims_out[1]-2);
            
            Rcpp::Rcout<<"\n Llegirem x per y posicions : "<<count[0]<< "--" <<count[1]<<"\n";
            Rcpp::Rcout<<"\n Començat a la posició x - y : "<<offset[0]<< "--" <<offset[1]<<"\n";
            
            
            // Real data set dimension
            NumericMatrix data( dims_out[0], 2 );
            
            
            // FALLA A LA LECTURA !!! ESTÀ LLEGINT EN ALGUNA POSICIÓ QUE NO TOCA???!!!
            // SEGUR QUE SI QUE ES AIXÒ !!!
            read_HDF5_matrix_subset(file, pdataset, offset, count, stride, block, REAL(data) );
            
            allAnnotations.resize(allAnnotations.rows() + data.nrow() , allAnnotations.cols());
            allAnnotations << allAnnotations, 
                              Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(data);
            
            
        }
            
        // Open file
        // Go to dataset
        // Get dimension dataset
        // Read 2-last columns
        // Result readed data
        
        
        //..// read_HDF5_matrix_subset()

        
        
    } catch( FileIException& error ) { // catch failure caused by the H5File operations
        pdataset->close();
        file->close();
        ::Rf_error( "c++ exception blockmult_hdf5 (File IException)" );
        return wrap(-1);
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        pdataset->close();
        file->close();
        ::Rf_error( "c++ exception blockmult_hdf5 (DataSet IException)" );
        return wrap(-1);
    } catch(std::exception &ex) {
        pdataset->close();
        file->close();
        Rcpp::Rcout<< ex.what();
        return wrap(-1);
    }

    pdataset->close();
    file->close();
    return(wrap(allAnnotations));
    
}