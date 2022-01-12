#include "include/CHdf5_writeDims.h"



// Write dimnames from Matrix RObjects to hdf5 data file
int write_hdf5_matrix_dimnames(H5File* file, std::string groupname, std::string datasetname, StringVector rownames, StringVector colnames )
{
    
    try{
        
        Exception::dontPrint();
        
        std::string strGroup = groupname + "/." + datasetname + "_dimnames";
        
        // Add rownames
        create_HDF5_groups_ptr(file, strGroup);
        
        if( rownames.length()>1 )
            write_hdf5_string_vector(file, strGroup + "/1" , rownames);
        // else
        //     Rcpp::Rcout<<"Info : no rownames to save";
        
        // Add colnames
        if( colnames.length()>1 )
            write_hdf5_string_vector(file, strGroup + "/2", colnames);
        // else
        //    Rcpp::Rcout<<"Info : no colnames to save";
        
        
    } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
        file->close();
        ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (File IException)" );
        return -1;
    } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
        file->close();
        ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (DataSet IException)" );
        return -1;
    } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
        file->close();
        ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (Group IException)" );
        return -1;
    } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
        file->close();
        ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (DataSpace IException)" );
        return -1;
    } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
        file->close();
        ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (Data TypeIException)" );
        return -1;
    }
    
    return(0);
}




//' Get dataset dimensions
//'
//' get dataset dimensions
//' 
//' @param filename, character array indicating the name of the file to create
//' @param group, Character array, indicating the group where the data set is stored
//' @param datasets, Character array, indicating the datasets name column and/or row names belong
//' @param rownames character vector with row names. 
//' @param colnames character vector with column names. 
//' @return no data returned.
//' @export
// [[Rcpp::export]]
void write_dimNames(std::string filename, std::string group, std::string dataset,
                    Rcpp::Nullable<Rcpp::StringVector> rownames, 
                    Rcpp::Nullable<Rcpp::StringVector> colnames,
                    Rcpp::Nullable<bool> force = false )
{
    
    H5File* file = nullptr;
    std::string finalElement;
    Rcpp::StringVector svrownames, svcolnames, dummy(1);
    
    try
        {
        // int res;
        
        bool bforce;
        
        if(force.isNull()) { bforce = false; } 
        else {   bforce = Rcpp::as<bool>(force); } 
        
        if(!ResFileExist(filename)){
            Rcpp::Rcout<<"File not exits, create file before write dataset row and/or column names ";
            return void();
        }
        
        file = new H5File( filename, H5F_ACC_RDWR );
        
        if( rownames.isNotNull() ){
            std::string strRowNames = group + "/." + dataset + "_dimnames/1";
            if(exists_HDF5_element_ptr(file, strRowNames) && bforce == false) {
                file->close();
                Rcpp::Rcout<<"Rownames allready exits, please set foce = TRUE if you want to overwrite row names";
                return void();
            }
            svrownames = Rcpp::as<Rcpp::StringVector>(rownames);
        } else {
            svrownames = dummy;
        }
        
        if( colnames.isNotNull() ){
            std::string strColNames = group + "/." + dataset + "_dimnames/2";
            if(exists_HDF5_element_ptr(file, strColNames) && bforce == false) {
                file->close();
                Rcpp::Rcout<<"Colnames allready exits, please set foce = TRUE if you want to overwrite column names";
                return void();
            }
            svcolnames = Rcpp::as<Rcpp::StringVector>(colnames);
        } else {
            svcolnames = dummy;
        }
        
        write_hdf5_matrix_dimnames( file, group, dataset, svrownames, svcolnames );
        
    } catch(FileIException& error) { // catch failure caused by the H5File operations
        file->close();
        ::Rf_error( "c++ exception write_dimNames (File IException)" );
        return void();
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
        file->close();
        ::Rf_error( "c++ exception write_dimNames (DataSet IException)" );
        return void();
    } catch(GroupIException& error) { // catch failure caused by the Group operations
        file->close();
        ::Rf_error( "c++ exception write_dimNames (Group IException)" );
        return void();
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
        file->close();
        ::Rf_error( "c++ exception write_dimNames (Data TypeIException)" );
        return void();
    }
    
    file->close();
    
    Rcpp::Rcout<< "Row names and/or column names for "<< dataset <<" has been written to a file\n";  
    return void();
    
}

