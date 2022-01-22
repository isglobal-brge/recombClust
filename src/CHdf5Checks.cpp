#include "include/CHdf5Checks.h"



bool exist_FileGroupDataset(std::string filename, std::string group, std::string dataset)
{
   
   H5File* file;

   try {
      
      if( ResFileExist_filestream(filename) ) {
         file = new H5File( filename, H5F_ACC_RDWR ); 
      } else {
         Rcpp::Rcout<<"\nFile not exits, create file before split dataset";
         return false;
      }
      
      if( group.compare("") != 0 ) {
         
         if( exists_HDF5_element_ptr(file, group) != 0 ) {
            
            if( dataset.compare("") != 0 ) {
               if( exists_HDF5_element_ptr(file, group + "/" + dataset ) == 0 ) {
                  file->close();
                  Rcpp::Rcout<<"Group not exists, create the input dataset before proceed";
                  return false;
               }
            }
            
         }  else { 
            
            file->close();
            Rcpp::Rcout<<"Group not exists, create the group and dataset before proceed";
            return false;
            
         }
      }
      
   } catch( FileIException& error ) { // catch failure caused by the H5File operations
      file->close();
      ::Rf_error( "c++ exception exist_FileGroupDataset (File IException)" );
      return false;
   } catch( GroupIException& error ) { // catch failure caused by the H5File operations
      file->close();
      ::Rf_error( "c++ exception exist_FileGroupDataset (Group IException)" );
      return false;
   } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
      file->close();
      ::Rf_error( "c++ exception exist_FileGroupDataset (DataSet IException)" );
      return false;
   } catch(std::exception& ex) {
      Rcpp::Rcout<< ex.what();
      return false;
   }
   
   file->close();
   return true;
   
}




int exist_File(std::string filename)
{
    int exists;
    
    try {
        
        if( ResFileExist_filestream(filename) ) {
            exists = 1;
        } else {
            exists = 0;
        }
        
    } catch( FileIException& error ) { // catch failure caused by the H5File operations
        ::Rf_error( "c++ exception exist_File (File IException)" );
        return 3;
    } catch(std::exception& ex) {
        Rcpp::Rcout<< ex.what();
        return 3;
    }
    
    return exists;
    
}





double prepare_outGroup(H5File* file, std::string outGroup, bool bforce)
{
  
  double res = 0;
  
  try {
    
    if( exists_HDF5_element_ptr(file, outGroup) == 0 ) {
      res = create_HDF5_groups_ptr(file, outGroup );
    } /*  Commented to prevent to remove complete group
          maybe there are other interesting datasets inside the group !!!
 
    else if (exists_HDF5_element_ptr(file, outGroup) != 0 && bforce == true)
    {
      res = remove_HDF5_element_ptr(file, outGroup);
      res = create_HDF5_group_ptr(file, outGroup );
    } else {
      throw std::range_error("Group exists, please set force = true if you want to rewrite data");
    }*/
    
    
  } catch( FileIException& error ) { // catch failure caused by the H5File operations
    file->close();
    ::Rf_error( "c++ exception prepare_outGroup (File IException)" );
    return (res);
  } catch( GroupIException& error ) { // catch failure caused by the H5File operations
    file->close();
    ::Rf_error( "c++ exception prepare_outGroup (Group IException)" );
    return (res);
  } catch(std::exception& ex) {
    Rcpp::Rcout<< ex.what();
    return (res);
  }
  
  return(res);
  
  
}


double prepare_outDataset(H5File* file, std::string outDataset, bool bforce)
{
    
    double res = 0;
    
    try {
        
        if( exists_HDF5_element_ptr(file, outDataset) != 0 && bforce == false) {
            Rcpp::Rcout<<"Output dataset exists, please set force = true if you want to rewrite data";
            //..// throw std::range_error("Output dataset exists, please set force = true if you want to rewrite data");
        } else if ( exists_HDF5_element_ptr(file, outDataset) !=0 && bforce == true) {
          
          Rcpp::Rcout<<"\n Hem entrat a la destrucció\n"<< exists_HDF5_element_ptr(file, outDataset)<<"\n Bforce = "<<bforce<<"\n Element"<<outDataset<<"\n";
            res = remove_HDF5_element_ptr(file, outDataset);
        }
        

    } catch( FileIException& error ) { // catch failure caused by the H5File operations
        file->close();
        ::Rf_error( "c++ exception prepare_outGroup (File IException)" );
        return (res);
    } catch( GroupIException& error ) { // catch failure caused by the H5File operations
        file->close();
        ::Rf_error( "c++ exception prepare_outGroup (Group IException)" );
        return (res);
    } catch(std::exception& ex) {
        Rcpp::Rcout<< ex.what();
        return (res);
    }
    
    return(res);
    
    
}



//' Test if element exits
//' 
//' Test if group name or dataset name exists inside the hdf5 data file
//'
//' @param filename string, path and file name to search element.
//' @param element string, group or dataset name to search inside hdf5 data file.
//' @return A boolean indicating whether or not the element being searched exists
//' @export
// [[Rcpp::export]]
bool existsHdf5Element(std::string filename, std::string element)
{
    bool bexists;
    
    H5File* file = new H5File( filename, H5F_ACC_RDWR ); 
    
    if( exists_HDF5_element_ptr( file, element ) == 0 ) {
        bexists =  false;
    } else {
        bexists =  true;
    }
    
    file->close();
    return(bexists);
}



//' Create a group 
//' 
//' Create a group inside Hdf5 data file
//'
//' @param filename string, path and file name to search element.
//' @param element string, group name to be created inside the hdf5 data file.
//' @param overwrite boolean, if true this function overwrites existing group.
//' @return integer, if 0 the process was successful and group was created in 
//' the file
//' @export
// [[Rcpp::export]]
bool setHdf5Group(std::string filename, std::string element, bool overwrite)
{
    int res = 0;
    H5File* file = new H5File( filename, H5F_ACC_RDWR ); 
    
    if( exists_HDF5_element_ptr(file, element) == 0 ) {
        Rcpp::Rcout<<"\nEstem per a crear";
        res = create_HDF5_groups_ptr(file, element );
    } else {
        
        if(overwrite == true) {
            res = remove_HDF5_element_ptr(file, element);
            res = create_HDF5_groups_ptr(file, element );
        } else {
            res = -1;
        }
        
    }
    
    file->close();
    return(res);
}


//' Create File
//' 
//' Create an empty Hdf5 data file
//'
//' @param filename string, path and file name to be created.
//' @param force boolean, if true this function overwrites existing group.
//' @return integer, if 0 the process was successful and group was created in 
//' the file
//' @export
// [[Rcpp::export]]
bool createEmptyHdf5File(std::string filename, Rcpp::Nullable<bool> force = false)
{
    H5File* file; //  = new H5File( filename, H5F_ACC_RDWR ); 
    int res;
    bool overwrite;
    
    if(force.isNull()) { overwrite = false; } 
    else {   overwrite = Rcpp::as<bool>(force); } 
    
    
    res = exist_File(filename);
    
    if( res == 0 ) {
        file = new H5File( filename, H5F_ACC_EXCL ); 
        file->close();
    } else {
      if( overwrite == true) {
        RemoveFile(filename);
        file = new H5File( filename, H5F_ACC_EXCL ); 
        file->close();
      } else {
        Rcpp::Rcout<<"File also exists, please set force = true if you want to overwrite file"; 
      }
    }
    
    return(res);
}
