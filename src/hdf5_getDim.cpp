#include "include/hdf5_getDim.h"



IntegerVector get_HDF5_dataset_size(DataSet dataset)
{
    // Get dataspace from dataset
    DataSpace dataspace = dataset.getSpace();
    IntegerVector dims;
    int ndims;
    
    // Get the number of dimensions in the dataspace.
    int rank = dataspace.getSimpleExtentNdims();
    
    // Get the dimension size of each dimension in the dataspace and
    // display them.
    hsize_t dims_out[2];
    ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
    
    if(rank==1)
        dims = IntegerVector::create( static_cast<int>(dims_out[0]), static_cast<int>(1));
    else
        dims = IntegerVector::create(static_cast<int>(dims_out[0]), static_cast<int>(dims_out[1]));
    
    return(dims);
}



//' Get dataset dimensions
//'
//' get dataset dimensions
//' 
//' @param filename, character array indicating the name of the file to create
//' @param element path to element, character array indicating the complete route to the element to query size (folder or dataset). 
//' @return dimension
//' @export
// [[Rcpp::export]]
Rcpp::RObject get_dimHdf5(std::string filename, std::string element)
{
    
    H5File* file = nullptr;
    DataSet dataset;
    IntegerVector dims_out(2), tmp_dims_out;
    
    try
    {
        // int res;
        
        if(!ResFileExist(filename))
            throw std::range_error("File not exits, create file before query datasett");
        
        file = new H5File( filename, H5F_ACC_RDWR );
        
        if(!exists_HDF5_element_ptr(file, element)) {
            file->close();
            throw std::range_error("Element not exits");
        } else{
            dataset = file->openDataSet(element);
            tmp_dims_out = get_HDF5_dataset_size(dataset);  
            dims_out[0] = tmp_dims_out[1];
            dims_out[1] = tmp_dims_out[0];
        }
        
        
    } catch(FileIException& error) { // catch failure caused by the H5File operations
        dataset.close();
        file->close();
        ::Rf_error( "c++ exception bdgetDim (File IException)" );
        return(wrap(-1));
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
        dataset.close();
        file->close();
        ::Rf_error( "c++ exception bdgetDim (DataSet IException)" );
        return(wrap(-1));
    } catch(GroupIException& error) { // catch failure caused by the Group operations
        dataset.close();
        file->close();
        ::Rf_error( "c++ exception bdgetDim (Group IException)" );
        return(wrap(-1));
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
        dataset.close();
        file->close();
        ::Rf_error( "c++ exception bdgetDim (Data TypeIException)" );
        return(wrap(-1));
    }
    
    dataset.close();
    file->close();
    return (dims_out);
    
}
