#include "include/CHdf5Functions.h"


// // Check if file exists
bool ResFileExist_filestream(std::string name) 
{
   
   bool exists = true;
   
   std::fstream fileStream;
   fileStream.open(name);
   
   if (fileStream.fail()) {
      exists = false;
   }
   
   return(exists);
   
}


// Create hdf5 data file with matrix data
int Create_hdf5_file(std::string filename)
{

   try {

      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();

      H5File file(filename, H5F_ACC_RDWR); // Open HDF5 file
      file.close();

   } catch(const H5::FileIException&) {

      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();

      H5File file(filename, H5F_ACC_TRUNC); // Create HDF5 file
      file.close();
   }

   return(0);
}


// Test if file path exists
bool pathExists(hid_t id, const std::string& path)
{
   return H5Lexists( id, path.c_str(), H5P_DEFAULT ) > 0;
}


// Returns if element (group or dataset) exists inside the HDF5 data file
bool exists_HDF5_element_ptr(H5File* file, const H5std_string element)
{
   bool bexists = false;
   try
   {
      H5::Exception::dontPrint();

      // Search dataset
      if(pathExists( file->getId(), element))
         bexists = true;

   } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception exists_HDF5_element_ptr (File IException)" );
      return -1;
   } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception exists_HDF5_element_ptr (DataSet IException)" );
      return -1;
   } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception exists_HDF5_element_ptr (Group IException)" );
      return -1;
   } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception exists_HDF5_element_ptr (DataSpace IException)" );
      return -1;
   } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception exists_HDF5_element_ptr (Data TypeIException)" );
      return -1;
   }
   return bexists;
}


/* Create a group in hdf5 file */
int create_HDF5_group_ptr( H5File* file, const H5std_string mGroup)
{
   try
   {
      Exception::dontPrint();

      std::string strgroup = mGroup;

      // if group no exists -> Create group to file
      if(!pathExists( file->getId(), strgroup))
         file->createGroup("/"+ mGroup);

   } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception create_HDF5_group_ptr (File IException)" );
      return -1;
   } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception create_HDF5_group_ptr (DataSet IException)" );
      return -1;
   } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception create_HDF5_group_ptr (Group IException)" );
      return -1;
   } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception create_HDF5_group_ptr (DataSpace IException)" );
      return -1;
   } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception create_HDF5_group_ptr (Data TypeIException)" );
      return -1;
   }

   return 0;
}

// Create multiple group in hdf5 data file, groups must be separated by "/"
int create_HDF5_groups_ptr( H5File* file, const H5std_string mGroup)
{
   try
   {
      Exception::dontPrint();

      std::string strgroup = mGroup;
      std::string results = "";
      std::vector<std::string> result;

      boost::split(result, mGroup, boost::is_any_of("/"));

      for (int i = 0; i < result.size(); i++) {
         if(!pathExists( file->getId(), results + result[i]))
            file->createGroup(results + result[i]);

         results = result[i] + "/";
      }

   } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception create_HDF5_groups_ptr (File IException)" );
      return -1;
   } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception create_HDF5_groups_ptr (DataSet IException)" );
      return -1;
   } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception create_HDF5_groups_ptr (Group IException)" );
      return -1;
   } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception create_HDF5_groups_ptr (DataSpace IException)" );
      return -1;
   } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception create_HDF5_groups_ptr (Data TypeIException)" );
      return -1;
   }

   return 0;
}



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
      else
         Rcpp::Rcout<<"Info : no rownames to save";

      // Add colnames
      if( colnames.length()>1 )
         write_hdf5_string_vector(file, strGroup + "/2", colnames);
      else
         Rcpp::Rcout<<"Info : no colnames to save";


   } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (File IException)" );
      return -1;
   } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (DataSet IException)" );
      return -1;
   } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (Group IException)" );
      return -1;
   } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (DataSpace IException)" );
      return -1;
   } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (Data TypeIException)" );
      return -1;
   }

   return(0);
}



// Create dataset from stringVector
int write_hdf5_string_vector(H5File* file, std::string datasetname, StringVector DatasetValues)
{

   try
   {

      typedef struct name {
         char chr[MAXSTRING];
      } name;

      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();

      // Create the data space for the dataset.
      hsize_t vectorsize;

      if (is<StringVector>(DatasetValues))
      {
         vectorsize = DatasetValues.length();

         // Define hdf5 dataspace size
         hsize_t dims[] = {vectorsize};
         DataSpace dataspace(RANK1, dims);

         // Create the memory datatype.
         H5::CompType mtype(sizeof(name));
         mtype.insertMember("chr", HOFFSET(name, chr), H5::StrType(H5::PredType::C_S1, MAXSTRING ));

         // Create the dataset.
         DataSet* dataset;
         dataset = new DataSet(file->createDataSet(datasetname, mtype, dataspace));

         // Get dataspace of the dataset.
         dataspace = dataset->getSpace();

         if(vectorsize > MAXSTRBLOCK) {

            // Number of blocks to process
            int iblocsks = vectorsize/MAXSTRBLOCK;

            for(int i=0; i<=iblocsks; i++)
            {

               // Gets block size to read
               hsize_t ilength = MAXSTRBLOCK;
               if(i == iblocsks){
                  ilength = vectorsize - (i * MAXSTRBLOCK);
               }

               // Convert Dataframe to range list
               name *names_list = new name[ilength];

               for(int row=0; row< ilength; row++ )
               {
                  String wchrom = as<StringVector>(DatasetValues)((i*MAXSTRBLOCK) + row);
                  std::string word = wchrom.get_cstring();

                  int j=0;
                  for( j=0; j < word.size() && j < (MAXSTRING-1); j++ ){
                     names_list[row].chr[j] = word[j]; }

                  names_list[row].chr[j] = '\0'; // insert hdf5 end of string
               }

               // HyperSlab position and length
               hsize_t start[1];
               start[0] = (i*MAXSTRBLOCK);
               hsize_t count[] = {ilength};

               DataSpace memspace(RANK1, count, NULL);

               // Get position and write data in dataset
               dataspace.selectHyperslab(H5S_SELECT_SET, count, start);
               dataset->write(names_list, mtype, memspace, dataspace);

               // Release resources
               delete[] names_list;
               memspace.close();
            }

         } else {

            int datarows = as<StringVector>(DatasetValues).size();

            // Convert Dataframe to range list
            name *names_list = new name[datarows];

            for(int i=0; i< datarows; i++ )
            {
               //..// name n;
               String wchrom = as<StringVector>(DatasetValues)(i);
               std::string word = wchrom.get_cstring();

               int j=0;
               for( j=0; j < word.size() && j < (MAXSTRING-1); j++ )
                  names_list[i].chr[j] = word[j];

               names_list[i].chr[j] = '\0'; // insert hdf5 end of string

            }

            dataset->write(names_list, mtype);
            delete[] names_list;
         }

         // Release resources
         dataspace.close();
         dataset->close();
         delete dataset;

      }

   }
   catch(H5::FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception write_hdf5_string_vector (File IException)" );
      return -1;
   } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception write_hdf5_string_vector (DataSet IException)" );
      return -1;
   } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception write_hdf5_string_vector (Group IException)" );
      return -1;
   } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_hdf5_string_vector (DataSpace IException)" );
      return -1;
   }

   return 0;

}


/* Create empty dataset in hdf5 file (pointer to file) */
int create_HDF5_dataset_ptr(H5File* file, const std::string CDatasetName,
                            const size_t rows, const size_t cols, std::string strdatatype)
{

   try
   {
      Exception::dontPrint();

      hsize_t     dimsf[2];              // dataset dimensions
      dimsf[0] = rows;
      dimsf[1] = cols;

      DataSpace dataspace( RANK2, dimsf );

      if( strdatatype == "int") {
         IntType datatype( PredType::NATIVE_INT );
         DataSet dataset = file->createDataSet( CDatasetName, datatype, dataspace );
         dataset.close();
      } else {
         IntType datatype( PredType::NATIVE_DOUBLE );
         DataSet dataset = file->createDataSet( CDatasetName, datatype, dataspace );
         dataset.close();
      }

      dataspace.close();

   } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_dataset_ptr (File IException)" );
      return -1;
   } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_dataset_ptr (DataSet IException)" );
      return -1;
   } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_dataset_ptr (Group IException)" );
      return -1;
   } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_dataset_ptr (DataSpace IException)" );
      return -1;
   } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_dataset_ptr (Data TypeIException)" );
      return -1;
   }

   return 0;
}



int write_HDF5_matrix_ptr(H5File* file, const std::string CDatasetName, RObject DatasetValues)
{
   try
   {
      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();

      // Create the data space for the dataset.
      std::vector<int> dims;

      if(is<NumericMatrix>(DatasetValues))
      {

         hsize_t dims[2];
         dims[0] = as<NumericMatrix>(DatasetValues).rows();
         dims[1] = as<NumericMatrix>(DatasetValues).cols();
         DataSpace dataspace(RANK2, dims);

         std::vector<double> matHiCValues = as<std::vector<double> >(as<NumericMatrix>(DatasetValues));

         DataSet dataset = file->createDataSet(CDatasetName,PredType::NATIVE_DOUBLE, dataspace);
         dataset = file->openDataSet(CDatasetName);

         dataset.write( &matHiCValues[0] , PredType::NATIVE_DOUBLE);

         dataset.close();
         dataspace.close();

      }
      else if( Rcpp::is<IntegerMatrix>(DatasetValues) )
      {
         hsize_t dims[2];
         dims[0] = as<IntegerMatrix>(DatasetValues).rows();
         dims[1] = as<IntegerMatrix>(DatasetValues).cols();
         DataSpace dataspace(RANK2, dims);

         std::vector<double> matHiCValues = as<std::vector<double> >(transpose(as<NumericMatrix>(DatasetValues)));

         DataSet dataset = file->createDataSet(CDatasetName, PredType::NATIVE_DOUBLE, dataspace);
         dataset = file->openDataSet(CDatasetName);
         dataset.write( &matHiCValues[0], PredType::NATIVE_DOUBLE);

         dataset.close();
         dataspace.close();
      } else if(is<NumericVector>(DatasetValues) || is<IntegerVector>(DatasetValues))
      {

         hsize_t vectorsize;

         if(is<IntegerVector>(DatasetValues) || is<LogicalVector>(DatasetValues))
            vectorsize = as<IntegerVector>(DatasetValues).length();
         else if (is<NumericVector>(DatasetValues))
            vectorsize = as<NumericVector>(DatasetValues).length();
         else
            vectorsize = 1;

         hsize_t dims[] = {vectorsize};
         DataSpace dataspace(RANK1, dims);


         if(is<IntegerVector>(DatasetValues) || is<LogicalVector>(DatasetValues) )
         {
            std::vector<int> vectHiCValues(dims[0]);
            //..// int vectHiCValues[dims[0]];
            for(int i=0;i<dims[0]; i++)
               vectHiCValues[i] = as<IntegerVector>(DatasetValues)(i);

            DataSet dataset = file->createDataSet(CDatasetName, PredType::NATIVE_INT, dataspace);
            dataset = file->openDataSet(CDatasetName);
            dataset.write( vectHiCValues.data(), PredType::NATIVE_INT);
            dataspace.close();
            dataset.close();
         }
         else if(is<NumericVector>(DatasetValues) )
         {
            std::vector<double> vectValues(dims[0]);
            //..// double vectValues[dims[0]];
            for(int i=0;i<dims[0]; i++)
               vectValues[i] = as<NumericVector>(DatasetValues)(i);

            DataSet dataset = file->createDataSet(CDatasetName, PredType::NATIVE_DOUBLE, dataspace);
            dataset = file->openDataSet(CDatasetName);
            //..// dataset.write(vectValues, PredType::NATIVE_DOUBLE);
            dataset.write( vectValues.data() , PredType::NATIVE_DOUBLE);
            dataspace.close();
            dataset.close();
         }

      }

   } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception write_HDF5_matrix_ptr (File IException)" );
      return -1;
   } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception write_HDF5_matrix_ptr (DataSet IException)" );
      return -1;
   } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception write_HDF5_matrix_ptr (Group IException)" );
      return -1;
   } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_ptr (DataSpace IException)" );
      return -1;
   } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_ptr (Data TypeIException)" );
      return -1;
   }

   return 0;
}



int write_HDF5_matrix_subset( H5File* file, DataSet* dataset,
                                 IntegerVector ivoffset, IntegerVector ivcount,
                                 IntegerVector ivstride, IntegerVector ivblock,
                                 RObject DatasetValues)
{
   try
   {
      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();

      hsize_t offset[2], count[2], stride[2], block[2];

      // Specify size and shape of subset to write
      offset[0] = ivoffset(0); offset[1] = ivoffset(1);
      count[0]  = ivcount(0); count[1]  = ivcount(1);
      stride[0] = ivstride(0); stride[1] = ivstride(1); // default 1
      block[0] = ivblock(0); block[1] = ivblock(1); // default 1


      // Create the data space for the dataset.
      std::vector<int> dims;
      if(is<NumericMatrix>(DatasetValues) || Rcpp::is<IntegerMatrix>(DatasetValues))
      {

         hsize_t dims[2];
         dims[0] = ivcount[0];
         dims[1] = ivcount[1];
         DataSpace dataspace(RANK2, dims);

         DataSpace memspace(RANK2, dims, NULL);

         dataspace = dataset->getSpace();
         dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);

         // Write a subset of data to the dataset
         std::vector<double> matdata = as<std::vector<double> >(transpose(as<NumericMatrix>(DatasetValues)));

         dataset->write(&matdata[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
         dataspace.close();

      }
      else if(is<NumericVector>(DatasetValues) || is<IntegerVector>(DatasetValues))
      {

         hsize_t dims[2];
         dims[0] = ivcount[0];
         dims[1] = ivcount[1];
         DataSpace dataspace(RANK2, dims);

         DataSpace memspace(RANK2, dims, NULL);

         dataspace = dataset->getSpace();
         dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);

         // Write a subset of data to the dataset
         std::vector<double> matdata = as<std::vector<double> >(as<NumericVector>(DatasetValues));

         dataset->write(&matdata[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
         dataspace.close();
      }

   } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception write_HDF5_matrix_subset_v2 (File IException)" );
      return -1;
   } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception write_HDF5_matrix_subset_v2 (DataSet IException)" );
      return -1;
   } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception write_HDF5_matrix_subset_v2 (Group IException)" );
      return -1;
   } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_subset_v2 (DataSpace IException)" );
      return -1;
   } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_subset_v2 (Data TypeIException)" );
      return -1;
   }

   return 0;
}


// Read dimnames from Hdf5 vector dataset
StringVector get_hdf5_matrix_dimnames(H5File* file, std::string groupname, std::string datasetname, int idim )
{
   StringVector dimnames;

   try{

      Exception::dontPrint();

      std::string strGroup = "." + datasetname + "_dimnames";

      // Get dataset
      std::string strDataset = strGroup + "/" + std::to_string(idim);
      DataSet* dataset = new DataSet(file->openDataSet(strDataset));

      // Define data types
      H5::DataSpace dataspace  = dataset->getSpace();
      H5::StrType   datatype   = dataset->getStrType();

      // Get array dimension !!!
      IntegerVector dims_out = get_HDF5_dataset_size(*dataset);

      for( size_t i = 0; i<dims_out[0]; i++ )
      {
         std::string field_value;
         dataset->read(field_value, datatype, dataspace);
         dimnames.push_back(field_value);

      }

      datatype.close();
      dataspace.close();


   } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception get_hdf5_matrix_dimnames (File IException)" );
      return -1;
   } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception get_hdf5_matrix_dimnames (DataSet IException)" );
      return -1;
   } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception get_hdf5_matrix_dimnames (Group IException)" );
      return -1;
   } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception get_hdf5_matrix_dimnames (DataSpace IException)" );
      return -1;
   } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception get_hdf5_matrix_dimnames (Data TypeIException)" );
      return -1;
   }
   return(dimnames);
}



// Read rhdf5 data matrix subset,
// input :
//      ivoffset : start position
//      ivcount : block size
//      ivstride :(1,1) by default.
//                Allows to sample elements along a dimension. A stride of one (or NULL) will select every element along a dimension,
//                a stride of two will select every other element, and a stride of three will select an element after every two elements.
//      ivblock : (1,1) by default.
//                The block array determines the size of the element block selected from a dataspace. If the block size is one or NULL
//                then the block size is a single element in that dimension.
// output :
//    rdatablock : matrix block
int read_HDF5_matrix_subset(H5File* file, DataSet* dataset,
                            IntegerVector ivoffset, IntegerVector ivcount,
                            IntegerVector ivstride, IntegerVector ivblock,
                            double* rdatablock)
{

   try
   {

      hsize_t offset[2], count[2], stride[2], block[2];
      offset[0] = ivoffset[0]; offset[1] = ivoffset[1];
      count[0] = ivcount[0]; count[1] = ivcount[1];
      stride[0] = ivstride[0]; stride[1] = ivstride[1];
      block[0] = ivblock[0]; block[1] = ivblock[1];

      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();

      // Define Memory Dataspace. Get file dataspace and select a subset from the file dataspace.
      hsize_t dimsm[2];
      dimsm[0] = count[0];
      dimsm[1] = count[1];

      DataSpace memspace(RANK2, dimsm, NULL);

      //  Get dataspace of the dataset.
      DataSpace dataspace = dataset->getSpace();
      dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);

      H5T_class_t type_class = dataset->getTypeClass();

      // Get class of datatype and print message if it's an integer.
      if( type_class == H5T_INTEGER )
         dataset->read( rdatablock, PredType::NATIVE_INT, memspace, dataspace );
      else if (type_class == H5T_FLOAT)
         dataset->read( rdatablock, PredType::NATIVE_DOUBLE, memspace, dataspace );

      dataspace.close();
   } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception read_HDF5_matrix_subset (File IException)" );
      return -1;
   } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception read_HDF5_matrix_subset (DataSet IException)" );
      return -1;
   } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception read_HDF5_matrix_subset (Group IException)" );
      return -1;
   } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception read_HDF5_matrix_subset (DataSpace IException)" );
      return -1;
   } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception read_HDF5_matrix_subset (Data TypeIException)" );
      return -1;
   }

   return 0;  // successfully terminated

}

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

bool remove_HDF5_element_ptr(H5File* file, const H5std_string element)
{
   
   bool bremok = true;
   
   try
   {
      Exception::dontPrint();
      
      int result = H5Ldelete(file->getId(), element.data(), H5P_DEFAULT);
      if(result<0)
         bremok = false;
      
   } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception remove_HDF5_element_ptr (File IException)" );
      return -1;
   } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception remove_HDF5_element_ptr (DataSet IException)" );
      return -1;
   } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception remove_HDF5_element_ptr (Group IException)" );
      return -1;
   } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception remove_HDF5_element_ptr (DataSpace IException)" );
      return -1;
   } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception remove_HDF5_element_ptr (Data TypeIException)" );
      return -1;
   }
   
   return(bremok);
   
}



int extend_HDF5_matrix_subset_ptr(H5File* file, DataSet* dataset, const size_t rows, const size_t cols)
{
   try
   {
      Exception::dontPrint();
      
      int rank, ndims;
      
      // Get dataspace from dataset
      DataSpace dataspace = dataset->getSpace();
      
      // Get the number of dimensions in the dataspace.
      rank = dataspace.getSimpleExtentNdims();
      
      // Get the dimension size of each dimension in the dataspace and
      hsize_t dims_out[2];
      ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
      
      // Create new dataset size from new dims and old dims
      hsize_t   newdims[2];
      newdims[0] = rows;
      newdims[1] = cols;
      
      hsize_t      size[2];
      size[0]   = dims_out[0] + newdims[0];
      size[1]   = dims_out[1] + newdims[1];
      
      dataset->extend( size );
      
   } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception extend_HDF5_matrix_subset_ptr (File IException)" );
      return -1;
   } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception extend_HDF5_matrix_subset_ptr (DataSet IException)" );
      return -1;
   } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception extend_HDF5_matrix_subset_ptr (Group IException)" );
      return -1;
   } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception extend_HDF5_matrix_subset_ptr (DataSpace IException)" );
      return -1;
   } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception extend_HDF5_matrix_subset_ptr (Data TypeIException)" );
      return -1;
   }
   return 0;
   
}