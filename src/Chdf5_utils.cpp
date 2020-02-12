# include "include/Chdf5_utils.h"


/** MOLT IMPORTANT !!!
 *    Documentació bàsica :   https://support.hdfgroup.org/HDF5/examples/intro.html#cpp
 **/

// Check if file exists
inline bool ResFileExist (const std::string& name) {
   struct stat buffer;   
   return (stat (name.c_str(), &buffer) == 0); 
}

/* Get unique values from string vector*/
std::vector<std::string> get_Unique_Values(StringVector vect)
{
   
   std::vector<std::string> v = as<std::vector<std::string> >(vect);
   std::vector<std::string>::iterator ip; 
   
   sort(v.begin(), v.end()); 
   ip = std::unique(v.begin(), v.begin() + v.size());
   
   //Resize vector
   v.resize(std::distance(v.begin(), ip));
   
   return(v); 
}


/* Get first position of each element in a string vector from a vector of unique element*/
std::vector<int> get_Position_Elements( std::vector<std::string> v, StringVector vtocount )
{
   std::vector<int> positions;
   std::vector<std::string> vsearchpos = as<std::vector<std::string> >(vtocount);
   
   for( size_t i=0; i<v.size(); i++)   
   {
      // Find element in vector
      std::vector<std::string>::iterator it = std::find(vsearchpos.begin(), vsearchpos.end(), v[i]);
      positions.push_back( distance(vsearchpos.begin(), it));
   }
   
   return(positions);
}



/* Get the number of elements that exist in an string vector from a vector of unique elements */
std::vector<int> count_Elements_Value( std::vector<std::string> v, StringVector vtocount )
{
   int vtcountlength = as<std::vector<std::string> >(vtocount).size();
   std::vector<std::string> vsearch = as<std::vector<std::string> >(vtocount);
   std::vector<int> counts;
   
   for( size_t i=0; i<v.size(); i++)   {
      //..// Rcpp::Rcout<<"\nSearching : "<<v[i];
      counts.push_back(std::count( vsearch.begin(), vsearch.end(), v[i]));
   }
   
   return(counts);
}



extern "C" {
   
   int create_HiCBrick_dataset(H5std_string filename, const std::string hiCDatasetName, RObject hiCDatasetValues)
   {
      // Try block to detect exceptions raised by any of the calls inside it
      try
      {
         // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
         Exception::dontPrint();
         
         // Open file
         H5File file(filename, H5F_ACC_RDWR);
         
         // Create the data space for the dataset.
         std::vector<int> dims;
         
         if(is<NumericMatrix>(hiCDatasetValues)) 
         {
            hsize_t dims[2];
            dims[0] = as<NumericMatrix>(hiCDatasetValues).rows();
            dims[1] = as<NumericMatrix>(hiCDatasetValues).cols();
            DataSpace dataspace(RANK2, dims);
            
            std::vector<double> matHiCValues = Rcpp::as<std::vector<double> >(hiCDatasetValues);
            /*
            double matHiCValues[dims[0]][dims[1]];
            
            for(int i=0;i<dims[0]; i++)
               for(int j=0; j<dims[1]; j++)
                  matHiCValues[i][j] = as<NumericMatrix>(hiCDatasetValues)(i,j);
            */
            
            DataSet dataset = file.createDataSet(hiCDatasetName, PredType::NATIVE_DOUBLE, dataspace);
            dataset = file.openDataSet(hiCDatasetName);
            // dataset.write( matHiCValues, PredType::NATIVE_DOUBLE);
            dataset.write( &matHiCValues[0], PredType::NATIVE_DOUBLE);
            
            
         } 
         else if(is<IntegerMatrix>(hiCDatasetValues)) 
         {
            hsize_t dims[2];
            dims[0] = as<IntegerMatrix>(hiCDatasetValues).rows();
            dims[1] = as<IntegerMatrix>(hiCDatasetValues).cols();
            DataSpace dataspace(RANK2, dims);
            
            std::vector<int> matHiCValues = Rcpp::as<std::vector<int> >(hiCDatasetValues);
            
            /*
            int matHiCValues[dims[0]][dims[1]];
            
            for(int i=0;i<dims[0]; i++)
               for(int j=0; j<dims[1]; j++)
                  matHiCValues[i][j] = as<IntegerMatrix>(hiCDatasetValues)(i,j);
            */
            
            DataSet dataset = file.createDataSet(hiCDatasetName, PredType::NATIVE_INT, dataspace);
            dataset = file.openDataSet(hiCDatasetName);
            dataset.write( &matHiCValues[0], PredType::NATIVE_INT);
            // dataset.write( matHiCValues, PredType::NATIVE_INT);
         } 
         else if(is<IntegerVector>(hiCDatasetValues) || is<NumericVector>(hiCDatasetValues) 
                    || is<LogicalVector>(hiCDatasetValues) || is<StringVector>(hiCDatasetValues) ) 
         {
            hsize_t vectorsize;
            
            if(is<IntegerVector>(hiCDatasetValues) || is<LogicalVector>(hiCDatasetValues)) 
               vectorsize = as<IntegerVector>(hiCDatasetValues).length();
            else if (is<NumericVector>(hiCDatasetValues))
               vectorsize = as<NumericVector>(hiCDatasetValues).length();
            else if (is<StringVector>(hiCDatasetValues))
               vectorsize = as<StringVector>(hiCDatasetValues).length();
            else 
               vectorsize = 1;
            
            hsize_t dims[] = {vectorsize};
            DataSpace dataspace(RANK1, dims);
            
            
            if(is<IntegerVector>(hiCDatasetValues) || is<LogicalVector>(hiCDatasetValues) ) 
            {
               int vectHiCValues[dims[0]];
               for(int i=0;i<dims[0]; i++)
                  vectHiCValues[i] = as<IntegerVector>(hiCDatasetValues)(i);
               
               DataSet dataset = file.createDataSet(hiCDatasetName, PredType::NATIVE_INT, dataspace);
               dataset = file.openDataSet(hiCDatasetName);
               dataset.write( vectHiCValues, PredType::NATIVE_INT);
            } 
            else if(is<NumericVector>(hiCDatasetValues) ) 
            {
               double vectHiCValues[dims[0]];
               for(int i=0;i<dims[0]; i++)
                  vectHiCValues[i] = as<NumericVector>(hiCDatasetValues)(i);
               
               DataSet dataset = file.createDataSet(hiCDatasetName, PredType::NATIVE_DOUBLE, dataspace);
               dataset = file.openDataSet(hiCDatasetName);
               dataset.write(vectHiCValues, PredType::NATIVE_DOUBLE);
            } 
            else if(is<StringVector>(hiCDatasetValues) ) 
            {
               
               typedef struct name {
                  char chr[MAXSTRING];
               } name;
               
               int datarows = as<StringVector>(hiCDatasetValues).size();
               
               // Convert Dataframe to range list
               name names_list[datarows]; 
               
               for(int i=0; i< datarows; i++ )
               {
                  name n;
                  String wchrom = as<StringVector>(hiCDatasetValues)(i);
                  std::string word = wchrom.get_cstring();
                  
                  int j=0;
                  for( j=0; j < word.size() && j < (MAXSTRING-1); j++ )
                     names_list[i].chr[j] = word[j];
                  
                  names_list[i].chr[j] = '\0'; // insert hdf5 end of string
                  
               }
               
               // Create the memory datatype.
               H5::CompType mtype(sizeof(name));
               mtype.insertMember("chr", HOFFSET(name, chr), H5::StrType(H5::PredType::C_S1, MAXSTRING ));
               
               // Create the dataset.
               DataSet* dataset;
               dataset = new DataSet(file.createDataSet(hiCDatasetName, mtype, dataspace));
               
               dataset->write(names_list, mtype);
               
               // Release resources
               delete dataset;
               
               
            }
         } 
         else if ( is<DataFrame>(hiCDatasetValues))
         {
            
            typedef struct stdata {
               std::string    chr;
               int value;
            } stdata;
            
            // Convert DataFrame to Struct
            size_t  lData = as<DataFrame>(hiCDatasetValues).nrows();
            stdata s1[lData];
            
            StringVector vchr = as<DataFrame>(hiCDatasetValues)[0];
            NumericVector vval = as<DataFrame>(hiCDatasetValues)[1];
            
            for (int i = 0; i< lData; i++) {
               s1[i].chr = vchr(i);
               s1[i].value = vval(i);
            }
            
            // Create the data space.
            hsize_t dim[] = {lData};   /* Dataspace dimensions */
            DataSpace dataspace( RANK1, dim );
            
            // Create the memory datatype.
            CompType dataType( sizeof(stdata) );
            dataType.insertMember( "chr", HOFFSET(stdata, chr), PredType::NATIVE_CHAR);
            dataType.insertMember( "value", HOFFSET(stdata, value), PredType::NATIVE_DOUBLE);
            
            // Create the dataset.
            DataSet* dataset;
            dataset = new DataSet(file.createDataSet(hiCDatasetName, dataType, dataspace));
            
            // Write data to the dataset;
            dataset->write( s1, dataType );

         }
      }  
      catch(FileIException error) { // catch failure caused by the H5File operations
         error.printErrorStack();
         return -1;
      } catch(DataSetIException error) { // catch failure caused by the DataSet operations
         error.printErrorStack();
         return -1;
      } catch(DataSpaceIException error) { // catch failure caused by the DataSpace operations
         error.printErrorStack();
         return -1;
      }
      
      return 0;
      
   }
   
   
   
   /* Creates a hdf5 dataset from dataframe with columns (String column - Numeric Column - Numeric Column) 
    * specific for bintable (ranges) with chromosome (name, start, end)
    */
   int create_HiCBrick_dataset_bintable(H5std_string filename, const std::string hiCDatasetName, RObject hiCDatasetValues)
   {
      
      const std::string member_chr("chr");
      const std::string member_start("start");
      const std::string member_end("end");
      
      /* Structure and dataset */
      typedef struct range {
         char chr[MAXSTRING];
         double start;
         double end;
      } range;
      
      
      StringVector chr;
      NumericVector start, end;
      int datarows = as<DataFrame>(hiCDatasetValues).nrows();
      
      
      if(is<DataFrame>(hiCDatasetValues))
      {
         if( Rf_isFactor(as<DataFrame>(hiCDatasetValues)[0]) )
            chr = as<DataFrame>(hiCDatasetValues)[0];
         else
            throw std::range_error("First column must be a string chromosome name");
         
         if(Rf_isNumeric(as<DataFrame>(hiCDatasetValues)[1]) && Rf_isNumeric(as<DataFrame>(hiCDatasetValues)[2]) )
         {
            start = as<DataFrame>(hiCDatasetValues)[1];
            end = as<DataFrame>(hiCDatasetValues)[2];
         }
         else
            throw std::range_error("Column 2 and 3 must be numerical start and end positions");
      }else
         throw std::range_error("bintable must be a dataframe");
      
      try
      {
         
         Exception::dontPrint();
         
         // Convert Dataframe to range list
         range ranges_list[datarows]; 
         
         for(int i=0; i< datarows; i++ )
         {
            range rng;
            String wchrom = chr(i);
            std::string word = wchrom.get_cstring();
            
            int j=0;
            for( j=0; j < word.size() && j < (MAXSTRING-1); j++ )
               ranges_list[i].chr[j] = word[j];
            
            ranges_list[i].chr[j] = '\0'; // insert hdf5 end of string
            
            ranges_list[i].start = start[i];
            ranges_list[i].end = end[i];
         }
         
         // the array of each length of multidimentional data.
         hsize_t dim[1];
         dim[0] = datarows;
         
         //Preparation of a dataset and a file.
         DataSpace dataspace(RANK1, dim);
         
         // Open file
         H5File* file = new H5File( filename, H5F_ACC_RDWR );
         
         // Create the memory datatype.
         H5::CompType mtype(sizeof(range));
         mtype.insertMember(member_chr, HOFFSET(range, chr), H5::StrType(H5::PredType::C_S1, MAXSTRING ));
         mtype.insertMember(member_start, HOFFSET(range, start), H5::PredType::NATIVE_DOUBLE);
         mtype.insertMember(member_end, HOFFSET(range, end), H5::PredType::NATIVE_DOUBLE);
         
         // Create the dataset.
         DataSet* dataset;
         dataset = new DataSet(file->createDataSet(hiCDatasetName, mtype, dataspace));
         
         dataset->write(ranges_list, mtype);
         
         // Release resources
         delete dataset;
         delete file;
         
      }  
      catch(FileIException error) { // catch failure caused by the H5File operations
         error.printErrorStack();
         return -1;
      } catch(DataSetIException error) { // catch failure caused by the DataSet operations
         error.printErrorStack();
         return -1;
      } catch(DataSpaceIException error) { // catch failure caused by the DataSpace operations
         error.printErrorStack();
         return -1;
      }
      
      return 0;
      
   }
   
   
   /* Create hdf5 dataset from 2 columns Dataframe (String column - Numeric Column) */
   int create_HiCBrick_dataset_CharNum(H5std_string filename, const std::string hiCDatasetName, RObject hiCDatasetValues)
   {
      
      const std::string member_chr("chr");
      const std::string member_double("value");
      
      /* Structure and dataset (range-data) i (start - value)*/ 
      typedef struct data {
         char chr[MAXSTRING];
         double value;
      } data;
      
      
      StringVector chr;
      NumericVector value;
      int datarows = as<DataFrame>(hiCDatasetValues).nrows();
      
      
      if(is<DataFrame>(hiCDatasetValues))
      {
         if( Rf_isFactor(as<DataFrame>(hiCDatasetValues)[0]) )
            chr = as<DataFrame>(hiCDatasetValues)[0];
         else
            throw std::range_error("First column must be a string chromosome name");
         
         if(Rf_isNumeric(as<DataFrame>(hiCDatasetValues)[1]) )
         {
            value = as<DataFrame>(hiCDatasetValues)[1];
         }
         else
            throw std::range_error("Column 2 must be numerical value");
      }else
         throw std::range_error("bintable must be a dataframe");
      
      try
      {
         
         Exception::dontPrint();
         
         // Convert Dataframe to data list
         data data_list[datarows]; 
         
         for(int i=0; i< datarows; i++ )
         {
            data dat;
            String wchrom = chr(i);
            std::string word = wchrom.get_cstring();
            
            int j=0;
            for( j=0; j < word.size() && j < (MAXSTRING-1); j++ )
               data_list[i].chr[j] = word[j];
            
            data_list[i].chr[j] = '\0'; // insert hdf5 end of string
            data_list[i].value = value[i];
         }
         
         // the array of each length of multidimentional data.
         hsize_t dim[1];
         dim[0] = datarows;
         
         //Preparation of a dataset and a file.
         DataSpace dataspace(RANK1, dim);
         
         // Open file
         H5File* file = new H5File( filename, H5F_ACC_RDWR );
         
         // Create the memory datatype.
         H5::CompType mtype(sizeof(data));
         mtype.insertMember(member_chr, HOFFSET(data, chr), H5::StrType(H5::PredType::C_S1, MAXSTRING ));
         mtype.insertMember(member_double, HOFFSET(data, value), H5::PredType::NATIVE_DOUBLE);
         
         // Create the dataset.
         DataSet* dataset;
         dataset = new DataSet(file->createDataSet(hiCDatasetName, mtype, dataspace));
         
         dataset->write(data_list, mtype);
         
         // Release resources
         delete dataset;
         delete file;
         
      }  
      catch(FileIException error) { // catch failure caused by the H5File operations
         error.printErrorStack();
         return -1;
      } catch(DataSetIException error) { // catch failure caused by the DataSet operations
         error.printErrorStack();
         return -1;
      } catch(DataSpaceIException error) { // catch failure caused by the DataSpace operations
         error.printErrorStack();
         return -1;
      }
      
      return 0;
      
   }
   
   
   /* Create a group in hdf5 file */
   int create_HDF5_HiCBrick_group(H5std_string filename, const H5std_string hiCGroup)
   {
      try
      {
         Exception::dontPrint();
         
         // Open file and write group
         H5File file(filename, H5F_ACC_RDWR);
         Group group(file.createGroup(hiCGroup));
         
      } // end of try block
      catch(FileIException error) { // catch failure caused by the H5File operations
         error.printErrorStack();
         return -1;
      } catch(GroupIException error) { // catch failure caused by the Group operations
         error.printErrorStack();
         return -1;
      }
      
      return 0;
   }
   
   
   
   /* Create single attribute to hdf5 group */
   int create_HDF5_HiCBrick_group_attribute(H5std_string HiCfilename, H5std_string HiCObject, 
                                      H5std_string HiCattrName, RObject HiCattr_data)
   {
      
      hsize_t dims[1] = { DIM1 };
      
      try
      {
         Exception::dontPrint();
         
         // Open an existing file and dataset.
         H5File file( HiCfilename, H5F_ACC_RDWR );
         
         // Create the data space for the attribute.
         DataSpace attr_dataspace = DataSpace (1, dims );
         
         Group group = file.openGroup( HiCObject);
         
         if( Rf_isString(HiCattr_data) )
         {
            // Prepare string data
            int strlength = as<std::string>(HiCattr_data).size();
            char stringdata[strlength+1];
            
            std::string word = as<std::string>(HiCattr_data).c_str();
            
            int j=0;
            for( j=0; j < word.size() && j < (strlength); j++ )
               stringdata[j] = word[j];
            
            stringdata[j] = '\0'; // insert hdf5 end of string
            
            // Create group attribute and write
            Attribute attribute = group.createAttribute(HiCattrName, 
                                                        H5::StrType(H5::PredType::C_S1, strlength ), 
                                                        attr_dataspace);
            
            attribute.write(H5::StrType(H5::PredType::C_S1, strlength), stringdata);
         } 
         else if( Rf_isInteger(HiCattr_data) ) 
         {
            int data[] =  {as<int>(HiCattr_data)};
            // Create group attribute and write
            Attribute attribute = group.createAttribute(HiCattrName, PredType::NATIVE_INT,attr_dataspace);
            attribute.write(PredType::NATIVE_INT, data);
         }
         else if( Rf_isNumeric(HiCattr_data) ) 
         {
            double data[] = {as<double>(HiCattr_data)};
            // Create group attribute and write
            Attribute attribute = group.createAttribute(HiCattrName, PredType::NATIVE_DOUBLE,attr_dataspace);
            attribute.write(PredType::NATIVE_DOUBLE, data);
         }
         
      }  // end of try block
      
      catch( DataSpaceIException error ) {
         error.printErrorStack();
         return -1;
      }catch( AttributeIException error ) { // catch failure caused by the H5File operations
         error.printErrorStack();
         return -1;
      }catch( FileIException error ) { // catch failure caused by the H5File operations
         error.printErrorStack();
         return -1;
      }catch( DataSetIException error ) { // catch failure caused by the DataSet operations
         error.printErrorStack();
         return -1;
      }
      
      return 0;  // successfully terminated
      
   }
}
