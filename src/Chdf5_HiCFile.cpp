#include "include/Chdf5_HiCFile.h"


//' Write data to hdf5 file with HiCBricks structure
//'
//' Creates a hdf5 file with HiCBrick structure with mandatory information,
//'   Base.matrices (group)
//'      Chromosome (group)
//'   Base ranges (group)
//'      Bintable (group)
//'
//' and the corresponding information.
//' 
//' @param filename, character array indicating the name of the file to create
//' @param mat numerical mattrix with correlation mattrix
//' @param bintable dataframe with chromosome name
//' \itemize{
//'  \item{"chrom"}{character, chromosome name}
//'  \item{"start"}{numerical, range start position}
//'  \item{"end"}{numerical, range end position}
//' }
//' @return none
//' @export
// [[Rcpp::export]]
int Create_HDF_HiCBricks_File(std::string filename, RObject mat, RObject bintable)
{
   int res;
   
   try
   {
      
      // Create HDF5 file
      H5File file(filename, H5F_ACC_TRUNC);
      
      // Create HiCBrick groups
      res = create_HDF5_HiCBrick_group(filename, "Base.matrices");
      res = create_HDF5_HiCBrick_group(filename, "Base.matrices/chromosome");
      res = create_HDF5_HiCBrick_group(filename, "Base.ranges");
      res = create_HDF5_HiCBrick_group(filename, "Base.ranges/Bintable");
      res = create_HDF5_HiCBrick_group(filename, "Base.metadata");
      
      // Create HiCBrick Datasets
      
      // Base.matrix Datasets
      double intUpperZero = (as<Eigen::MatrixXd>(mat).array()>0).eval().count();
      double intZero = (as<Eigen::MatrixXd>(mat).array()==0).eval().count();
      
      NumericVector bincoverage = wrap(intUpperZero / (as<Eigen::MatrixXd>(mat).cols() * as<Eigen::MatrixXd>(mat).rows()));
      NumericVector rowsums = wrap(as<Eigen::MatrixXd>(mat).rowwise().sum());
      NumericVector sparsity = wrap(intZero / (as<Eigen::MatrixXd>(mat).cols() * as<Eigen::MatrixXd>(mat).rows()));
      
      res = create_HiCBrick_dataset(filename, "Base.matrices/chromosome/matrix", mat);
      res = create_HiCBrick_dataset(filename, "Base.matrices/chromosome/bin.coverage", bincoverage);
      res = create_HiCBrick_dataset(filename, "Base.matrices/chromosome/row.sums", rowsums);
      res = create_HiCBrick_dataset(filename, "Base.matrices/chromosome/sparsity", sparsity);
      
      // Create HiCBrick Attributes
      
      int issparse = 0;
      double maxdistancediag = 0;
      
      // Get data for attributes sparse and maximum distance to diagonal
      {
         
         Eigen::Map<Eigen::MatrixXd> tmat = as<Eigen::Map<Eigen::MatrixXd> >(mat);
         
         if( intZero > (tmat.rows()*tmat.cols()/2) ) issparse = 1;
         
         Eigen::MatrixXd::Index index;
         Eigen::VectorXd v = tmat.diagonal();
         
         // maxdistancediag = (tmat.colwise() - v).colwise().squaredNorm().maxCoeff(&index);
         maxdistancediag = (tmat.colwise() - v).colwise().norm().maxCoeff(&index);
      }
      
      res = create_HDF5_HiCBrick_group_attribute(filename, "Base.matrices", "Filename", wrap(filename));
      res = create_HDF5_HiCBrick_group_attribute(filename, "Base.matrices", "Min", wrap(as<Eigen::MatrixXd>(mat).minCoeff()));
      res = create_HDF5_HiCBrick_group_attribute(filename, "Base.matrices", "Max", wrap(as<Eigen::MatrixXd>(mat).maxCoeff()));
      res = create_HDF5_HiCBrick_group_attribute(filename, "Base.matrices", "sparsity", wrap(issparse));
      res = create_HDF5_HiCBrick_group_attribute(filename, "Base.matrices", "distance", wrap(maxdistancediag));
      res = create_HDF5_HiCBrick_group_attribute(filename, "Base.matrices", "Done", wrap(1));
      
      
      // Base.ranges Datasets (Bintable)
      
      if(is<DataFrame>(bintable))
      {
         if( Rf_isFactor(as<DataFrame>(bintable)[0]) && Rf_isNumeric(as<DataFrame>(bintable)[1]) && Rf_isNumeric(as<DataFrame>(bintable)[2]) )
         {
            std::vector<std::string> uniquechr = get_Unique_Values( as<DataFrame>(bintable)[0] );
            std::vector<int> noffsets = get_Position_Elements(uniquechr, as<DataFrame>(bintable)[0]);
            std::vector<int> nchr = count_Elements_Value(uniquechr, as<DataFrame>(bintable)[0]);
            
            DataFrame offsets = DataFrame::create( Named("chr") = wrap(uniquechr) , Named("n") = wrap(noffsets) );
            DataFrame lengths = DataFrame::create( Named("chr") = wrap(uniquechr) , Named("n") = wrap(nchr) );
            
            res = create_HiCBrick_dataset_bintable(filename, "Base.ranges/Bintable/ranges", bintable);
            res = create_HiCBrick_dataset_CharNum(filename, "Base.ranges/Bintable/offsets", offsets);
            res = create_HiCBrick_dataset_CharNum(filename, "Base.ranges/Bintable/lengths", lengths);
            res = create_HiCBrick_dataset(filename, "Base.ranges/Bintable/chr.names", wrap(uniquechr));
         }else
            throw std::range_error("Column 1 must be character, chr name and Column 2 and 3 must be numerical start and end positions");
         
      }else
         throw std::range_error("bintable must be a dataframe");
      
   }
   catch( FileIException error ) { // catch failure caused by the H5File operations
      error.printErrorStack();
      return -1;
   }
   
   return(0);
}


/*** R





*/
