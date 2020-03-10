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
   int iresolution, ilengths, idimensions; 
   std::string strchromosome;
   
   try
   {
      if ( mat.sexp_type()==0 || bintable.sexp_type()==0  )
         throw std::range_error("Data matrix and bintable must exsits and mustn't be null");

      // Create HDF5 file
      H5File file(filename, H5F_ACC_TRUNC);
      
      // Create HiCBrick main groups
      res = create_HDF5_HiCBrick_group(filename, "Base.matrices");
      res = create_HDF5_HiCBrick_group(filename, "Base.ranges");
      res = create_HDF5_HiCBrick_group(filename, "Base.ranges/Bintable");
      res = create_HDF5_HiCBrick_group(filename, "Base.metadata");
      

      StringVector chrom;
      // CharacterVector chromUnique;
      if(is<DataFrame>(bintable))
      {
         // get the different chromosomes to create the subgroups in Base.matrices
         chrom = as<DataFrame>(bintable)[0];
      }else
         throw std::range_error("bintable must be a dataframe");
      
      
      // -- Base.matrices --  Subgroups with matrix Datasets
      
      Eigen::MatrixXd matcor = as<Eigen::MatrixXd>(mat);
      std::vector<int> indexs;
      
      // Test if dim(matcor) == number of bintable rows, if not, add rows and columns until bintable ranges
      if(matcor.rows() < as<DataFrame>(bintable).rows() )
         addRCtoCorrelationMatrix( &matcor, as<DataFrame>(bintable).rows() - matcor.rows());

      std::map<std::string, std::vector<double> > mapChr = MapVectorIndexPosition(chrom);
      std::map<std::string, std::vector<double> >::iterator it = mapChr.begin();
      
      // Get the number of chromosomes in file - only one chromosome is accepted (pro tempora)
      if(mapChr.size()==1 )
      {
         for( std::map<std::string,std::vector<double> >::const_iterator it=mapChr.begin(); it!=mapChr.end(); it++) 
         {
            std::string strsubgroup = "Base.matrices/" + it->first;
            res = create_HDF5_HiCBrick_group(filename, strsubgroup );
            
            strsubgroup = strsubgroup + "/" + it->first;
            res = create_HDF5_HiCBrick_group(filename, strsubgroup );
            
            res = create_HiCBrick_dataset(filename, strsubgroup+ "/matrix", wrap(matcor));
            

            // Base.matrix Datasets
            double intUpperZero = (matcor.array()>0).eval().count();
            double intZero = (matcor.array()==0).eval().count();
            
            //..// NumericVector rbincoverage = wrap(intUpperZero / as<Eigen::MatrixXd>(mat).cols());
            NumericVector rbincoverage = wrap( notZeroinRows(matcor)/matcor.cols() );
            
            NumericVector rowsums = wrap(matcor.rowwise().sum());
            NumericVector sparsity = wrap(intZero / (matcor.cols() * matcor.rows()));
            
            //..// NumericVector cbincoverage = wrap(intUpperZero / as<Eigen::MatrixXd>(mat).rows());
            NumericVector cbincoverage = wrap( notZeroinCols(matcor) / matcor.rows() );
            
            
            NumericVector colsums = wrap(matcor.colwise().sum());
            colsums.attr("dim") = Dimension( colsums.length() ,1);
            
            
            res = create_HiCBrick_dataset(filename, strsubgroup + "/chr1_bin_coverage", rbincoverage);
            res = create_HiCBrick_dataset(filename, strsubgroup + "/chr1_row_sums", rowsums);
            
            res = create_HiCBrick_dataset(filename, strsubgroup + "/chr2_bin_coverage", cbincoverage);
            res = create_HiCBrick_dataset(filename, strsubgroup + "/chr2_col_sums", colsums);
            
            res = create_HiCBrick_dataset(filename, strsubgroup + "/sparsity", sparsity);
            

            // Create atributes
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
            
            // res = create_HDF5_HiCBrick_group_attribute(filename, strsubgroup, "Filename", wrap(filename));
            res = create_HDF5_HiCBrick_group_attribute(filename, strsubgroup, "filename", wrap(getFileName(filename,true,'/')));
            res = create_HDF5_HiCBrick_group_attribute(filename, strsubgroup, "min", 
                                                       wrap(std::to_string(as<Eigen::MatrixXd>(mat).minCoeff())));
            res = create_HDF5_HiCBrick_group_attribute(filename, strsubgroup, "max", 
                                                       wrap(std::to_string(as<Eigen::MatrixXd>(mat).maxCoeff())));
            res = create_HDF5_HiCBrick_group_attribute(filename, strsubgroup, "sparsity", 
                                                       wrap(std::to_string(issparse)));
            res = create_HDF5_HiCBrick_group_attribute(filename, strsubgroup, "distance", 
                                                       wrap(std::to_string(maxdistancediag)));
            res = create_HDF5_HiCBrick_group_attribute(filename, strsubgroup, "done", 
                                                       wrap("1"));
            
            
            
            
         }
      }else {
         throw std::range_error("Only the option of one cromosome is contempled");
         // TODO : Possibility of more than one chromosome
      }

      
            
      // -- Base.metadata --  Subgroups with matrix Datasets

      // Create dataframe with chromosome information
      {
         std::vector<std::string> chr;
         std::vector<double> nrow;
         std::vector<double> size;
         //size_t ielement = 0;
         
         //size_t col=0;
         for( std::map<std::string,std::vector<double> >::const_iterator it=mapChr.begin(); it!=mapChr.end(); it++)
         {
            chr.push_back(it->first );
            nrow.push_back(it->second.size());
            size.push_back(vecmax(as<NumericVector>(as<DataFrame>(bintable)[2]))); // What is size?? - Value is wrong
               // - vecmin(as<NumericVector>(as<DataFrame>(bintable)[2])));
         }

         DataFrame chromoinfo = DataFrame::create( Named("chr") = wrap(chr),
                                                   Named("nrow") = wrap(nrow),
                                                   Named("size") = wrap(size) );
         
         res = create_HiCBrick_dataset_chrominfo(filename, "Base.metadata/chromoinfo", chromoinfo);
         
      }
      
      
      
      // -- Base.ranges --  Subgroups with matrix Datasets
      
      // Create dataframe with chromosome information
      {
         if(is<DataFrame>(bintable))
         {
            if( Rf_isFactor(as<DataFrame>(bintable)[0]) && Rf_isNumeric(as<DataFrame>(bintable)[1]) && Rf_isNumeric(as<DataFrame>(bintable)[2]) )
            {
               std::vector<std::string> uniquechr = get_Unique_Values( as<DataFrame>(bintable)[0] );
               // std::vector<int> noffsets = get_Position_Elements(uniquechr, as<DataFrame>(bintable)[0]);
               std::vector<int> noffsets = get_Position_Elements(uniquechr, as<DataFrame>(bintable)[0], 1);
               std::vector<int> nchr = count_Elements_Value(uniquechr, as<DataFrame>(bintable)[0]);
               
               DataFrame offsets = DataFrame::create( Named("chr") = wrap(uniquechr) , Named("n") = wrap(noffsets) );
               DataFrame lengths = DataFrame::create( Named("chr") = wrap(uniquechr) , Named("n") = wrap(nchr) );
               
               res = create_HiCBrick_dataset_bintable(filename, "Base.ranges/Bintable/ranges", bintable);
               res = create_HiCBrick_dataset(filename, "Base.ranges/Bintable/offset", wrap(noffsets));
               res = create_HiCBrick_dataset(filename, "Base.ranges/Bintable/lengths", wrap(nchr));
               res = create_HiCBrick_String(filename, "Base.ranges/Bintable/chr.names", wrap(uniquechr[0]));
               
               std::vector<int> rsup =  as<DataFrame>(bintable)[2];
               std::vector<int> rinf =  as<DataFrame>(bintable)[1];
               iresolution = ( rsup[0] - rinf[0] ) + 1;
               ilengths = *std::max_element(std::begin(rsup), std::end(rsup));
               strchromosome = uniquechr[0];
               idimensions = nchr[0];
               
            }else
               throw std::range_error("Column 1 must be character, chr name and Column 2 and 3 must be numerical start and end positions");
            
         }else
            throw std::range_error("bintable must be a dataframe");
      }
      
      // Write json configuration file
      std::string strfilename = "HiCBricks_builder_config.json";
      
      Rcpp::List jsonparam = List::create(Named("resolution") = wrap(iresolution),
                                          Named("length") = wrap(ilengths),
                                          Named("dimensions") = wrap(idimensions),
                                          Named("chromosomes") = wrap(strchromosome),
                                          Named("path") = wrap(getPath(filename)),
                                          Named("HDF5File") = wrap( getFileName(filename,true,'/'))
                                          );
      
      Rcpp::Rcout<<"\n"<<getPath(filename) << strfilename<<".json\n";
      Rcpp::Rcout<<"Cromo : "<<strchromosome<<" - dimensions : "<<idimensions<<" - Resolution : "
                 <<iresolution<<" - lengths : "<<ilengths<<"\n";
      res = write_jsonBrick_file( getPath(filename) + strfilename, jsonparam );
      
   }
   catch( FileIException error ) { // catch failure caused by the H5File operations
      error.printErrorStack();
      return -1;
   }
   
   return(0);
}


/*** R

library(HiCBricks)
library(SeqArray)
library(recombClust)

# RecombClust
res <- CRecombClust("inst/extdata/example.vcf", 7, "54301673", "54306217"  )

# direct from results : 
mat.c <- do.call(cbind, lapply(res$models, `[[`, "r1"))
annot.c <- data.frame(t(sapply(res$models, function(x) x$annot)))


# a test.gds file is created the file contains :
# +    [  ]
# \--+ LDMixtureModels   [  ]
#    |--+ models   { Float64 60x153, 71.7K }
#    |--+ annotations   { Float64 153x2, 2.4K }
#     \--+ haplosname   { Str8 60, 360B }
#  
# where :
#     matrix = models
#     annot = annotations
#     
# but it has to be improved


## Only RecombProb without correlation matrix 
X <- cGetRecombProb(as.matrix(mat.c), as.matrix(annot.c), 200)

## Recombprob with correlation matrix
sortidaf <- CgetCorrRecombClust(as.matrix(mat.c), as.matrix(annot.c), 200)

# Creates HiCBricks compatible file
Create_HDF_HiCBricks_File("test/Fitxer_hdf5.brick", sortidaf$corrmat, sortidaf$recombmat$dfchunk)

# HiCBricks functions 
My_BrickContainer <- load_BrickContainer(project_dir = "test/")
My_BrickContainer
Brick_vizart_plot_heatmap(File = "test/HiCBrick_cot.pdf",
                          Bricks = list(My_BrickContainer),
                          x_coords = "chr1:54301673:54306073",
                          y_coords = "chr1:54301673:54306073",
                          resolution = 200,
                          palette = "Reds",
                          rotate = TRUE)

*/
