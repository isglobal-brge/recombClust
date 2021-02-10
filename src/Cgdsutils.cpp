#include "include/Cgdsutils.h"

using namespace Rcpp;


// Check if file exists
inline bool FileExist (const std::string& name) {
   struct stat buffer;   
   return (stat (name.c_str(), &buffer) == 0); 
}

/*

// Return file extension
std::string getFileExtension(std::string filePath)
{
   // Find the last position of '.' in given string
   std::size_t pos = filePath.rfind('.');
   
   // If last '.' is found
   if (pos != std::string::npos) 
      return filePath.substr(pos);

   // If no extension return empty string
   return "";
}


// Return file name with or without extension
std::string getFileName(std::string filePath, bool wext, char separator )
{
   // Get dot position (last)
   std::size_t dpos = filePath.rfind('.');
   std::size_t spos = filePath.rfind(separator);
   
   if(spos != std::string::npos)
      return filePath.substr(spos + 1, filePath.size() - (wext || dpos != std::string::npos ? 1 : dpos) );

   return "";
}


// Return file name and path without extension
std::string getPath(std::string filePath)
{
   // Find lenght of filename
   std::size_t pos = getFileName(filePath, true, '/').length();

   // If point of extension is found
   if (pos != std::string::npos)
      return filePath.substr(0, filePath.length() - pos ); // return filename
   
   return filePath;
}

*/


//' Covnert input data file, .vcf or .bed to gds 
//'
//' @param filename, string with route to file with .vcf, .bed or .gds data type
//' @return converted gds filename
//' @export
// [[Rcpp::export]]
std::string CGetDatafromFile( std::string file) 
{
   
   std::string path;
   std::string fileWithoutExt;
   std::string outputfile;
   
   try
   {
   
      if(FileExist (file))
      {
         Rcpp::Rcout<<"\n File exists\n";
         // Get filepath decomposition
         path = getPath(file);
         std::string filename = getFileName(file, false, '/');
         std::string extension =  getFileExtension(file);
         fileWithoutExt = filename.substr(0, filename.length()-extension.length() );

         if( extension !=".bed" && extension != ".vcf" && extension!=".gds")
            throw std::range_error("Unknown file format");

         // Obtaining namespace of SeqArray package
         Environment pkg = Environment::namespace_env("SeqArray");

         if ( extension ==".bed") {
            // Picking up seqBED2GDS() function from SeqArray package
            Function bed2gds = pkg["seqBED2GDS"];
            bed2gds(file, path + fileWithoutExt + ".fam", path + fileWithoutExt + ".bim", path + fileWithoutExt + ".gds");
            outputfile = path + fileWithoutExt + ".gds";
            
         }else if ( extension ==".vcf"){
            // Picking up seqVCF2GDS() function from SeqArray package
            Function vcf2gds = pkg["seqVCF2GDS"];
            vcf2gds(file, path + fileWithoutExt + ".gds");
            outputfile = path + fileWithoutExt + ".gds";
         } else if (extension == ".gds"){
            outputfile = file;
         }
            
         
      } else {
         Rcpp::Rcout<<"\nERROR\n";
         throw std::range_error("File doesn't exist");
      }
   } catch (std::exception &E) {
      forward_exception_to_r(E);
   }
   
   return(outputfile);

}



//' Get's SNP-block pair
//'
//' @param filename, string with route to file with .vcf, .bed or .gds data type
//' @return converted gds filename
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix CgdsSNPpairMatrix( Rcpp::StringMatrix filteredsnp) 
{
   
   int nrows = filteredsnp.nrow(), ncols = filteredsnp.ncol();
   
   Rcpp::NumericMatrix genosnpmatrix( nrows*2, ncols );
   Rcpp::StringVector rn_SNPmats(nrows*2);
   Rcpp::StringVector rnames = rownames(filteredsnp);

   for( int i = 0; i<nrows; i++ )
   {  
      int ipar = 2*i, inon = (ipar)+1;
    
      rn_SNPmats[ipar] = Rcpp::as<std::string>(rnames[i]) + "_1";
      rn_SNPmats[inon] = Rcpp::as<std::string>(rnames[i]) + "_2";

      for( int j = 0; j<ncols; j++ )
      {
         std::string gab = Rcpp::as<std::string>(filteredsnp(i,j));
         
         if( gab.compare("A/A")==0) {
            genosnpmatrix(ipar, j) = 0;
            genosnpmatrix(inon, j) = 0;
         }else if(gab.compare("A/B")==0) {
            genosnpmatrix(ipar, j) = 0;
            genosnpmatrix(inon, j) = 1;
         }else if(gab.compare("B/A")==0) {
            genosnpmatrix(ipar, j) = 1;
            genosnpmatrix(inon, j) = 0;
         }else if (gab.compare("B/B")==0) {
            genosnpmatrix(ipar, j) = 1;
            genosnpmatrix(inon, j) = 1;
         } else {
            genosnpmatrix(ipar, j) = NA_REAL;
            genosnpmatrix(inon, j) = NA_REAL;
         }
      }
   }
   
   rownames(genosnpmatrix) = rn_SNPmats;
   colnames(genosnpmatrix) = colnames(filteredsnp);
   
   return(genosnpmatrix);
   
}



//' Get's SNP-block pair
//'
//' @param filename, string with route to file with .vcf, .bed or .gds data type
//' @return converted gds filename
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix CvcfSNPpairMatrix( Rcpp::StringMatrix filteredsnp) 
{
   
   int nrows = filteredsnp.nrow(), ncols = filteredsnp.ncol();
   
   Rcpp::NumericMatrix genosnpmatrix( nrows*2, ncols );
   Rcpp::StringVector rn_SNPmats(nrows*2);
   Rcpp::StringVector rnames = rownames(filteredsnp);
   
   for( int i = 0; i<nrows; i++ )
   {  
      int ipar = 2*i, inon = (ipar)+1;
      
      rn_SNPmats[ipar] = Rcpp::as<std::string>(rnames[i]) + "_1";
      rn_SNPmats[inon] = Rcpp::as<std::string>(rnames[i]) + "_2";
      
      for( int j = 0; j<ncols; j++ )
      {
         std::string gab = Rcpp::as<std::string>(filteredsnp(i,j));
         
         if( gab.compare("0|0")==0) {
            genosnpmatrix(ipar, j) = 0;
            genosnpmatrix(inon, j) = 0;
         }else if(gab.compare("0|1")==0) {
            genosnpmatrix(ipar, j) = 0;
            genosnpmatrix(inon, j) = 1;
         }else if(gab.compare("1|0")==0) {
            genosnpmatrix(ipar, j) = 1;
            genosnpmatrix(inon, j) = 0;
         }else if (gab.compare("1|1")==0) {
            genosnpmatrix(ipar, j) = 1;
            genosnpmatrix(inon, j) = 1;
         } else {
            genosnpmatrix(ipar, j) = NA_REAL;
            genosnpmatrix(inon, j) = NA_REAL;
         }
      }
   }
   
   rownames(genosnpmatrix) = rn_SNPmats;
   colnames(genosnpmatrix) = colnames(filteredsnp);
   
   return(genosnpmatrix);
   
}



//' Convert 3D-Matrix to 2D-Matrix
//'
//' @param 3D-matrix with allele, samples and SNPs information [allele,samples,SNPs]
//' @return 2D-matrix with one row with each sample-allel [sample-allele, SNP]
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix CTransformtoSampleAlleles(Rcpp::NumericVector  x, Rcpp::StringVector allele, 
                                              Rcpp::StringVector sample, Rcpp::StringVector variable) 
{
   
   Rcpp::NumericVector dim;
   Rcpp::StringVector rnames( sample.length()*2 );

   if(x.hasAttribute("dim"))  {
      dim = x.attr("dim");
   }else { 
      throw std::range_error("Error no dimensions found in array");
   }
   
   if( dim.length()!=3) throw std::range_error("Error array dimensions, expecting 3d-array");
   
   // Number of rows for new matrix
   int rows = dim[0]*dim[1], cols = dim[2];
   
   Rcpp::IntegerMatrix mat2d( rows, cols);

   // Convert 3d-array to 2d-array (matrix)
   for( int i=0; i< dim[0] ; i++)
   {
      for(int j=0; j<dim[1]; j++)
      {
         rnames[2*j+i] = Rcpp::as<std::string>(sample[j]) + Rcpp::as<std::string>(allele[i]);
         for(int k=0; k<dim[2]; k++)
         {
            mat2d(2*j+i,k) = x( 2*j+i+ k*(dim[1]*2));
         }
      }
   }
   
   // rownames and colnames for new matrix
   colnames(mat2d) = variable;
   rownames(mat2d) = rnames;

   return(mat2d);
   
}



/*** R


library(SeqArray)

snpgdsBED2GDS()
gdsfmt::openfn.gds()

file_output_bed <- CGetDatafromFile("inst/extdata/colorectal.bed")




file_output_vcf <- CGetDatafromFile("inst/extdata/example.vcf")




*/
