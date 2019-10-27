#ifndef Cgdsutils
#define Cgdsutils

   #include <Rcpp.h>
   #include <sys/stat.h>
   #include <R_GDS_CPP.h>
   #include <R_GDS.h>
   #include <R_GDS2.h>
   #include <dTrait.h>
   #include <dType.h>
   #include <CoreDEF.h>


   typedef Rcpp::ListOf<Rcpp::NumericMatrix> MatList;
   
   inline bool FileExist (const std::string& name);
   std::string getFileExtension(std::string filePath);
   std::string getFileName(std::string filePath, bool wext, char separator );
   std::string getPath(std::string filePath);
   std::string CGetDatafromFile( std::string file) ;
   Rcpp::NumericMatrix CgdsSNPpairMatrix( Rcpp::StringMatrix filteredsnp) ;
   Rcpp::NumericMatrix CvcfSNPpairMatrix( Rcpp::StringMatrix filteredsnp) ;
   //..// void Transform(Rcpp::NumericVector x);
   Rcpp::IntegerMatrix CTransformtoSampleAlleles(Rcpp::NumericVector  x, Rcpp::StringVector allele, 
                                                 Rcpp::StringVector sample, Rcpp::StringVector variable);

#endif