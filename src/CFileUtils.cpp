#include "include/CFileUtils.h"

using namespace Rcpp;

//' Write data to file
//'
//' @param columns string vector with column names
//' @param results Numerical vector with data
//' @param outputfile relative path to output data file
//' @return none
//' @export
// [[Rcpp::export]]
void CWriteResults( StringVector columns, NumericVector results, std::string outputfile)
{
   std::ofstream outfile;
   // std::string file = "R1_models.txt";
   if( !ResFileExist(outputfile) ) {
      outfile.open(outputfile);
      outfile<<columns<<std::endl;
   }else {
      outfile.open(outputfile, std::ios_base::app);
   }
   outfile<<results<<std::endl;
   
}



// Check if file exists
inline bool ResFileExist (const std::string& name) {
   struct stat buffer;   
   return (stat (name.c_str(), &buffer) == 0); 
}
