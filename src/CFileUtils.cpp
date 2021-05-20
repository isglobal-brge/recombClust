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
   
   // Remove path from filePath
   filePath = filePath.substr(getPath(filePath).length(), filePath.length());

   // Get dot position (last)
   std::size_t dpos = filePath.rfind('.');
   std::size_t spos = filePath.rfind(separator);
   
   if(spos != std::string::npos)
      return filePath.substr(spos + 1, filePath.size() - (wext || dpos != std::string::npos ? 1 : dpos) );
   else
      return filePath.substr(0, filePath.size() - (wext || dpos != std::string::npos ? 1 : dpos) );
   
   return "";
}


// // Return file name and path without extension
// std::string getPath(std::string filePath)
// {
//    // Find lenght of filename
//    std::size_t pos = getFileName(filePath, true, '/').length();
//    
//    // If point of extension is found
//    if (pos != std::string::npos)
//       return filePath.substr(0, filePath.length() - pos ); // return filename
//    
//    return filePath;
// }



// Return file name and path without extension
std::string getPath(std::string filePath)
{
   char sep = '/';
   
#ifdef _WIN32
   sep = '\\';
#endif
   
   size_t i = filePath.rfind(sep, filePath.length());
   if (i != std::string::npos) {
      return(filePath.substr(0, i));
   }
   
   return("");
}