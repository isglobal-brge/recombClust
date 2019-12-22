#ifndef CFileUtils
#define CFileUtils

   #include <RcppEigen.h>
   #include <sys/stat.h>
   #include <fstream>

   using namespace Rcpp;

   void CWriteResults( StringVector columns, NumericVector results, std::string outputfile);
   inline bool ResFileExist (const std::string& name);

#endif