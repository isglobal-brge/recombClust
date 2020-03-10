#ifndef CWriteJson
#define CWriteJson

   #include <RcppEigen.h>
   #include <iostream>
   #include <fstream>
   #include <boost/filesystem.hpp>
   
   #include "CFileUtils.h"


   using namespace Rcpp;

   int write_jsonBrick_file(std::string strfilename, RObject params);
   
#endif
