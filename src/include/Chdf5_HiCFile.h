#ifndef Chdf5_HiCFile
#define Chdf5_HiCFile

   #include "Chdf5_utils.h"
   #include <RcppEigen.h>
   #include <iostream>
   #include <string>
   #include "H5Cpp.h"
   #include <sys/stat.h>
   #include "VectorUtils.h"
   #include "CMapUtils.h"
   #include "CWriteJson.h"

   int Create_HDF_HiCBricks_File(std::string filename, RObject mat, RObject bintable);

#endif