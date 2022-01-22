#ifndef CUtils_str
#define CUtils_str

   #include<RcppEigen.h>

   using namespace Rcpp;

   char *convert(const std::string & s);
   bool in_array(const std::string &value, const std::vector<std::string> &array);


#endif