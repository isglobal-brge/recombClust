#include "include/CUtils_str.h"
using namespace Rcpp;


// Returns true if value is in string vector
bool in_array(const std::string &value, const std::vector<std::string> &array)
{
   return std::find(array.begin(), array.end(), value) != array.end();
}