#include "include/CUtils_str.h"
using namespace Rcpp;


// Convert std::vector to char* or std::char<char*>
char *convert(const std::string & s)
{
    char *pc = new char[s.size()+1];
    std::strcpy(pc, s.c_str());
    return pc; 
}


// Returns true if value is in string vector
bool in_array(const std::string &value, const std::vector<std::string> &array)
{
   return std::find(array.begin(), array.end(), value) != array.end();
}