#include "include/CgetMatrixFromMap.h"
using namespace Rcpp;


// Get Complete Matrix AA from Mapped data.
void getMatrixfromMap( NumericMatrix *AA, StringVector *names, NumericVector *bb, std::map<std::string, double> datamap, int ipos)
{
   
   double imaxf = 999999;
   int i = ipos;
   /*** std::unordered_map<std::string, int>sdatamap = sortMapbyFreqs(datamap); ***/
   std::multimap<double, std::string, std::greater<double>>sdatamap = sortMultimapbyValue(datamap);
   
   bb->at(ipos) = 1;
   
   for (auto x : sdatamap)
   {
      if(i == ipos)  
         imaxf = x.first ;
      
      AA->at(i, ipos) = ( (double)x.first * (-1))/ imaxf ;
      AA->at(ipos,i) = 1;
      
      names->at(i) = x.second;
      
      i++;
   }
   
}