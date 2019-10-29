#ifndef CMapUtils
#define CMapUtils

   #include <RcppEigen.h>
   #include <map>
   #include <unordered_map>

   std::map<std::string, double> VectortoOrderedMap(Rcpp::NumericVector v);
   Rcpp::NumericVector getNumericVectorfromStringVector(std::map<std::string, double> mapv, Rcpp::StringVector strvalues );
   std::unordered_map<std::string, double> sortMapbyFreqs( std::map<std::string, double>sNPFreq );
   std::multimap<double, std::string, std::greater <double>> sortMultimapbyValue( std::map<std::string, double> sNPFreq );

#endif