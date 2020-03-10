#ifndef CMapUtils
#define CMapUtils

   #include <RcppEigen.h>
   #include <map>
   #include <unordered_map>

   std::map<std::string, double> VectortoOrderedMap(Rcpp::NumericVector v);
   std::map<double, std::vector<double>> d2MatrixtoOrderedMap(Eigen::MatrixXd m);
   std::map<double, std::vector<int>> VectorstoOrderedMap(Eigen::VectorXd key, Eigen::VectorXd val);
   Rcpp::NumericVector getNumericVectorfromStringVector(std::map<std::string, double> mapv, Rcpp::StringVector strvalues );
   std::unordered_map<std::string, double> sortMapbyFreqs( std::map<std::string, double>sNPFreq );
   std::multimap<double, std::string, std::greater <double>> sortMultimapbyValue( std::map<std::string, double> sNPFreq );
   std::map<std::string, std::vector<double> > MapVectorIndexPosition(Rcpp::StringVector v);
   
#endif