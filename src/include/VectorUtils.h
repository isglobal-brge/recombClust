#ifndef VectorUtils
#define VectorUtils

   #include <RcppEigen.h>

   using namespace Rcpp;

   NumericVector generate_seq (int start, int end, int inc);
   double vecmin(NumericVector x);
   double vecmax(NumericVector x);
   std::vector<std::string> vectorUniqueStrings(std::vector<std::string> v);
   std::vector<int> getIndexes(std::vector<std::string> &input, CharacterVector searched);

   // Matrix Related
   Eigen::VectorXd notZeroinCols( Eigen::MatrixXd mat);
   Eigen::VectorXd notZeroinRows( Eigen::MatrixXd mat);
      
#endif
