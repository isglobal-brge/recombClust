#include "include/VectorUtils.h"
using namespace Rcpp;


// Aquí també es podria canviar per un push i no cal posar tamany.
NumericVector generate_seq (int start, int end, int inc)
{
   double ilen = std::ceil((end - start)/inc) + 2;
   NumericVector seq(ilen);
   double val = start;
   
   seq[0] = start;
   for (size_t i= 1; val<=end; i++)
   {
      val = val + inc;
      seq[i] = val;
   }
   return(seq);
}

double vecmin(NumericVector x) {
   // Rcpp supports STL-style iterators
   NumericVector::iterator it = std::min_element(x.begin(), x.end());
   // we want the value so dereference 
   return *it;
}


double vecmax(NumericVector x) {

   NumericVector::iterator it = std::max_element(x.begin(), x.end());
   return *it;
}


// Returns unique values from string vector
std::vector<std::string> vectorUniqueStrings(std::vector<std::string> v)
{
   std::sort(v.begin(), v.end()); 
   auto last = std::unique(v.begin(), v.end());
   v.erase(last, v.end());
   return v;
}


// Returns indexes from searched value in vector
std::vector<int> getIndexes(std::vector<std::string> &input, CharacterVector tosearch) 
{
   std::vector<int> v;
   
   for (int i = 0; i < input.size(); i++) {
      if (input[i] == as<std::string>(tosearch)) {
         v.push_back(i);
      }
   }
   
   return v;
}

// Return the number of 0's in rows
Eigen::VectorXd notZeroinRows( Eigen::MatrixXd mat)
{
   Eigen::VectorXd rcounts(mat.rows());
   
   for ( int i; i < mat.rows(); i++ )
      rcounts(i) =  ( mat.row(i).array()!=0.0 ).eval().count();
   
   return rcounts;
}



// Return the number of 0's in columns
Eigen::VectorXd notZeroinCols( Eigen::MatrixXd mat)
{
   Eigen::VectorXd ccounts(mat.cols());
   
   for ( int i; i < mat.cols(); i++ )
      ccounts(i) = ( mat.col(i).array()!=0.0 ).eval().count();
   
   return ccounts;
}
