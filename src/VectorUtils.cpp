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
   // Rcpp supports STL-style iterators
   NumericVector::iterator it = std::max_element(x.begin(), x.end());
   // we want the value so dereference 
   return *it;
}

