#include <Rcpp.h>
using namespace Rcpp;

//' timesTwo
//' 
//' This function performs ...
//' 
//' @param x ...
//' @return 2x ...
//' @export
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

/*** R
timesTwo(42)

*/
