#include "include/CgetProbs.h"

using namespace Rcpp;

//' Compute cluster Recomb freq by mean of voting
//'
//' @param mat with initial probabilities
//' @param sel overlap ranges
//' @return NumericVector with recombination probabilites
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector getProbs( Rcpp::RObject mat, Rcpp::RObject sel )
{
   Eigen::MatrixXd dmat = Rcpp::as<Eigen::MatrixXd>(mat);
   Eigen::VectorXd dsel = Rcpp::as<Eigen::VectorXd>(sel);
   
   Eigen::MatrixXd dmatsel = Eigen::MatrixXd::Zero(dmat.rows(), dsel.size());

   for( int i=0; i<dsel.size(); i++){
      dmatsel.col(i) = dmat.col(dsel[i]-1);
   }

   return(wrap(dmatsel.rowwise().mean()));
  
}
