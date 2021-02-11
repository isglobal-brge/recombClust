#include "include/CgetCorrelationMatrix.h"

using namespace Rcpp;


//' Compute correlation matrix
//' 
//' Compute correlation matrix from probabilitie matrix. 
//'
//' @param matrix, Numerical matrix with probabilites
//' @return matrix, Correlation matrix
//' @export
// [[Rcpp::export]]
Rcpp::RObject getCorrelationMatrix(Rcpp::RObject mat,  Rcpp::Nullable<bool> absval = R_NilValue) 
{
   
   Eigen::MatrixXd nmat = Rcpp::as<Eigen::MatrixXd>(mat);
   Eigen::MatrixXd cov;
   bool babsval = false;
   
   if( absval.isNull()) {
      babsval = false;
   } else {
      babsval = Rcpp::as<bool> (absval);
   }
   
   Eigen::MatrixXd centered = nmat.rowwise() - nmat.colwise().mean();
   
   if( nmat.rows()>1) {
      cov = (centered.adjoint() * centered) / double(nmat.rows() - 1);
   }
   
   if( babsval == true) {
      cov = cov.cwiseAbs();
   }
   
   return(wrap(cov));
}


