#include "include/matrix_utilities.h"



Rcpp::NumericVector flatmatrm(Rcpp::NumericMatrix x)
{
   double xcol = x.ncol();
   double xrow = x.nrow();
   double dim = xrow * xcol;
   
   Rcpp::NumericVector fr( dim );
   
   for( double i=0; i< xrow; i++)
   {
      for( double j=0; j<xcol; j++)
      {
         fr[i*xcol+j] = x(i,j);
      }
   }
   
   return(fr);
}


Rcpp::NumericVector flatmatcm(Rcpp::NumericMatrix x)
{
   double xcol = x.ncol();
   double xrow = x.nrow();
   double dim = xrow * xcol;
   
   Rcpp::NumericVector fc( dim );
   
   for( double i=0; i< xrow; i++)
   {
      for( double j=0; j<xcol; j++)
      {
         fc[i+j*xrow] = x(i,j);
      }
   }
   return(fc);
}

/*
 Eigen::MatrixXd RcppNormalize_Data ( Eigen::MatrixXd  X )
 {
 Eigen::RowVectorXd mean = X.colwise().mean();
 Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
 return (X.rowwise() - mean).array().rowwise() / std.array();
 }
 */


Eigen::MatrixXd RcppNormalize_Data ( Eigen::MatrixXd  X, bool bc, bool bs )
{
   Eigen::MatrixXd rX;
   
   if( bc==true && bs==true )  {
      
      Eigen::RowVectorXd mean = X.colwise().mean();
      Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
      rX = (X.rowwise() - mean).array().rowwise() / std.array();
      
   }   else if (bc == true)   {
      
      Eigen::RowVectorXd mean = X.colwise().mean();
      rX = (X.rowwise() - mean);
      
   }  else if ( bs == true)   {
      
      Eigen::RowVectorXd mean = X.colwise().mean();
      Eigen::RowVectorXd std = (X.array().square().colwise().sum() / (X.rows() - 1)).sqrt();
      rX = X.array().rowwise() / std.array();
   } 
   
   return(rX);
}




Rcpp::NumericMatrix RcppNormalize_Data_r ( Rcpp::NumericMatrix  x )
{
   Eigen::MatrixXd X = Rcpp::as<Eigen::MatrixXd> (x);
   Eigen::RowVectorXd mean = X.colwise().mean();
   Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
   return Rcpp::wrap((X.rowwise() - mean).array().rowwise() / std.array());
}



/***R
m <- 10
n <- 10 
x <- matrix(rnorm(m*n), nrow=m, ncol=n)

Normalize_Data(x, bscale = TRUE, bcenter = FALSE )

scale(x, scale = TRUE, center = FALSE)

Dx <- DelayedArray(x)
Normalize_Data(Dx)


*/

