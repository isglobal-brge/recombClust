#include "include/CFlatMat.h"
using namespace Rcpp;


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