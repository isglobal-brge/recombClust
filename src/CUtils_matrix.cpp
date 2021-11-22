#include "include/CUtils_matrix.h"

using namespace Rcpp;

//' @title Remove columns and Rows from squared matrix
//' @description removeMatrixColsandRows removes the rows and columns indicated in vector indices vIndex from a matrix
//' @param Mat Matrix to remove index columns and rows
//' @param vIndex Rows and column numbers to be removed from a matrix
//' 
//' <output>
//' @return  The original matrix without rows and column index .
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject removeMatrixColsandRows( Rcpp::RObject Mat, Rcpp::RObject vIndex)
{
   Rcpp::NumericMatrix Original = as<Rcpp::NumericMatrix>(Mat);
   
   
   
   // Create index sequence for all columns
   std::vector<int> v(Original.rows()) ; // vector with 100 ints.
   std::iota (std::begin(v), std::end(v), 0);
   
   
   
   
   
      
   
   
   // Una idea seria generar la seqüencia de 1 a nfiles i ncolumnes
   // un cop creada la sequencia : --> Eliminar els que hi ha a vindex
   //                              --> Fer la sel·lecció + assignació
   // 
   //    Matrix M1;
   //    std::vector<int> rk, ck; // or VectorXi
   //    Matrix M2 = M1(rk,ck);
   
   
   //.. WORKING IN PROGRESS ..// Original = 
   
}