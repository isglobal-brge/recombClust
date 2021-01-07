#ifndef CCLQD_mod
#define CCLQD_mod

   #include<RcppEigen.h>
   #include"CUtils_str.h"
   #include"CUtils.h"
   #include"CGraphUtils.h"
   #include <Eigen/Dense>



   // C++ Functions (accessible)
   NumericMatrix arrange_matrix_by_cutoff(Eigen::MatrixXd CorMat, double cutoff);
      
   // R functions (accessible)
   Rcpp::RObject CLQD_mod( Eigen::MatrixXd OCM, 
                           double CLQcut, int clstgap, std::string hrstType, 
                           int hrstParam, std::string CLQmode);
   
   
#endif