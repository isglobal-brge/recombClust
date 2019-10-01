#include "include/CQRSolve.h"
#include<RcppEigen.h>



Eigen::VectorXd CQR_Solve(Eigen::MatrixXd A, Eigen::VectorXd b)
{
   Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
   return(x);
}
