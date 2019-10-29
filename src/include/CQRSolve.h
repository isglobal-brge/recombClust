#ifndef CQRSolve
#define CQRSolve

   #include<RcppEigen.h>

   Eigen::VectorXd CQR_Solve(Eigen::MatrixXd A, Eigen::VectorXd b);

#endif