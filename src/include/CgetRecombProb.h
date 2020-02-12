#ifndef CgetRecombProb
#define CgetRecombProb

#include <RcppEigen.h>
#include "VectorUtils.h"
#include "CMapUtils.h"
#include "matrix_utilities.h"

using namespace Rcpp;

// void cGetRecombProb( NumericMatrix probmat, NumericMatrix annot, NumericMatrix range, int window );
List cGetRecombProb( NumericMatrix probmat, NumericMatrix annot, int window );

#endif
