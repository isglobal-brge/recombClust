#include "include/CUpdateModel.h"
using namespace Rcpp;


//' Update model parameters
//'
//' @param params List with current parameter values
//' @param funcProps1 Function used to update block frequencies for population 1
//' @param funcProps2 Function used to update block frequencies for population 2
//' @return List with the updated frequencies
List CupdateModel( List par)
{

   
   
   StringVector inds = as<StringVector>(par["inds"]);
   int blocksize = as<int>(par["blocksize"]);
   
   // compute responsabilities
   NumericVector npr1 = as<double>(par["prob0"]) * getNumericVectorfromStringVector( VectortoOrderedMap( as<NumericVector>(par["props1"]) ), inds );
   NumericVector npr2 = (1-as<double>(par["prob0"])) * getNumericVectorfromStringVector( VectortoOrderedMap( as<NumericVector>(par["props2"]) ), inds );

   NumericVector R1 = npr1/(npr1 + npr2);
   
   NumericVector R2 = npr2/(npr1 + npr2);
   
   // compute the new probability of NO inversion;
   double npprob0 = mean(R1);

   NumericVector npprops1 = CRecombFreq( R1, inds, blocksize );
   NumericVector npprops2 = CLinkageFreq( R2, inds, blocksize );
   
   return  List::create(Named("r1") = npr1,  Named("r2") = npr2,  
                        Named("props1") = npprops1, Named("props2") = npprops2,
                        Named("prob0") = npprob0, Named("inds") = inds, Named("blocksize") = blocksize);
}

