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
   
   // compute new frequencies for each haplotype in each population (no-inv and inv)
   //..// NumericVector npprops1 = CRecombFreq(inds, blocksize);
   //..// NumericVector npprops2 = CLinkageFreq(inds, blocksize);
   
   NumericVector npprops1 = CRecombFreq( R1, inds, blocksize );
   NumericVector npprops2 = CLinkageFreq( R2, inds, blocksize );
   
   return  List::create(Named("r1") = npr1,  Named("r2") = npr2,  
                        Named("props1") = npprops1, Named("props2") = npprops2,
                        Named("prob0") = npprob0, Named("inds") = inds, Named("blocksize") = blocksize);
}


/***
updateModel <- function(params, funcProps1, funcProps2){
   
   newparams <- params
   
#compute responsabilities
   newparams$r1 <- params$prob0*params$props1[params$inds]
   newparams$r2 <- (1 - params$prob0)*params$props2[params$inds]
   
   R1 <- newparams$r1/(newparams$r1 + newparams$r2)
      R2 <- newparams$r2/(newparams$r1 + newparams$r2)
      
#Compute the new probability of NO inversion;
      newparams$prob0 <- mean(R1)
         
#compute new frequencies for each haplotype in each population (no-inv and inv)
         newparams$props1 <- funcProps1(R1, newparams$inds)
            newparams$props2 <- funcProps2(R2, newparams$inds)
            newparams
}
***/