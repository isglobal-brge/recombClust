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
   
   // Rcpp::Rcout<<"\n Dins cUpdateModel";
   
   
   // Rcpp::Rcout<<"\n\bDades amb les que treballarem : \n\n";
   
   // Rcpp::Rcout<<"\n r1 : \n"<<as<NumericVector>(par["r1"])<<"\n";
   // Rcpp::Rcout<<"\n r2 : \n"<<as<NumericVector>(par["r2"])<<"\n";
   // Rcpp::Rcout<<"\n Props1 : \n"<<as<NumericVector>(par["props1"])<<"\n";
   // Rcpp::Rcout<<"\n Props2 : \n"<<as<NumericVector>(par["props2"])<<"\n";
   // Rcpp::Rcout<<"\n inds : \n"<<as<StringVector>(par["inds"])<<"\n";
   // Rcpp::Rcout<<"\n prob0 : \n"<<as<double>(par["prob0"])<<"\n";
   // Rcpp::Rcout<<"\n blocksize : \n"<<as<int>(par["blocksize"])<<"\n";
   
   
   
   StringVector inds = as<StringVector>(par["inds"]);
   int blocksize = as<int>(par["blocksize"]);
   
   // compute responsabilities
   // Rcpp::Rcout<<"\n Calculem responsabilitats npr1";
   NumericVector npr1 = as<double>(par["prob0"]) * getNumericVectorfromStringVector( VectortoOrderedMap( as<NumericVector>(par["props1"]) ), inds );
   // Rcpp::Rcout<<"\n Calculat npr1 -->  calculem npr2";
   NumericVector npr2 = (1-as<double>(par["prob0"])) * getNumericVectorfromStringVector( VectortoOrderedMap( as<NumericVector>(par["props2"]) ), inds );
   // Rcpp::Rcout<<"\n Calculat npr2 -->  calculem R1";
   
   NumericVector R1 = npr1/(npr1 + npr2);
   
   // Rcpp::Rcout<<"\n Calculat R1 -->  calculem R2";
   
   NumericVector R2 = npr2/(npr1 + npr2);
   
   // Rcpp::Rcout<<"\n Calculat R2 -->  calculem npprob0";
   // compute the new probability of NO inversion;
   double npprob0 = mean(R1);
   
   // compute new frequencies for each haplotype in each population (no-inv and inv)
   //..// NumericVector npprops1 = CRecombFreq(inds, blocksize);
   //..// NumericVector npprops2 = CLinkageFreq(inds, blocksize);
   
   // Rcpp::Rcout<<"\n Calculat npprob0 -->  calculem npprops1";
   
   NumericVector npprops1 = CRecombFreq( R1, inds, blocksize );
   
   // Rcpp::Rcout<<"\n Calculat npprops1 -->  calculem npprops2";
   
   NumericVector npprops2 = CLinkageFreq( R2, inds, blocksize );
   
   // Rcpp::Rcout<<"\n Calculat npprops2 -->  Enviem les dades";
   
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