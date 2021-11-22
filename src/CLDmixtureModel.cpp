#include "include/CLDmixtureModel.h"

using namespace Rcpp;



//' @title LDmixture model to a pair of SNP-blocks
//' 
//' @description Runs LDmixture model to a pair of SNP-blocks
//' 
//' @param dat Matrix with the genotype data
//' @param resfilename string, path and file name to store results. By default
//' this runLDmixture stores data in recombClustResults.hdf5, if resfilename 
//' is empty then function works on memory. To work with big regions it is 
//' recommended to work with data on disk. Only use resfilename = "" when regions 
//' are small in order to prevent memory overflows and low performance.
//' @param resgroup string, folder inside a file where we want to store results. By default
//' runLDmixture stores data inside group "AllModels" in resfilename file.
//' @param overwrite boolean, (optional) either a logical value indicating whether
//' the output file can be overwritten or not, by default files are not overwritten.
//' @param grstart Numerical genomic region start
//' @param grend Numerical genomic region end
//' @param maxSteps Numerical with the maximum number of iterations run by the EM algorithm
//' @param prob0 Initial mixture probability.
//' @return A list with the LDmixture results 
//' \itemize{
//'  \item{"logMix"}{Log-likelihood of mixture model}
//'  \item{"logLD"}{Log-likelihood of linkage model}
//'  \item{"logNoLD"}{Log-likelihood of recomb model}
//'  \item{"BIC"}{BIC of the mixture vs the base model}
//'  \item{"prob"}{Proportion of chromosomes belonging to recomb model}
//'  \item{"steps"}{Number of iterations until converge of the EM algorithm}
//'  \item{"pval"}{P-value of the Chi-square test}
//'  \item{"r1"}{Responsibilities for recomb population of each chromosomes. It is
//'  only available for selected models (BIC > 10, pval > 0.05)}
//' }
//' @export
// [[Rcpp::export]]
List LDmixtureModel( RObject dat, 
                     std::string resfilename, std::string resgroup,
                     std::string grchr, double grstart, double grend,
                     Nullable<int> maxSteps = R_NilValue, 
                     Nullable<double> prob0 = R_NilValue,
                     Nullable<int> blocksize = R_NilValue) 
{
/* dat conté : 
              [,1] [,2]
    HG00100_1 "00" "00"
    HG00100_2 "11" "11"
    HG00101_1 "00" "00"
    HG00101_2 "00" "00"
    HG00102_1 "11" "11"
    HG00102_2 "00" "01"
    ...
*/
   
   int imaxSteps, iblocksize;
   double dprob0;
   CharacterMatrix cdat = transpose(as<CharacterMatrix>(dat));


   // Variable initialization default values
   if(maxSteps.isNotNull())  imaxSteps = as<int> (maxSteps);
   else    imaxSteps = 100;

   if(prob0.isNotNull())  dprob0 = as<double> (prob0);
   else    dprob0 = 0.5;
   
   if(blocksize.isNotNull())  iblocksize = as<int>(blocksize);
   else    iblocksize = 2;

   int nSNP = cdat(0,0).size();
   StringVector cnames = colnames(cdat);
   StringVector inds = concatenate( cdat(0,_), cdat(1,_), "+" );
   
   NumericVector r1(inds.length()); 
   r1.fill(1);

   NumericVector propsRecomb = CRecombFreq(r1, inds, iblocksize);
   NumericVector propsLink = CLinkageFreq(r1, inds, iblocksize);
  
   r1 = getNumericVectorfromStringVector( VectortoOrderedMap(propsRecomb), inds );
   NumericVector r2 = getNumericVectorfromStringVector( VectortoOrderedMap(propsLink), inds );

   r1.attr("names") = inds;
   r2.attr("names") = inds;
   
   double LoglikeRecomb = sum(log(r1));
   double LoglikeLinkage = sum(log(r2));
   
   r1 = dprob0*r1;
   r2 = (1-dprob0)*r2;
   
   
   List params =  List::create(Named("r1") = r1,  Named("r2") = r2,  
                               Named("props1") = propsRecomb, Named("props2") = propsLink,
                               Named("prob0") = dprob0, Named("inds") = inds, Named("blocksize") = iblocksize);

   // EM loop
   double tol =  1;
   int steps = 1;
   double MINTOL = .00001;

   while(tol > MINTOL & steps <= imaxSteps)
   {
      List newparams = CupdateModel(params);
      
      NumericVector parprops1 = as<NumericVector>(params["props1"]);
      NumericVector nparprops1 = as<NumericVector>(newparams["props1"]);
      
      NumericVector parprops2 = as<NumericVector>(params["props2"]);
      NumericVector nparprops2 = as<NumericVector>(newparams["props2"]);
      
      NumericVector difprops1 = (parprops1 - nparprops1);
      NumericVector difprops2 = (parprops2 - nparprops2);

      tol = sqrt( std::inner_product(difprops1.begin(), difprops1.end(), difprops1.begin(), 0.0) +
                  std::inner_product(difprops2.begin(), difprops2.end(), difprops2.begin(), 0.0) +
                  abs(as<double>(params["prob0"]) - as<double>(newparams["prob0"])) );
      
      params = clone(newparams);
      
      steps = steps+1;

      // If one of the populations have all samples, leave the loop
      if (as<double>(newparams["prob0"]) == 1 || as<double>(newparams["prob0"]) == 0) {
         break;
      }
   }

   //get last values to compute likelihood of the complete inversion model
   r1 = as<double>(params["prob0"]) * getNumericVectorfromStringVector( VectortoOrderedMap(as<NumericVector>(params["props1"])), inds );
   r2 = (1 - as<double>(params["prob0"])) * getNumericVectorfromStringVector( VectortoOrderedMap(as<NumericVector>(params["props2"])), inds );
   double LoglikeMix = sum(log(r1 + r2));
   
   NumericVector R1 = r1/(r1 + r2);
   
   R1.attr("names") = inds;
   
   if(resfilename.compare("") == 0)
   {
      return  List::create(   Named("logNoLD") = LoglikeRecomb,  
                              Named("prob") = as<double>(params["prob0"]),
                              Named("r1") = R1,
                              Named("annot") = NumericVector::create(_["start"] = grstart, 
                                    _["end"] = grend)
      );
      
   } else {
      // Write data to disk
      writeResultModel( resfilename, resgroup, LoglikeRecomb, as<double>(params["prob0"]), R1, grchr, grstart, grend );
      
      return(NULL);
   }

   
   
}   




/*** R



library(microbenchmark)
res <- microbenchmark( C <- CrunRecombClust(snpMat[, 1:20], annot = GRsnps[1:20]),
                       R <- runRecombClust(snpMat[, 1:20], annot = GRsnps[1:20]),
                       times = 1L, unit = "s")
print(summary(res)[, c(1:7)],digits=3)

stopifnot(all.equal(C$class, R$class))
stopifnot(all.equal(C$pc, R$pc))
stopifnot(all.equal(C$mat, R$mat))
stopifnot(all.equal(C$models[[5]], R$models[[5]]))


# Prova amb fitxers : 

Cres.file <- CrunRecombClust("inst/extdata/example.vcf",1, 20)

res <- microbenchmark( C.file <- CrunRecombClust("inst/extdata/example.vcf",1, 20),
                       R <- runRecombClust(snpMat[, 1:20], annot = GRsnps[1:20]),
                       times = 1L, unit = "s")

print(summary(res)[, c(1:7)],digits=3)

C.file$class
R$class
C.file$mat[1:5,1:12]
R$mat[1:5,1:12]

*/
