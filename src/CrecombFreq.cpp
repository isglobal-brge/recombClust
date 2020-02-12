#include "include/CrecombFreq.h"
using namespace Rcpp;

//' Reestimate recomb block frequencies
//'
//' @param Resp Numerical with the responsibilities for recomb population
//' @param Block Character with the blocks genotypes
//' @return Numerical with the new block frequencies for recomb population
Rcpp::NumericVector CRecombFreq( NumericVector Resp, Rcpp::StringVector Block, int nSNP )
{
   
   std::map<std::string, double> freqLeft, freqRight;
   std::map<std::string, double> freqs;
   Rcpp::StringVector combs;
   int ncombs = pow(2,nSNP);
   int sizeblock = Block(0).size();
   
   // Get frequencies
   freqs = getFreqSum(Resp, Block, 0, sizeblock);

   // Get combinations and frequency sorted by freq descending 
   freqLeft = getFreqAllCombs(Resp, Block, 0, nSNP);
   freqRight = getFreqAllCombs(Resp, Block, nSNP+1, nSNP);
   
   // solve linear algebra to find new frequency values
   NumericMatrix AA (2*ncombs, 2*ncombs);
   NumericMatrix bb (1, 2*ncombs);
   
   StringVector names(2*ncombs);
   
   AA.fill_diag(1);
   getMatrixfromMap(&AA, &names, &bb, freqLeft, 0);
   getMatrixfromMap(&AA, &names, &bb, freqRight, ncombs);
   
   Eigen::VectorXd props = CQR_Solve( as<Eigen::MatrixXd>(AA), as<Eigen::VectorXd>(bb));

   Eigen::MatrixXd ans = props.head(ncombs) * props.tail(ncombs).adjoint();
   
   // final touch up, round and vector names
   StringVector namesf = concatenate(getNamesbyBlocks(names[ Range(0,ncombs-1)] ), 
                                     getNamesbyWord(names[ Range(ncombs, names.length()-1) ] ), "+");
   NumericVector ansf = flatmatcm( wrap(ans));
   ansf.attr("names") = namesf;
   for(NumericVector::iterator i = ansf.begin(); i != ansf.end(); ++i)
      *i = round( *i * 1000000000000000.0 ) / 1000000000000000.0;


   return ansf;
      
}


