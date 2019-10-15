#include "include/CLinkageFreq.h"

using namespace Rcpp;


//' Reestimate linkage block frequencies
//'
//' @param Resp Numerical with the responsibilities for linkage population
//' @param Block Character with the blocks genotypes
//' @return Numerical with the new block frequencies for linkage population
Rcpp::NumericVector CLinkageFreq( Rcpp::NumericVector Resp, Rcpp::StringVector Block, int nSNP )
{
   int ncombs = pow( 2, nSNP );
   int sizeblock = Block(0).size();
   
   // Get frequencies
   std::map<std::string, double> freqs = getFreqSum(Resp, Block, 0, sizeblock);
   std::multimap<double, std::string, std::greater<double>>sfreqs = sortMultimapbyValue(freqs);

   // Make all combination with 0 and 1
   int arr[nSNP];

   Rcpp::StringVector leftlevs, rightlevs, leftLevsvec, rightLevsvec;
   std::map<std::string, double>selAlleles;

   getAllBinaryStrings(nSNP, arr, 0, &leftlevs);
   rightLevsvec = leftlevs;
   rightlevs = rightLevsvec; 
   leftLevsvec = leftlevs;

   std::unordered_map<std::string,double>::iterator it;
   std::multimap<double, std::string>::iterator itm;

   
   while( sfreqs.size()>0 && leftLevsvec.size()>0 )
   {

      itm = sfreqs.begin();
      std::string allele = itm->second;
      
      if ( (std::find(leftLevsvec.begin(), leftLevsvec.end(), allele.substr(0,nSNP)) != leftLevsvec.end())  &&
           (std::find(rightLevsvec.begin(), rightLevsvec.end(), allele.substr(nSNP+1,nSNP)) != rightLevsvec.end()) )
      {

         selAlleles.insert( selAlleles.end(), std::pair<std::string, double>(allele, 0));
         
         auto itl = std::find(leftLevsvec.begin(), leftLevsvec.end(), allele.substr(0,nSNP));
         if (itl != leftLevsvec.end()) leftLevsvec.erase(itl);
         
         auto itr = std::find(rightLevsvec.begin(), rightLevsvec.end(), allele.substr(nSNP+1,nSNP));
         if (itr != rightLevsvec.end()) rightLevsvec.erase(itr);
         
      }
      sfreqs.erase(itm);
      
   }
   
   if(leftLevsvec.length()>0) {
      StringVector leftright = concatenate(leftLevsvec, rightLevsvec, "+");

      for( int i=0; i<leftright.length() ; i++)
         selAlleles.insert( selAlleles.end(), std::pair<std::string, double>((std::string)leftright[i], 0));
   }

   for (std::pair<std::string, double> x : selAlleles){
      std::map<std::string,double>::iterator ito;
      ito = freqs.find( x.first );
      if (ito != freqs.end())  
         selAlleles[x.first] = ito->second ;
   }
   

   // solve linear algebra to find new frequency values
   NumericMatrix AA (ncombs, ncombs);
   NumericMatrix bb (1,ncombs);
   int fullcombs = pow(2,2*nSNP);
   
   StringVector names(ncombs);
   
   AA.fill_diag(1);
   getMatrixfromMap(&AA, &names, &bb, selAlleles, 0);
   

   Eigen::VectorXd props = Eigen::VectorXd::Zero(fullcombs);
   props.segment(0,ncombs) = CQR_Solve( as<Eigen::MatrixXd>(AA), as<Eigen::VectorXd>(bb));
   props = props.unaryExpr([](double v) { return v>1e-5? v : 1e-5; });
   props = props /props.array().sum();
   
   StringVector namesf;
   
   Rcpp::StringVector combs;
   getAllBinaryStrings( 2*nSNP, arr, 0, &combs);
   splitAllBinaryStrings( &combs, nSNP, '+' );

   for( int i=0; i<fullcombs; i++) 
   {
      if(i<names.size()) {
         namesf.push_back(names[i]) ;
         combs.erase(std::remove(combs.begin(), combs.end(), names[i]), combs.end());
      } else {
         namesf.push_back(combs[i-names.size()]);
      }
   }

   Eigen::Map<Eigen::MatrixXd> ans(props.data(), props.rows(), props.cols());
   
   NumericVector ansf = flatmatcm( wrap(ans));
   ansf.attr("names") = namesf;
   for(NumericVector::iterator i = ansf.begin(); i != ansf.end(); ++i)
      *i = round( *i * 1000000000000000.0 ) / 1000000000000000.0;
   
   return(ansf);


}

