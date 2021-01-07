#include "include/CUtils.h"

using namespace Rcpp;

// Concatenate two string vectors position by position with par string in the middle.
CharacterVector concatenate( Rcpp::StringVector x, Rcpp::StringVector y, std::string par)
{
   CharacterVector res(x.size());
   
   try {
      
      for (int i=0; i < x.size(); i++)
         res[i] = ((Rcpp::as<std::string>(x[i]) + par) + Rcpp::as<std::string>(y[i]));
      
      return (res);
      
   } catch(std::exception &ex) {	
      forward_exception_to_r(ex);
   } catch(...) { 
      ::Rf_error("c++ exception (unknown reason)"); 
   }
   
}

// Duplicate array content by block of x vector
Rcpp::StringVector getNamesbyBlocks( Rcpp::StringVector x )
{
   int vlen = x.length();
   int vlength = pow(2, x.length());
   Rcpp::StringVector vnames;
   
   for( int i = 0; i<vlength; i++)
      vnames.push_back( x[i%vlen] );
   
   return(vnames);
}

// Duplicate n times each position of array x 
Rcpp::StringVector getNamesbyWord( Rcpp::StringVector x )
{
   int vlen = x.length();
   int vlength = pow(2, x.length());
   Rcpp::StringVector vnames;
   
   for( int i = 0; i<vlength; i++)
      vnames.push_back( x[(int)(i/vlen)] );

   return(vnames);
}



// Create map with elements and freqs from StringVector from start and start+length 
std::map<std::string, double> getFreq(Rcpp::StringVector Block, int start, int length)
{
   std::map<std::string, double> mapfreq; 
   int nblocks = Block.length();
   
   for (int i = 0; i < nblocks; i++) 
      mapfreq[ as<std::string>(Block(i)).substr(start,length) ]++; 

   return(mapfreq);
}




// Create map with elements and freqs from StringVector from start and start+length 
std::map<std::string, double> getFreqSum(Rcpp::NumericVector R, Rcpp::StringVector Block, int start, int length)
{
   std::map<std::string, double> mapfreq; 
   int nblocks = Block.length();
   
   for (int i = 0; i < R.length(); i++) 
      mapfreq[ as<std::string>(Block(i)).substr(start,length) ] += R[i]; 
   
   return(mapfreq);
}




// Create map with elements and freqs from StringVector from start and start+length 
std::map<std::string, double> getFreqAllCombs(Rcpp::NumericVector R, Rcpp::StringVector Block, int start, int ilength)
{
   std::map<std::string, double> mapfreq; 
   int nblocks = Block.length();
   
   int arr[ilength];
   Rcpp::StringVector combs;
   getAllBinaryStrings(ilength, arr, 0, &combs);
   
      
   for (int i = 0; i < pow(2, ilength); i++) 
      mapfreq[as<std::string>(combs(i))] = 0; 
   
   for (int i = 0; i < R.length(); i++) 
      mapfreq[ as<std::string>(Block(i)).substr(start,ilength) ] += R[i];
      
      // mapfreq[ as<std::string>(Block(i)).substr(start,ilength) ]++; 
   
   return(mapfreq);
}



// Generate all binary strings of lenght n - recursive
void getAllBinaryStrings(int n, int arr[], int i, Rcpp::StringVector *k)
{
   if (i == n) {
      std::string word;
      for (int i = 0; i < n; i++) 
         word += std::to_string(arr[i]);
      k->push_back(word);
      return;
   }
   
   // First assign "0"
   arr[i] = 0;
   getAllBinaryStrings(n, arr, i + 1, k);
   
   // Then "1"
   arr[i] = 1;
   getAllBinaryStrings(n, arr, i + 1, k);
}


// Split allele with a symbol 
void splitAllBinaryStrings( Rcpp::StringVector *combs, int nSNP, char symbol )
{
   int clenght = combs->length();
   
   for(int i=0; i<clenght; i++)
   {
      std::string comb = (std::string)combs->at(i);
      comb =  comb.substr(0,nSNP) + symbol + comb.substr(nSNP, nSNP);
      combs->at(i) = comb;
   }
      
}

// Check if all elements are equal to zero
bool All_Elements_Zero( Rcpp::RObject obj)
{
   bool allzero = false;

   if( Rf_isMatrix(obj) ){
      Eigen::MatrixXd tmpObj = as<Eigen::MatrixXd>(obj);
      Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(tmpObj.rows(), tmpObj.cols());
      if (tmpObj.isApprox(tmp)){
         Rcpp::Rcout<<"same thing"<<std::endl;
         allzero = true;
      }else
         Rcpp::Rcout << "NOT same thing"<<std::endl;

   }
   else if(Rf_isVector(obj)){
      Eigen::VectorXd tmpObj = as<Eigen::VectorXd>(obj);
      Eigen::VectorXd tmp = Eigen::VectorXd::Zero(tmpObj.size());
      if (tmpObj.isApprox(tmp)){
         Rcpp::Rcout<<"same thing"<<std::endl;
         allzero = true;
      }else
         Rcpp::Rcout << "NOT same thing"<<std::endl;
   }

   return(allzero);
}