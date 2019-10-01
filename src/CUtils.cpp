#include "include/CUtils.h"

using namespace Rcpp;

// Concatenate two string vectors position by position with par string in the middle.
CharacterVector concatenate( Rcpp::CharacterVector x, Rcpp::CharacterVector y, std::string par)
{
   CharacterVector res(x.size());
   for (int i=0; i < x.size(); i++)
   {
      res[i] = ((Rcpp::as<std::string>(x[i]) + par) + Rcpp::as<std::string>(y[i]));
   }
   return (res);
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


/***
// Sort mapped data by freq descending
std::unordered_map<std::string, int> sortMapbyFreqs( std::map<std::string, double> sNPFreq )
{
   
   // Type of Predicate that accepts 2 pairs and return a bool
   typedef std::function<bool(std::pair<std::string, int>, std::pair<std::string, int>)> Comparator;
   
   // Lambda function to compare two pairs
   Comparator compFunctor =
      [](std::pair<std::string, int> elem1 ,std::pair<std::string, int> elem2)
      {
         return elem1.second <= elem2.second;
      };
      
   // Declaring a set that will store the pairs using above comparision logic
   std::set<std::pair<std::string, int>, Comparator> setOfWords(
            sNPFreq.begin(), sNPFreq.end(), compFunctor);
   
   std::unordered_map<std::string, int> sortedsNPFreq;
   for (std::pair<std::string, int> element : setOfWords)
      sortedsNPFreq[ element.first ] = element.second;

   return sortedsNPFreq;
}
***/


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
