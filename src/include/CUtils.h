#ifndef CUtils
#define CUtils

   #include <Rcpp.h>
   #include <map>
   #include <unordered_map>

   // [[Rcpp::plugins(cpp11)]]

   Rcpp::CharacterVector concatenate( Rcpp::CharacterVector x, Rcpp::CharacterVector y, std::string par);
   
   std::map<std::string, double> getFreq(Rcpp::StringVector Block, int start, int length);
   std::map<std::string, double> getFreqSum(Rcpp::NumericVector R, Rcpp::StringVector Block, int start, int length);
   std::map<std::string, double> getFreqAllCombs(Rcpp::NumericVector R, Rcpp::StringVector Block, int start, int ilength);
      
   int getCombinations (int size, Rcpp::StringVector symbol);
   void getAllBinaryStrings(int n, int arr[], int i, Rcpp::StringVector *k);
   Rcpp::StringVector getNamesbyBlocks( Rcpp::StringVector x );
   Rcpp::StringVector getNamesbyWord( Rcpp::StringVector x );
   void splitAllBinaryStrings( Rcpp::StringVector *combs, int nSNP, char symbol );

      
#endif