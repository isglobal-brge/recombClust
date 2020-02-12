#include "include/CMapUtils.h"
using namespace Rcpp;


// Convert NumericVector with names to map (key:names - value: v.value)
std::map<std::string, double> VectortoOrderedMap(NumericVector v)
{
   try 
   {
      StringVector namesv = v.names();
      std::map<std::string, double> mapv;
      
      for (size_t i = 0; i < namesv.size(); ++i)  
         mapv[(std::string)namesv[i]] = v[i];
      
      return mapv;
      
   } catch(std::exception &ex) {	
      forward_exception_to_r(ex);
   } catch(...) { 
      ::Rf_error("c++ exception (unknown reason)"); 
   }
   
}

// Convert 2dMatrix like 
//  (1,1,1,1,1,2,2,2,2,2)(1,2,3,4,5,1,2,3,4)
// to map of vectors (key:int - value: vector of values)
std::map<double, std::vector<double>> d2MatrixtoOrderedMap(Eigen::MatrixXd m)
{
   try 
   {
      std::map<double, std::vector<double>> mapv;
      
      for( size_t i=0; i< m.rows(); ++i )  {
         mapv[(double)m(i,0)].push_back(m(i,1));
      }
      
      //mapv[(double)actelem] = vd;
      
      return mapv;
      
   } catch(std::exception &ex) {	
      forward_exception_to_r(ex);
   } catch(...) { 
      ::Rf_error("c++ exception (unknown reason)"); 
   }
   
}


// Convert 2 Vectors with same length to map of vectors.
//  key : (1,1,1,1,1,2,2,2,2,2)  values: (1,2,3,4,5,1,2,3,4,5)
// to map of vectors (key:int - value: vector of values)
std::map<double, std::vector<int>> VectorstoOrderedMap(Eigen::VectorXd key, Eigen::VectorXd val)
{
   try 
   {
      std::map<double, std::vector<int>> mapv;
      
      for( size_t i=0; i< key.size(); ++i )  {
         mapv[(double)key(i)].push_back(val(i));
      }
      
      return mapv;
      
   } catch(std::exception &ex) {	
      forward_exception_to_r(ex);
   } catch(...) { 
      ::Rf_error("c++ exception (unknown reason)"); 
   }
   
}


// Assign numerical values to vector vased on map
NumericVector getNumericVectorfromStringVector(std::map<std::string, double> mapv, StringVector strvalues )
{
   try 
   {
      int vsize = strvalues.length();
      NumericVector v(vsize);

      for(int i=0; i < vsize; i++) {
         double valor = mapv[(std::string)strvalues[i]];
         v[i] = valor;
      }
         
      return v;
   
   } catch(std::exception &ex) {	
      forward_exception_to_r(ex);
   } catch(...) { 
      ::Rf_error("c++ exception (unknown reason)"); 
   }
   
}


// Sort mapped data by freq descending
std::unordered_map<std::string, double> sortMapbyFreqs( std::map<std::string, double> sNPFreq )
{
   
   // Type of Predicate that accepts 2 pairs and return a bool
   typedef std::function<bool(std::pair<std::string, double>, std::pair<std::string, double>)> Comparator;
   
   // Lambda function to compare two pairs
   Comparator compFunctor =
      [](std::pair<std::string, double> elem1 ,std::pair<std::string, double> elem2)
      {
         return elem1.second <= elem2.second;
      };
      
      // Declaring a set that will store the pairs using above comparision logic
      std::set<std::pair<std::string, double>, Comparator> setOfWords(
            sNPFreq.begin(), sNPFreq.end(), compFunctor);
      
      std::unordered_map<std::string, double> sortedsNPFreq;
      for (std::pair<std::string, double> element : setOfWords)
         sortedsNPFreq[ element.first ] = element.second;
      
      return sortedsNPFreq;
}


// flip < key, value > and sort map by value --> < value(sorted), key >
std::multimap<double, std::string, std::greater <double>> sortMultimapbyValue( std::map<std::string, double> sNPFreq )
{

   std::multimap< double, std::string, std::greater <double> > mmsortFreq;
   for (std::pair<std::string, double> element : sNPFreq)
      mmsortFreq.insert(std::pair<double, std::string>(element.second, element.first));
   
   return mmsortFreq;
}