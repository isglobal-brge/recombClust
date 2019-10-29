#ifndef CLDmixtureModel
#define CLDmixtureModel

   #include <RcppEigen.h>
   #include "CUtils.h"
   #include "CrecombFreq.h"
   #include "CLinkageFreq.h"
   #include "CMapUtils.h"
   #include "CUpdateModel.h"
   #include "CMapUtils.h"


   using namespace Rcpp;
   
   
   // define a custom template binary functor
   template<typename ch> struct reducehaplotype {
      // EIGEN_EMPTY_STRUCT_CTOR(reducehaplotype);
      String operator()(const ch& h1, const ch& h2) const { return String( h1 + h2 ); }
   };
   
   
   struct Functor {
      std::string
      operator()(const std::string& lhs, const internal::string_proxy<STRSXP>& rhs) const
      {
         return lhs + rhs;
      }
   };

   Rcpp::List cLDmixtureModel( Rcpp::RObject dat, int maxSteps, double prob0, int blocksize );
      
#endif