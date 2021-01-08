#include "include/CGraphUtils.h"

using namespace Rcpp;

//' @title Get graphs data from addjacency matrix
//' @description \code{Big_LD} returns adjacency matrix graph data to test graph complexity
//' @param OCM Correlation matrix
//' @param hrstParam A numeric value of threshold for the correlation value |r|, between 0 to 1.
//' 
//' <output>
//' @return  A list with , cores, highcore and local_cores
//' \itemize{
//'  \item{"adjMatrix"}{adjacent matrix graph}
//'  \item{"cores"}{a maximal subgraph in which each vertex has at least degree k. }
//'  \item{"highcore"}{Number of cores with "cores">hrstParam}
//'  \item{"local_cores"}{local_cores}
//' }
//' 
//' @importFrom igraph graph_from_adjacency_matrix coreness
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject get_graph_matrix_data( Rcpp::RObject OCM, int hrstParam )
{
   
   NumericVector cores, bothhighSNPs; 
   std::map<int, std::vector<int>> local_cores;
   int highcore, local_hrstParam;
   
   
      try {
      
         Rcpp::Environment i_graph("package:igraph");
         Rcpp::Function get_graph_from_adjacency_matrix = i_graph["graph_from_adjacency_matrix"];
         Function get_coreness = i_graph["coreness"];
         
         SEXP g = get_graph_from_adjacency_matrix(as<NumericMatrix>(OCM), Named("mode", "undirected"),  
                                                   Named("weighted", true), Named("diag", false), 
                                                   Named("add.colnames",R_NaN)  );

         cores = get_coreness(Named("graph", g));
         local_cores = get_local_cores(cores, hrstParam);
         local_hrstParam = local_cores.rbegin()->first;
         bothhighSNPs = local_cores[static_cast<char>(local_hrstParam)];
         
         /* DEBUG ONLY :
            Rcpp::Rcout<<"\n CoresData : "<< cores;
            Rcpp::Rcout<<"\n local_hrstParam : "<< local_hrstParam;
          */
            Rcpp::Rcout<<"\n bothhighSNPs : "<< bothhighSNPs;
         


   } catch(std::exception &ex) {	
      forward_exception_to_r(ex);
   } catch(...) { 
      ::Rf_error("c++ exception (unknown reason)"); 
   }
   
   return  List::create(Named("cores") = cores,  
                        Named("local_cores") = get_localcores_map_as_vector(local_cores),
                        Named("local_hrstParam") = local_hrstParam,
                        Named("bothhighSNPs") = bothhighSNPs);
   
}




// Get a resume of local cores in a vector
//..// Rcpp::RObject get_local_cores(NumericVector cores, int hrst)
std::map<int, std::vector<int>> get_local_cores(NumericVector cores, int hrst)
{
   
   //..// NumericVector local_c;
   std::map<int, std::vector<int>> counts;
   
   try{
      
      //..// CharacterVector names_local;   
      int position = 0;
      
      for (NumericVector::iterator i = cores.begin() ; i != cores.end(); ++i) {
         if( *i >= hrst ){
            counts[ *i ].push_back(position) ;
         }
         Rcpp::Rcout<<"Posicio val : "<<position<<"\n";
         position++;
      }
      
      /* DEBUG ONLY :
         // Print data in map (sorted??)
         Rcpp::Rcout << "\nMap content : \n";
         for(auto elem : counts) {
            Rcpp::Rcout << elem.first << " :\n";
            for(auto it2 = elem.second.begin(); it2 != elem.second.end(); ++it2)
               Rcpp::Rcout <<"\t"<< *it2 ;
            Rcpp::Rcout <<"\n";
         }
      */
      
   } catch(std::exception &ex) {
      forward_exception_to_r(ex);
   } catch(...) {
      ::Rf_error("c++ exception (unknown reason)");
   }
   
   //..// return(local_c);
   return(counts);
   
}


// Convert std::map<int, std::vector<int>> map to vector with numbers of values for each key
NumericVector get_localcores_map_as_vector(std::map<int, std::vector<int>> localCores)
{
   CharacterVector names_local;
   NumericVector local_c;
   
   for( std::map<int, std::vector<int>>::iterator it = localCores.begin(); it != localCores.end(); ++it ) {
      local_c.push_back( it->second.size() );
      names_local.push_back( std::to_string(it->first) );
   }
   
   local_c.names() = names_local;
   
   return(local_c);
}


/*** R
a <- matrix(c(1,2,3,4,5,1,2,0,0,0,1,0,0,2,5,1,2,4,4,5,1,0,0,0,0), byrow = TRUE, nrow = 5)
# OCM <- R
OCM <- a
hrstParam <- 1 # valor provisional
res <- get_graph_matrix_data(OCM, hrstParam)

***/