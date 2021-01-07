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
   
   // g <- igraph::graph_from_adjacency_matrix(OCM1, mode="undirected", weighted=TRUE, diag=FALSE, add.colnames=NA)
   // cores <- igraph::coreness(g)
   // highcore <- sum(cores>=hrstParam)
   // local_cores <- table(cores[cores>=(hrstParam)])
   // local_cores <- local_cores[local_cores>=(hrstParam)]
   
   
   NumericVector cores, local_cores;
   int highcore;
   
   
      try {
      
         Rcpp::Environment i_graph("package:igraph");
         Rcpp::Function get_graph_from_adjacency_matrix = i_graph["graph_from_adjacency_matrix"];
         Function get_coreness = i_graph["coreness"];
         
         SEXP g = get_graph_from_adjacency_matrix(as<NumericMatrix>(OCM), Named("mode", "undirected"),  
                                                   Named("weighted", true), Named("diag", false), 
                                                   Named("add.colnames",R_NaN)  );

         cores = get_coreness(Named("graph", g));
         local_cores = get_local_cores(cores, hrstParam);


   } catch(std::exception &ex) {	
      forward_exception_to_r(ex);
   } catch(...) { 
      ::Rf_error("c++ exception (unknown reason)"); 
   }
   
   return  List::create(Named("cores") = cores,  
                        Named("local_cores") = local_cores);
   
}




// Get a resume of local cores in a vector
Rcpp::RObject get_local_cores(NumericVector cores, int hrst)
{
   
   NumericVector local_c;
   
   
   try{
      
      CharacterVector names_local;   
      std::map<int, int> counts;
      for (NumericVector::iterator i = cores.begin() ; i != cores.end(); ++i) {
         if( *i >= hrst )
            ++counts[ *i ];
      }
      
      // Convert map to vector
      for( std::map<int, int>::iterator it = counts.begin(); it != counts.end(); ++it ) {
         local_c.push_back( it->second );
         names_local.push_back( static_cast<char>(it->first) );
      }
      
      local_c.names() = names_local;
      
   } catch(std::exception &ex) {
      forward_exception_to_r(ex);
   } catch(...) {
      ::Rf_error("c++ exception (unknown reason)");
   }
   
   return(local_c);
   
}




// 
// // Get graph from adjacent matrix
// Rcpp::RObject get_graph_adjacent_matrix(OCM)
// {
//    Eigen::MatrixXd inputmat = Rcpp::as<Eigen::MatrixXd>(OCM);
//    
//    try{
//       
//       
//       
//    } catch(std::exception &ex) {	
//       forward_exception_to_r(ex);
//    } catch(...) { 
//       ::Rf_error("c++ exception (unknown reason)"); 
//    }
//    
//    
//    
//    
// }


/*** R

OCM <- R
hrstParam <- 1 # valor provisional
res <- get_graph_matrix_data(OCM, hrstParam)

***/