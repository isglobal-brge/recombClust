#include "include/CgetRecombProb.h"

using namespace Rcpp;



//' Get the probability of chromosome recombination
//'
//' @param mat Matrix with model probabilities - recombClust output.
//' @param annot Matrix with initial and final coordinades of the model.
//' @param range range where the function is executed.
//' @param window window size.
//' @return A vector with chromosome probabilities to be in recomb population
//' @export
// [[Rcpp::export]]
List cGetRecombProb( NumericMatrix probmat, NumericMatrix annot, int window ) 
{
   try
   {
      Eigen::Map<Eigen::MatrixXd> eigmat (Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(probmat));

      // Obtaining namespace of GenomicRanges and IRanges packages
      Environment pkgGenomicRanges = Environment::namespace_env("GenomicRanges");
      Environment pkgIRanges = Environment::namespace_env("IRanges");
      // Picking functions from GenomicRanges and IRanges
      Function makeRangesdf = pkgGenomicRanges["makeGRangesFromDataFrame"];
      Function fOverlaps = pkgIRanges["findOverlaps"]; ;
      
      std::string chrom = "chr1";
      chrom.push_back('\0');
      
      // convert annot to dataframe and add column with chromosome 'chr1'
      DataFrame dfannot = DataFrame::create( Named("chr") = chrom,
                                             Named("start") = annot( _, 0),
                                             Named("end") = annot( _, 1)); 
      
      double intpart;
      double min = vecmin(annot( _, 0));
      double max = vecmax(annot( _, 1));
      modf( (max - min)/window, &intpart);
      double maxr = min + ( intpart * window); // gets only the range between min and max for complete windows

      NumericVector starts = generate_seq( min, maxr, window);
      
      DataFrame dfchunks = DataFrame::create( Named("chr") = chrom,
                                              Named("start") = starts,
                                              Named("end") = (starts+(window-1))); 
      
      RObject chunks = makeRangesdf(dfchunks);
      RObject grannot = makeRangesdf(dfannot);
      RObject overLaps = fOverlaps(chunks, grannot, Rcpp::_["type"] = "within");

      // get map from overlaps   
      std::map<double, std::vector<int>> mapoverlap = VectorstoOrderedMap(Rcpp::as<Eigen::VectorXd>(overLaps.slot("from")),
                                                                          Rcpp::as<Eigen::VectorXd>(overLaps.slot("to")));

      Eigen::MatrixXd res(eigmat.rows(), mapoverlap.size());
      
      size_t j=0;
      for (std::pair<double, std::vector<int>> element : mapoverlap)
      {
         // Accessing VALUE from element.
         std::vector<int> values = element.second;
         Eigen::MatrixXd tmp( eigmat.rows(), values.size() );

         for(size_t i=0; i< values.size(); i++){
            tmp.col(i) = eigmat.col( (values[i]-1) ); 
         }
         
         res.col(j) = tmp.rowwise().mean();
         j++;
      }
       
      return  List::create(Named("res") = wrap(res),
                           Named("annot") = wrap(grannot),
                           Named("overlaps") = overLaps,
                           Named("mapoverlaps") = wrap(mapoverlap),
                           Named("dfchunks") = wrap(dfchunks)
                           );
      
   } catch(std::exception &ex) {	
      forward_exception_to_r(ex);
   } catch(...) { 
      ::Rf_error("c++ exception (unknown reason)"); 
   }
   
}

/*** R



## DEBUG
#..# ######################
#..# R --debugger=lldb
#..# run
#..# ######################



*/