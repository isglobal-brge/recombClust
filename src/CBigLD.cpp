#include "include/CBigLD.h"

using namespace Rcpp;


// I have to MODIFY THIS INFORMATION WITH RECOMBCLUST SPECIFICS !!!!!!

//' @title Estimation of LD block regions
//' @description returns the estimation of LD block regions of given data.
//' @param geno A data frame or matrix of additive genotype data, each column is additive genotype of each SNP.
//' @param SNPinfo A data frame or matrix of SNPs information.  1st column is rsID and 2nd column is bp position.
//' @param CLQcut A numeric value of threshold for the correlation value |r|, between 0 to 1.
//' @param clstgap  An integer value to specifying the threshold of physical distance (bp) between two consecutive SNPs
//' which do not belong to the same clique, i.e., if a physical distance between two consecutive SNPs in a clique
//' greater than \code{clstgap}, then the algorithm split the cliques satisfying each
//' clique do not contain such consecutive SNPs
//' @param leng  An integer value to specify the number of SNPs in a preceding and a following region
//' of each sub-region boundary, every SNP in a preceding and every SNP in a following region need to be in weak LD.
//' @param MAFcut An numeric value to specifying the MAF threshold. 
//' @param subSegmSize  An integer value to specify the upper bound of the number of SNPs in a one-take sub-region.
//' @param appendRare If \code{appendRare = TRUE}, the algorithm append rare SNPs (MAF<MAFcut) to the constructed LD blocks or add a new LD blocks
//' @param checkLargest If \code{checkLargest = TRUE}, the algorithm use heuristic procedure to reduce runtime of CLQ-D execution
//' @param hrstType, string value with hrstType=c("near-nonhrst", "fast", "nonhrst")
//' 
//' <output>
//' @return  A data frame of block estimation result.
//' Each row of data frame shows the starting SNP and end SNP of each estimated LD block.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject BigLD( Rcpp::RObject corMat = R_NilValue , Nullable<double> CLQcut = 0.5, 
                     Nullable<int> clstgap = 40000, std::string CLQmode = "maximal", // Tot i que el paràmetre CLQmode no cal.
                     std::string hrstType = "fast", Nullable<int> hrstParam = 200,
                     std::string chrN = "chr1")
{
   
   std::vector<std::string> vhrstType {"near-nonhrst", "fast", "nonhrst"};
   // std::string str_hrstTypex = Rcpp::as<std::string>(hrstType);
   
   if (!in_array(hrstType, vhrstType)) {
      throw std::range_error("Inadmissible value for variable hrstType, possible values are : 'near-nonhrst', 'fast' or 'nonhrst'");
   }
   
   if (!in_array(hrstType, vhrstType)) {
      throw std::range_error("Inadmissible value for variable hrstType");
   }
   
   //       hrstType= c("near-nonhrst", "fast", "nonhrst")
   
   
   /***
   subbinvec <- CLQD_mod(OCM=corMat, CLQcut=CLQcut, clstgap=clstgap,
                         hrstType=hrstType, hrstParam = hrstParam, CLQmode = "maximal")
   cat('CLQ done!\r')
   if(all(is.na(subbinvec) == TRUE)) next;
   bins <- seq_len(max(subbinvec[which(!is.na(subbinvec))]))
      clstlist <- lapply(bins, function(x) which(subbinvec  == x))
      clstlist <- lapply(clstlist, sort)  ###
   clstlist <- clstlist[order(vapply(clstlist, min, 1))]  ###
   nowLDblocks <- constructLDblock(clstlist)
      nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop = FALSE]
   cat('constructLDblock done!\r')
      
      start.index <- nowLDblocks[,1]
   end.index <- nowLDblocks[,2]
   start.bp <- gsub("-.*", "", colnames(corMat)[nowLDblocks[,1]])
      end.bp <- gsub(".*-", "", colnames(corMat)[nowLDblocks[,2]])
      recombBlocks <- data.frame(chr = chrN, start = start.bp, end = end.bp, start.idx = start.index, end.idx = end.index)
      
      recombBlocksGR <- makeGRangesFromDataFrame(recombBlocks, keep.extra.columns = TRUE)
      message("\nBigLD done!")
      return(recombBlocksGR)
   ***/
   return(wrap(1));
 
}



/*** R




BigLD <- function(corMat=NULL,  cutByForce=NULL, CLQcut=0.5, clstgap=40000, CLQmode=c("density", "maximal"),
                  leng=200, subTaskSize=1500, MAFcut=0.05, appendRare=FALSE,
                  hrstType=c("near-nonhrst", "fast", "nonhrst"), hrstParam=200,
                  chrN="chr1", startbp=-Inf, endbp=Inf)
{
   
   ##  Possibles valors : 
   ##       hrstType= c("near-nonhrst", "fast", "nonhrst")
   
   
   skipRatio=0.0
   CLQmode <- match.arg(CLQmode)
   
   
   subbinvec <- CLQD_mod(OCM=corMat, CLQcut=CLQcut, clstgap=clstgap,
                         hrstType=hrstType, hrstParam = hrstParam, CLQmode = "maximal")
   cat('CLQ done!\r')
   if(all(is.na(subbinvec) == TRUE)) next;
   bins <- seq_len(max(subbinvec[which(!is.na(subbinvec))]))
   clstlist <- lapply(bins, function(x) which(subbinvec  == x))
   clstlist <- lapply(clstlist, sort)  ###
   clstlist <- clstlist[order(vapply(clstlist, min, 1))]  ###
   nowLDblocks <- constructLDblock(clstlist)
   nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop = FALSE]
   cat('constructLDblock done!\r')
   
   start.index <- nowLDblocks[,1]
   end.index <- nowLDblocks[,2]
   start.bp <- gsub("-.*", "", colnames(corMat)[nowLDblocks[,1]])
   end.bp <- gsub(".*-", "", colnames(corMat)[nowLDblocks[,2]])
   recombBlocks <- data.frame(chr = chrN, start = start.bp, end = end.bp, start.idx = start.index, end.idx = end.index)
   
   recombBlocksGR <- makeGRangesFromDataFrame(recombBlocks, keep.extra.columns = TRUE)
   message("\nBigLD done!")
   return(recombBlocksGR)
}
*/
