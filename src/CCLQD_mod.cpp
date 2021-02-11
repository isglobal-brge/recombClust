#include "include/CCLQD_mod.h"

using namespace Rcpp;




// I have to MODIFY THIS INFORMATION WITH RECOMBCLUST SPECIFICS !!!!!!

//' @title Estimation of LD block regions
//' @description \code{Big_LD} returns the estimation of LD block regions of given data.
//' @param OCM Correlation matrix
//' @param CLQcut A numeric value of threshold for the correlation value |r|, between 0 to 1.
//' @param clstgap  An integer value to specifying the threshold of physical distance (bp) between two consecutive SNPs
//' which do not belong to the same clique, i.e., if a physical distance between two consecutive SNPs in a clique
//' greater than \code{clstgap}, then the algorithm split the cliques satisfying each
//' clique do not contain such consecutive SNPs
//' @param hrstType, string value with hrstType=c("near-nonhrst", "fast", "nonhrst")
//' @param leng  An integer value to specify the number of SNPs in a preceding and a following region
//' of each sub-region boundary, every SNP in a preceding and every SNP in a following region need to be in weak LD.
//' @param MAFcut An numeric value to specifying the MAF threshold. 
//' @param subSegmSize  An integer value to specify the upper bound of the number of SNPs in a one-take sub-region.
//' @param appendRare If \code{appendRare = TRUE}, the algorithm append rare SNPs (MAF<MAFcut) to the constructed LD blocks or add a new LD blocks
//' @param checkLargest If \code{checkLargest = TRUE}, the algorithm use heuristic procedure to reduce runtime of CLQ-D execution
//' @param CLQmode A character string to specify the way to give priority among detected cliques, the algorithm gives priority to the largest clique
//' 
//' <output>
//' @return  A data frame of block estimation result.
//' Each row of data frame shows the starting SNP and end SNP of each estimated LD block.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject CLQD_mod( Rcpp::RObject CorMat,
                        double CLQcut = 0.2, int clstgap = 40000, 
                        std::string hrstType = "fast", int hrstParam = 200,
                        std::string CLQmode = "maximal")
{
   
   int iskipRatio = 0;
   std::vector<std::string> vhrstType {"near-nonhrst", "fast", "nonhrst"};
   std::vector<std::string> vCLQmode {"mmaximal"};
   Rcpp::NumericMatrix OCM = as<NumericMatrix>(CorMat);
   bool heuristic = false;
   
   // std::vector<std::string> OCMColnames = colnames(OCM);
   
   int binnum = 1;
   
   
   if (!in_array(hrstType, vhrstType)) {
      throw std::range_error("Inadmissible value for variable hrstType, possible values are : 'near-nonhrst', 'fast' or 'nonhrst'"); 
   }
   if (!in_array(hrstType, vhrstType)) {
      throw std::range_error("Inadmissible value for variable hrstType");
   }
   
   // Group windows higher than CLQcut
   OCM = arrange_matrix_by_cutoff(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(OCM), CLQcut);
   
   // Create numericvector with na values  (New rcppversion : NumericVector m(n, NA) ;)
   NumericVector binvector(OCM.cols());
   std::fill( binvector.begin(), binvector.end(), NumericVector::get_na() ) ;

   Rcpp::NumericMatrix OCM1 = OCM;
   // Aquí seguirem endavant o en un altre cas, retornarem binvector
   if( !All_Elements_Zero(OCM) ){
      
      Rcpp::Rcout<<"\nSeguim endavant\n";
      
      if( hrstType.compare("near-nonhrst") == 0)
      {
         heuristic = true;
         int heuristicNum = 0;
         while(heuristic == true)
         {
            List coresdata = as<List>(get_graph_matrix_data(OCM1, hrstParam));
            NumericVector cores = coresdata["cores"];
            Rcpp::NumericVector local_cores = coresdata["local_cores"];
            int local_hrstParam = coresdata["local_hrstParam"];
            NumericVector bothhighSNPs = coresdata["bothhighSNPs"];
            
            Rcpp::Rcout<<"bothhighSNPs : \n"<<bothhighSNPs.length();

            if( local_cores.length()>0 ){
            
               if( heuristicNum == 0 ){
                  heuristicNum = heuristicNum + 1;
                  Rcpp::Rcout<<"\nUse near-CLQ heuristic procedure!!\t - \t heuristic loop"<<heuristicNum<<"\n";
               }else{
                  heuristicNum = heuristicNum+1;
                  Rcpp::Rcout<<"heuristic loop"<<heuristicNum<<"\n";
               }
               
////// WORKING IN PROGRESS !!!!               
/***              
               NumericVector SNPset1 = testNAIndex(binvector, bothhighSNPs);
               

               // NumericVector SNPset1 = 
               
               //                      SNPset1 <- which(is.na(binvector))[bothhighSNPs]
               //                      
               //                      // Fins aquí, el que fan aquestes línies no es més que mirar els nodes on els local_cores és màxim.
               //                      només caldria trobar les posicions .... revisar-ho en general perquè potser no caldria fer-ho així
               //                      i seria fins i tot més senzill..... 
               //                      

                              

               Rcpp::NumericMatrix nowOCM = removeMatrixColsandRows(OCM,SNPset1);
               Rcpp::NumericMatrix nowOCM = removeMatrixColumns(OCM,SNPset1);
               nowOCM = removeMatrixRows(nowOCM,SNPset1);
 ***/
////// WORKING IN PROGRESS !!!!

               
               //                   nowOCM <- OCM[SNPset1, SNPset1]
               //                   heuristicBins <- heuristicCLQ(nowOCM, hrstParam)
               //                      binvector[SNPset1[heuristicBins]]<-binnum
               //                      binnum <- binnum+1
               //                   OCM1 <- OCM[is.na(binvector), is.na(binvector)]
               
               
               heuristic = false;
            } else {
               heuristic = false;
            }
            

            
//             while(heuristic == TRUE){
//                g <- igraph::graph_from_adjacency_matrix(OCM1, mode="undirected", weighted=TRUE, diag=FALSE, add.colnames=NA)
//                cores <- igraph::coreness(g)
//                highcore <- sum(cores>=hrstParam)
//                local_cores <- table(cores[cores>=(hrstParam)])
//                local_cores <- local_cores[local_cores>=(hrstParam)]
//                if(length(local_cores>0)){
//                   if(heuristicNum==0){
//                      message("Use near-CLQ heuristic procedure!!")
//                      heuristicNum <- heuristicNum+1
//                      message(paste("hueristic loop", heuristicNum))
//                   }else{
//                      heuristicNum <- heuristicNum+1
//                      message(paste("hueristic loop", heuristicNum))
//                   }
// # find dense region
// # local_cores <- table(cores[cores>=(hrstParam)])
// # local_cores <- local_cores[local_cores>=(hrstParam)]
//                   local_hrstParam <- names(local_cores[which(as.integer(names(local_cores))==max(as.integer(names(local_cores))))])
//                      local_hrstParam <- as.numeric(local_hrstParam)
//                      bothhighSNPs <- which(cores == local_hrstParam)
//                      SNPset1 <- which(is.na(binvector))[bothhighSNPs]
//                   nowOCM <- OCM[SNPset1, SNPset1]
//                   heuristicBins <- heuristicCLQ(nowOCM, hrstParam)
//                      binvector[SNPset1[heuristicBins]]<-binnum
//                      binnum <- binnum+1
//                   OCM1 <- OCM[is.na(binvector), is.na(binvector)]
//                }else{
//                   heuristic <- FALSE
//                }
//             } #end while
            
         }
      }
    
   }
      

   
   
   return(binvector);
}


// Modify correlation matrix attending to : 
//    Diagonal matrix = 0
//    NA values = 0
//    values lower than cutoff = 0
//    values higher or equal to cutoff = 1
NumericMatrix arrange_matrix_by_cutoff(Eigen::MatrixXd CorMat, double cutoff)
{
   
   CorMat.diagonal() = Eigen::VectorXd::Zero(CorMat.diagonalSize());

   Eigen::MatrixXd mask0 = Eigen::MatrixXd::Zero(CorMat.rows(), CorMat.cols());
   Eigen::MatrixXd mask1 = Eigen::MatrixXd::Constant(CorMat.rows(), CorMat.cols(), 1);

   CorMat = (CorMat.array().real() < cutoff).select(mask0 ,mask1 ); 

   return(wrap(CorMat));
   
}



/*** R

library(recombClust)
a <- matrix(c(1,2,3,4,5,1,2,0,3,0,1,0,0,2,5,1,2,4,4,5,0,0,0,4,5), byrow = TRUE, nrow = 5)
CLQD_mod(a, hrstType = "near-nonhrst",  hrstParam = 1)


CLQD_mod <- function(OCM, CLQcut=0.2, clstgap=40000, hrstType=c("near-nonhrst", "fast", "nonhrst"), hrstParam=200,
                     CLQmode=c("density", "maximal"))
{
   #..# #'#######################################################################################################
   #..# skipRatio <- 0
   #..# CLQmode <- match.arg(CLQmode)
   #..# hrstType <- match.arg(hrstType)
   
   #..# ## Group windows higher than CLQcut
   #..# diag(OCM) <- 0
   #..# OCM[is.na(OCM)] <- 0
   #..# OCM[abs(OCM) < CLQcut] <- 0
   #..# OCM[abs(OCM) >= CLQcut] <- 1
   
   #..# # SNPinfo<-SNPinfo[!allNASNPs, ]
   
   #..# # Main Function
   #..# SNPbps <- as.numeric(gsub("-.*", "", colnames(OCM)))
   #..# binvector <- rep(NA, dim(OCM)[2])
   #..# binnum <- 1
   #..# # re.SNPbps <- SNPbps
   #..# if(all(OCM==0)) return(binvector)
   
   ## Inicialització per al testeig : 
         OCM <- a
         # hrstParam <- 1
         # skipRatio <- 0
         # 
         # ## Group windows higher than CLQcut
         # diag(OCM) <- 0
         # OCM[is.na(OCM)] <- 0
         # OCM[abs(OCM) < CLQcut] <- 0
         # OCM[abs(OCM) >= CLQcut] <- 1
         # 
         # # SNPinfo<-SNPinfo[!allNASNPs, ]
         # 
         # # Main Function
         # SNPbps <- as.numeric(gsub("-.*", "", colnames(OCM)))
         # binvector <- rep(NA, dim(OCM)[2])
         # binnum <- 1
         # re.SNPbps <- SNPbps
         if(all(OCM==0)) return(binvector)
         
   ## Fi Inicialització per al testeig : ##
   
   # test graph complexity
   OCM1 <- OCM
   # bothhigh <- (degs>=400 & cores>=400)
   if(hrstType == "nonhrst"){
      heuristic <- FALSE
   }else if(hrstType == "near-nonhrst"){
      heuristic <- TRUE
      heuristicNum <- 0
      while(heuristic == TRUE){
         g <- igraph::graph_from_adjacency_matrix(OCM1, mode="undirected", weighted=TRUE, diag=FALSE, add.colnames=NA)
         cores <- igraph::coreness(g)
         highcore <- sum(cores>=hrstParam)
         local_cores <- table(cores[cores>=(hrstParam)])
         local_cores <- local_cores[local_cores>=(hrstParam)]
         if(length(local_cores)>0){     ######### REVISAR AQUESTA CONDICIÓ A L'ALGORITME ORIGINAL !!! AQUÍ ESTAVA MALAMENT I JO L'HE MODIFICAT !!!
            if(heuristicNum==0){
               message("Use near-CLQ heuristic procedure!!")
               heuristicNum <- heuristicNum+1
               message(paste("hueristic loop", heuristicNum))
            }else{
               heuristicNum <- heuristicNum+1
               message(paste("hueristic loop", heuristicNum))
            }
            # find dense region
            # local_cores <- table(cores[cores>=(hrstParam)])
            # local_cores <- local_cores[local_cores>=(hrstParam)]
            local_hrstParam <- names(local_cores[which(as.integer(names(local_cores))==max(as.integer(names(local_cores))))])
            local_hrstParam <- as.numeric(local_hrstParam)
            bothhighSNPs <- which(cores == local_hrstParam)
            SNPset1 <- which(is.na(binvector))[bothhighSNPs]
            nowOCM <- OCM[SNPset1, SNPset1]
            heuristicBins <- heuristicCLQ(nowOCM, hrstParam)
            binvector[SNPset1[heuristicBins]]<-binnum
            binnum <- binnum+1
            OCM1 <- OCM[is.na(binvector), is.na(binvector)]
         }else{
            heuristic <- FALSE
         }
      } #end while
   }else if(hrstType == "fast"){
      checkLargest <- TRUE
      r2Mat <- OCM
      re.SNPbps <- as.numeric(gsub("-.*", "", colnames(OCM)))
      # firstTerm = T
      while(checkLargest == TRUE){
         g <- igraph::graph_from_adjacency_matrix(r2Mat, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NA)
         compo = igraph::components(g)
         componum = which(compo$csize==max(compo$csize))[1]
         compov = which(compo$membership==componum)
         compadjM = OCM1[compov, compov]
         cg = igraph::graph_from_adjacency_matrix(compadjM, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NA)
         if((median(igraph::coreness(cg))>80 & max(igraph::coreness(cg))>100)| (quantile(igraph::coreness(cg), 0.75)>100 & max(igraph::coreness(cg))>100)){
            message("use fast heuristic procedure!")
            # if(quantile(degrees, 0.7) == 1) break
            # message(head((quantile(degrees, 0.7))))
            degrees = apply(r2Mat, 1, sum)
            maxdegv = which(degrees >=max(quantile(degrees, 0.7), 80))
            # if(length(maxdegv)>=1){
            maxdegvs = maxdegv
            edgeDens = NULL
            for(maxdegv in maxdegvs){
               Bignbds = which(r2Mat[maxdegv,, drop = FALSE]>0, arr.ind = TRUE)
               Bignbds.c = unique(Bignbds[,2])
               newr2Mat = r2Mat[Bignbds.c,Bignbds.c]
               EdgeDen = sum(newr2Mat)/((dim(newr2Mat)[1])*(dim(newr2Mat)[1]-1))
               edgeDens = c(edgeDens, EdgeDen)
            }
            maxdegvs = maxdegvs[order(edgeDens, decreasing = TRUE)]
            edgeDens = edgeDens[order(edgeDens, decreasing = TRUE)]
            degv = maxdegvs[1]
            edgeD = edgeDens[1]
            Bignbds = which(r2Mat[degv,, drop = FALSE]>0, arr.ind = TRUE)
            Bignbds.c = unique(Bignbds[,2])
            # maxiC = maximal.cliques(g, min = dim(OCM1)[1]*0.9)
            # largestOneRange = range(Bignbds.c)
            # largestSNPn = diff(largestOneRange)
            # largestCsize = length(Bignbds.c)
            nowSNPsbp = re.SNPbps[Bignbds.c]
            nowSNPsbploca = match(nowSNPsbp, SNPbps)
            binvector[nowSNPsbploca] <- binnum
            binnum = binnum + 1
            r2Mat <- r2Mat[-Bignbds.c, -Bignbds.c, drop = FALSE]
            OCM1 <- OCM1[-Bignbds.c, -Bignbds.c, drop = FALSE]
            re.SNPbps <- re.SNPbps[-Bignbds.c]
            # message("case2")
            checkLargest = TRUE
            if(length(re.SNPbps)<500)  checkLargest = FALSE
         }else{
            checkLargest = FALSE
         }
         
      }
   }
   # binvector splitting by clstgap
   if(all(is.na(binvector))==FALSE){
      binveclist <- lapply(seq_len(max(binvector, na.rm = TRUE)),
                           function(x) SNPbps[which(binvector == x)])
      binvecbplist <- newSplitCliques(binveclist, clstgap)
      binvecbpindex <- lapply(binvecbplist, function(x) match(x, SNPbps))
      binvecbpindex <- binvecbpindex[vapply(binvecbpindex, length, c(1))>1]
      binvector <- rep(NA, length(binvector))
      for(i in seq_len(length(binvecbpindex))){
         binvector[binvecbpindex[[i]]]<-i
      }
   }
   # take Toooo Big block First!
   r2Mat <- OCM[is.na(binvector), is.na(binvector)]
   if(sum(is.na(binvector))<=1 | (sum(r2Mat)==0)) return(binvector)
   g <- igraph::graph_from_adjacency_matrix(r2Mat, mode="undirected", weighted=TRUE, diag=FALSE, add.colnames=NA)
   max.cliques <- igraph::max_cliques(g, min=2)
   if(length(max.cliques)==0) stop("max.cliques is empty")
   re.SNPbps <- SNPbps[is.na(binvector)]
   bp.cliques <- lapply(max.cliques, function(x) re.SNPbps[x])
   split.bp.cliques <- newSplitCliques(bp.cliques, clstgap)
   # reduce the number of maximal clique? (candidate) or
   # modify density function. narrow SNP distance, so small cliques are chosen preferencely.
   repeat {
      # message(binnum)
      if (all(is.na(binvector) == FALSE)) {
         break
      }
      if(length(split.bp.cliques)==0) break
      now.split.bp.cliques <- split.bp.cliques
      if(CLQmode =="density"){
         density.v <- vapply(now.split.bp.cliques, function(x) ((length(x)))/(diff(range(x))/1000), 1)
      }else{
         density.v <- vapply(now.split.bp.cliques, length, 1)
      }
      
      max.d <- which(density.v == max(density.v))
      max.cluster <- now.split.bp.cliques[max.d]
      if (length(max.cluster) > 1) {
         # if there are two bins of same density, then we choose the bigger one.
         max.cluster <- max.cluster[order(vapply(max.cluster, length, 1), decreasing=TRUE)]
      }
      max.cluster <- max.cluster[[1]]
      max.cluster.od <- match(max.cluster, re.SNPbps)
      # if (codechange == TRUE) {
      #   max.cluster.od <- ChooseMaximal(max.cluster.od, CLQcut, OCM)
      #   max.cluster <- re.SNPbps[max.cluster.od]
      # }
      ## excluding all SNPs in max.cluster from re.SNPbps
      split.bp.cliques <- lapply(split.bp.cliques, function(x) setdiff(x, max.cluster))
      split.bp.cliques <- unique(split.bp.cliques)
      split.bp.cliques <- split.bp.cliques[which(vapply(split.bp.cliques, length, 1) > 1)]
      binvector[match(max.cluster, SNPbps)] <- binnum
      binnum=binnum + 1
      # r2Mat <- r2Mat[-max.cluster.od, -max.cluster.od]
      # OCM <- OCM[-max.cluster.od, -max.cluster.od]
      # re.SNPbps <- setdiff(re.SNPbps, max.cluster)
      if (length(re.SNPbps) < 2) {
         break
      }
      if(length(split.bp.cliques)==0) break
      # message(sum(is.na(binvector)))
   }  ##end repeat
   return(binvector)
}

*/
