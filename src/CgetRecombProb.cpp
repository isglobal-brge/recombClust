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
      
      
      // std::vector<std::string> chr(annot(_,0).size(), "chr1");
      
      std::string chrom = "chr1";
      chrom.push_back('\0');
      
      // convert annot to dataframe and add column with chromosome 'chr1'
      DataFrame dfannot = DataFrame::create( Named("chr") = chrom,
                                             Named("start") = annot( _, 0),
                                             Named("end") = annot( _, 1)); 
      
      double min = vecmin(annot( _, 0));
      double max = vecmax(annot( _, 1));
      NumericVector starts = generate_seq( min, max, window);
      
      DataFrame dfchunks = DataFrame::create( Named("chr") = chrom,
                                              Named("start") = starts,
                                              Named("end") = (starts+(window-1))); 
      
      RObject chunks = makeRangesdf(dfchunks);
      RObject grannot = makeRangesdf(dfannot);
      RObject overLaps = fOverlaps(chunks, grannot);
      
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
               tmp.col(i) = eigmat.col( values[i] ); 
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


## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# OBTENIR LES MATRIUS RESULTATS DE RECOMBCLUST
#####
library(SeqArray)
library(recombClust)
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Github/recombClust-C")
# save.image("~/Library/Mobile Documents/com~apple~CloudDocs/UVic/TFM/Treballant/recombClust/VarSortidaCorr1.RData")
# load("~/Library/Mobile Documents/com~apple~CloudDocs/UVic/TFM/Treballant/recombClust/VarSortidaCorr1.RData")


# Fitxer vcf : "inst/extdata/example.vcf"
# Fitxer gds : filename <- "inst/extdata/example.gds"
# Per si s'han de tancar tots els fitxers .gds : 
#     showfile.gds(closeall=TRUE) 

res <- CRecombClust("inst/extdata/example.vcf", 7, "54301673", "54306217"  )

# OBTENIR LA PROBABILITAT DE RECOMBINACIÓ DE LES VARIABLES
#####

# Desllistem
    mat <- res$models[[1]]$r1
    annot <- res$models[[1]]$annot
    for (i in 2:length(res$models))
    {
       mat <- rbind(mat, res$models[[i]]$r1)
       annot <- rbind(annot, res$models[[i]]$annot)
    }

# sortida <- cGetRecombProb( as.matrix(t(res$mat)), as.matrix(annot),2)
# dim(sortida$res)


# OBTENIR LA MATRIU DE CORRELACIONS DE LA MATRIU DE PROBABILITATS DE RECOMBINACIÓ
#####

## RecombProb + Matriu de correlacions 
sortidaf <- CgetCorrRecombClust(as.matrix(t(mat)), as.matrix(annot),2)

## Write data to hdf5 file
Create_HDF_HiCBricks_File("test/Fitxer_hdf5.h5", sortidaf$corrmat, sortidaf$recombmat$dfchunk)

### FINS AQÚ ESCRIPTURA FITXER !! --> FALTA ACABAR DE FER PROVES AMB DADES BONES ABANS DE PUJAR-HO
### TOT CAP A GITHUB --> PREVIST FER-HO DEMÀ A LA TARDA.


## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# https://slack-redir.net/link?url=https%3A%2F%2Fbioconductor.org%2Fpackages%2Frelease%2Fbioc%2Fvignettes%2FHiCBricks%2Finst%2Fdoc%2FIntroductionToHiCBricks.html
# https://bioconductor.org/packages/release/bioc/html/HiCBricks.html
# https://bioconductor.org/packages/release/bioc/manuals/HiCBricks/man/HiCBricks.pdf
# 




## TRY HiCBricks

write.table(sortidaf$corrmat, "test/mattrix.txt", append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = TRUE)

write.table(sortidaf$recombmat$dfchunk, "test/bintable.txt" , append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = TRUE)


## DADES DEL FITXER bintable.txt i mattrix.txt

# Load bintable
Bintable_path <- file.path("test","bintable.txt")
out_dir <- file.path(tempdir(), "HiCBricks_test1")
dir.create(out_dir)
Create_many_Bricks(BinTable = Bintable_path, 
                   bin_delim=" ", output_directory = out_dir, 
                   file_prefix = "HiCBricks_test1", remove_existing=TRUE, 
                   experiment_name = "HiCBricks test1", resolution = 100000)


# Load mattrix
BrickContainer_dir <- file.path(tempdir(), "HiCBricks_test1")
My_BrickContainer <- load_BrickContainer(project_dir = BrickContainer_dir)
Example_dataset_dir <- system.file("extdata", package = "HiCBricks")

Chromosomes <- c("chr1")
for (chr in Chromosomes) {
   
#   Matrix_file <- file.path(Example_dataset_dir,
#                            paste(paste("mattrix.txt", 
#                                        chr, sep = "_"), "txt.gz", sep = "."))
   
   Matrix_file <- file.path("test","mattrix.txt")
   Brick_load_matrix(Brick = My_BrickContainer,
                     chr1 = chr,
                     chr2 = chr,
                     resolution = 100000,
                     matrix_file = Matrix_file,
                     delim = " ",
                     remove_prior = TRUE)
}


Brick_get_bintable(My_BrickContainer, resolution = 100000)
Brick_return_region_position(Brick = My_BrickContainer,
                             region = "chr1:54301675:54305898", resolution = 100000)


# Query mattrix data
Brick_list_matrix_mcols()

MCols.dat <- Brick_get_matrix_mcols(Brick = My_BrickContainer,
                                    chr1 = "chr1",
                                    chr2 = "chr1",
                                    resolution = 100000,
                                    what = "chr1_row_sums")
head(MCols.dat, 100)

Brick_matrix_isdone(Brick = My_BrickContainer,
                    chr1 = "chr1",
                    chr2 = "chr1",
                    resolution = 100000)

Brick_matrix_issparse(Brick = My_BrickContainer,
                      chr1 = "chr1",
                      chr2 = "chr1",
                      resolution = 100000)

Brick_matrix_maxdist(Brick = My_BrickContainer,
                     chr1 = "chr1",
                     chr2 = "chr1",
                     resolution = 100000)

# chr1_bin_coverage quantifies the proportion of non-zero rows
# chr2_bin_coverage quantifies the proportion of non-zero cols
# chr1_row_sums quantifies the total signal value of any row
# chr2_col_sums quantifies the total signal value of any col




# plot data

Brick_vizart_plot_heatmap(File = file.path(tempdir(), 
                                           "chr3R-1-10MB-normal.pdf"),
                          Bricks = list(My_BrickContainer),
                          x_coords = "chr1:54301673:54305898",
                          y_coords = "chr1:54301673:54305898",
                          resolution = 100000,
                          palette = "Reds",
                          width = 10,
                          height = 11,
                          return_object=TRUE)


Failsafe_log10 <- function(x){
   x[is.na(x) | is.nan(x) | is.infinite(x)] <- 0
   return(log10(x+1))
}

Brick_vizart_plot_heatmap(File = file.path(tempdir(),
                                           "chr3R-1-10MB-normal-colours-log10-rotate.pdf"),
                          Bricks = list(My_BrickContainer),
                          x_coords = "chr1:54301673:54305898",
                          y_coords = "chr1:54301673:54305898",
                          resolution = 100000,
                          FUN = Failsafe_log10,
                          value_cap = 0.99,
                          distance = 60,
                          legend_title = "Log10 Hi-C signal",
                          palette = "Reds",
                          width = 10,
                          height = 11,
                          rotate = TRUE,
                          return_object = TRUE)








# EXEMPLE DE LA VIGNETTE

Bintable_path <- system.file(file.path("extdata",
                                       "Bintable_100kb.bins"), package = "HiCBricks")

out_dir <- file.path(tempdir(), "HiCBricks_vignette_test")
dir.create(out_dir)
Create_many_Bricks(BinTable = Bintable_path, 
                   bin_delim=" ", output_directory = out_dir, 
                   file_prefix = "HiCBricks_vignette_test", remove_existing=TRUE, 
                   experiment_name = "HiCBricks vignette test", resolution = 100000)

BrickContainer_dir <- file.path(tempdir(), "HiCBricks_vignette_test")
My_BrickContainer <- load_BrickContainer(project_dir = BrickContainer_dir)
Example_dataset_dir <- system.file("extdata", package = "HiCBricks")

Chromosomes <- c("chr2L", "chr3L", "chr3R", "chrX")
for (chr in Chromosomes) {
   Matrix_file <- file.path(Example_dataset_dir,
                            paste(paste("Sexton2012_yaffetanay_CisTrans_100000_corrected", 
                                        chr, sep = "_"), "txt.gz", sep = "."))
   Brick_load_matrix(Brick = My_BrickContainer,
                     chr1 = chr,
                     chr2 = chr,
                     resolution = 100000,
                     matrix_file = Matrix_file,
                     delim = " ",
                     remove_prior = TRUE)
}






*/