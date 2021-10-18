#' Main function to run all the recombClust pipeline
#'
#' This function takes as input phased chromosomes and returns the responsibilities
#' PCA and the cluster classification of the samples.
#'
#' @export
#'
#' @param filename File with genotypic data (SNPs and samples) the data file type must be .vcf (Variant Call Format),
#' .gds (Genomic Data Structure) or .bed (Browser Extensible Data).
#' @param annot GenomicRanges with the SNPs annotation
#' @param range Region GenomicRanges
#' @param resfilename string, path and file name to store results. By default
#' this runLDmixture stores data in recombClustResults.hdf5, if resfilename 
#' is empty then function works on memory. To work with big regions it is recommended 
#' to work with data on disk. Only use resfilename = "" when regions are small in order to 
#' prevent memory overflows and low performance.
#' @param resgroup string, folder inside a file where we want to store results. By default
#' runLDmixture stores data inside group "AllModels" in resfilename file.
#' @param overwrite boolean, (optional) either a logical value indicating whether
#'  the output file can be overwritten or not, by default files are not overwritten.
#' @param window only to debug, by default window = 500bp
#' @param ... Further arguments passed to runLDmixtureModel
#' @return probTab: Summarized probabilities by window
computeRecombProbabilities <- function(haplos, annot, range, 
                                       resfilename = "recombClustResults.hdf5",
                                       resgroup = "AllModels",
                                       overwrite = FALSE,
                                       window = 500, ...) 
{

    # Get models
    models <- runLDmixtureModel(haplos, annot, range, resfilename, resgroup, overwrite, ...)
    
    ## Remove failed models
    goodModels <- vapply(models, class, character(1)) == "list"

    if(length(goodModels)>0) {
        
        if(resfilename == "") {
            # Get matrix
            matProb <- createMatrixfromModels(models[goodModels], range)
            # Summarize recombClust probabilities by windows
            matWindProb <- getRecombProb(matProb$mat, matProb$annot, range, window)
            rownames(matWindProb) <- rownames(haplos)
            
        } else {
            getRecombProb_hdf5(resfilename, resgroup, range, window)
        }
        
    } else {
        matWindProb <- NA
    }
        
    
    
    # # Return models and matrix with summarized probabilities by widow
    return(probTab = matWindProb)
}
