library(recombClust)

## Packages to load SNP data
library(VariantAnnotation)
library(GenomicRanges)


## Function to load a VCF file, convert it to a snpMatrix and removing SNPs with a low MAF
getVCFmatrixChr <- function(range = NULL, samples = NULL, snps.names = NULL, minmaf = 0.1, Remove.Granges = NULL, ...){
   vcf <- loadVCFrange(range, samples, ...)
   
   snpsVCF <- genotypeToSnpMatrix(vcf)
   snpsVCF$genotypes <- snpsVCF$genotypes[, !snpsVCF$map$ignore]
   sums <- col.summary(snpsVCF$genotypes)
   snps <- colnames(snpsVCF$genotypes)[sums$MAF > minmaf]
   
   vcf <- vcf[snps, ]
   geno <- geno(vcf)$GT
   phase <- lapply(1:ncol(geno), function(x){
      chr1 <- as.numeric(sapply(geno[, x], substring, 1, 1))
      chr2 <- as.numeric(sapply(geno[, x], substring, 3, 3))
      matrix(c(chr1, chr2), nrow = 2, byrow = TRUE)
   })
   phase <- Reduce(function(...) rbind(...), phase)
   rownames(phase) <- paste(rep(colnames(geno), each = 2), 1:2, sep = "_")
   colnames(phase) <- rownames(geno)
   snpsVCF <- list(genotypes = new("SnpMatrix", 2*phase + 1), map = data.frame(name = colnames(phase)))
   ## Conversion from VCF to SNP matrix produces some SNPs to be NA (multiallelic or bigger than 1)
   snpsVCF$map$position <- start(rowRanges(vcf))
   snpsVCF$map$chromosome <- as.character(GenomicRanges::seqnames(rowRanges(vcf)))
   snpsVCF$map$allele.2 <- unlist(snpsVCF$map$allele.2)
   rownames(snpsVCF$map) <- rownames(geno)
   snpsVCF
}

## Load VCF selecting samples 
loadVCFrange <- function(range = NULL, samples = NULL, vcffile){
   vcfsamples <- samples(scanVcfHeader(vcffile))
   if (!is.null(samples)){
      samples <- vcfsamples[vcfsamples %in% samples]
   } else{
      samples <- vcfsamples
   }
   if (!is.null(range)){
      param <- ScanVcfParam(samples = samples, which = range)
   } else {
      param <- ScanVcfParam(samples = samples)
   }
   vcf <- readVcf(vcffile, "hg19", param)
   vcf
}




## Get vcf path
vcf_file <- system.file("extdata", "example.vcf", package = "recombClust")
   
## Load vcf file
snps <- getVCFmatrixChr(vcf = vcf_file)
   
# You should convert them to a GenomicRanges prior passing them to recombClust:
GRsnps <- makeGRangesFromDataFrame(snps$map, start.field = "position", 
                                      end.field = "position")
   
## Create matrix with snps
snpMat <- as(snps$genotypes, "numeric")/2 ## Divide by 2 to have 0 as reference and 1 as alternative
   



# En la parte superior -> obtención datos snpMat y GRsnps necesarios para ejecutar el runRecombClust inicial


res <- microbenchmark( C.file <- CrunRecombClust("inst/extdata/example.vcf",1, 20),
                       R <- runRecombClust(snpMat[, 1:20], annot = GRsnps[1:20]),
                       times = 3L, unit = "s")

print(summary(res)[, c(1:7)],digits=3)

all.equal(C.file$class,R$class )
C.file$class
R$class

C.file$mat[1:5,1:12]
R$mat[1:5,1:12]

all.equal(C.file$mat[1:5,1:12], R$mat[1:5,1:12])
