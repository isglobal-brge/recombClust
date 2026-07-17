## hdf5omics — embedded, extraction-ready omics layer (see docs/DESIGN.md 6.3).
## Ingestion of phased genomic data into an HDF5 haplotype matrix. Depends only on
## SeqArray/gdsfmt (parsing) and BigDataStatMeth (HDF5); never on recombClust internals.

#' Convert a genomic file to GDS if needed (internal).
#' @param file Path to a .vcf(.gz), .bed or .gds file.
#' @return Path to a GDS file (a temp file for vcf/bed, or \code{file} if already gds).
#' @keywords internal
#' @noRd
.omics_to_gds <- function(file) {
  if (grepl("\\.gds$", file, ignore.case = TRUE)) return(file)
  out <- tempfile(fileext = ".gds")
  if (grepl("\\.vcf(\\.gz)?$", file, ignore.case = TRUE)) {
    SeqArray::seqVCF2GDS(file, out, verbose = FALSE)
  } else if (grepl("\\.bed$", file, ignore.case = TRUE)) {
    SeqArray::seqBED2GDS(file, out, verbose = FALSE)
  } else {
    stop("hdf5omics_import: unsupported file type: ", file)
  }
  out
}

#' Import phased genomic data into an HDF5 haplotype matrix
#'
#' Reads a VCF/GDS/BED file with SeqArray, filters by region, samples and (optionally)
#' biallelic SNVs, flattens the phased genotypes into a haplotype matrix
#' (sample-alleles x SNP, 0/1) and writes it to an HDF5 file via BigDataStatMeth. The
#' returned descriptor carries the SNP positions needed downstream (e.g. by
#' \code{\link{recombClust}}).
#'
#' This is the ingestion of the embedded \code{hdf5omics} layer; MAF filtering and
#' imputation are provided separately (reusing BigDataStatMeth's omics functions).
#'
#' @param file Input .vcf(.gz), .bed or .gds file.
#' @param out_file Output HDF5 file to create/write.
#' @param out_group,out_dataset Location of the haplotype matrix inside \code{out_file}.
#' @param range Optional list(chr, start, end) restricting to a genomic region.
#' @param samples Optional character vector of sample ids to keep.
#' @param snv_only Keep only biallelic single-nucleotide variants (default TRUE).
#' @param minmaf Optional minor-allele-frequency threshold. Applied on the genotypes
#'   with SeqArray (\code{seqSetFilterCond}), which is the correct place for MAF:
#'   it is a per-variant genotype quantity and filtering here keeps the SNP positions
#'   in sync automatically. (BigDataStatMeth's \code{filter_maf} assumes 0/1/2 dosage
#'   coding and is not applicable to the 0/1 haplotype matrix; see docs/DESIGN.md 6.3.)
#' @param overwrite Overwrite an existing output dataset (default TRUE).
#' @return A list descriptor: \code{file}, \code{group}, \code{dataset},
#'   \code{positions}, \code{chromosome}, \code{samples} (haplotype row names),
#'   \code{nchrom} (= 2 x samples), \code{nsnp}.
#' @export
hdf5omics_import <- function(file, out_file,
                             out_group = "haplotypes", out_dataset = "H",
                             range = NULL, samples = NULL,
                             snv_only = TRUE, minmaf = NULL, overwrite = TRUE) {
  gds_file <- .omics_to_gds(file)
  f <- SeqArray::seqOpen(gds_file, readonly = TRUE)
  on.exit(SeqArray::seqClose(f), add = TRUE)

  ## Full-variant chromosome via gdsfmt: SeqArray's "chromosome" accessor segfaults
  ## in some SeqArray/gdsfmt builds; gdsfmt reads the node directly and safely.
  chr_all <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(f, "chromosome"))
  pos_all <- SeqArray::seqGetData(f, "position")

  if (!is.null(samples))
    SeqArray::seqSetFilter(f, sample.id = samples, verbose = FALSE)

  if (!is.null(range)) {
    sel <- chr_all == as.character(range$chr) &
           pos_all >= range$start & pos_all <= range$end
    SeqArray::seqSetFilter(f, variant.sel = sel, verbose = FALSE)
  }

  if (isTRUE(snv_only)) {
    parts <- strsplit(SeqArray::seqGetData(f, "allele"), ",", fixed = TRUE)
    is_snv <- vapply(parts, function(a) length(a) == 2 && all(nchar(a) == 1), logical(1))
    vid <- SeqArray::seqGetData(f, "variant.id")
    SeqArray::seqSetFilter(f, variant.id = vid[is_snv], verbose = FALSE)
  }

  ## MAF filtering on the genotypes (correct per-variant MAF; keeps positions in
  ## sync). NOT BigDataStatMeth::filter_maf, which assumes 0/1/2 dosage coding.
  if (!is.null(minmaf))
    SeqArray::seqSetFilterCond(f, maf = minmaf, verbose = FALSE)

  ## Extract phased genotypes and flatten to sample-allele x SNP.
  ## seqGetData("genotype") is [allele, sample, variant] (column-major); a plain
  ## matrix() reshape to 2*nSample rows reproduces CTransformtoSampleAlleles exactly
  ## (row 2*j+i = sample j, allele i), so no C++ is needed for the flatten.
  geno <- SeqArray::seqGetData(f, "genotype")
  d <- dim(geno); nS <- d[2]; nV <- d[3]
  haplos <- matrix(as.double(geno), nrow = 2L * nS, ncol = nV)

  sid <- SeqArray::seqGetData(f, "sample.id")
  rn <- as.vector(rbind(paste0(sid, "_1"), paste0(sid, "_2")))
  pos <- SeqArray::seqGetData(f, "position")
  vid2 <- SeqArray::seqGetData(f, "variant.id")
  chr <- chr_all[vid2]

  BigDataStatMeth::hdf5_create_matrix(
    out_file, paste0(out_group, "/", out_dataset),
    data = haplos, dtype = "double", overwrite = overwrite)
  BigDataStatMeth::hdf5_close_all()

  list(file = out_file, group = out_group, dataset = out_dataset,
       positions = pos, chromosome = chr, samples = rn,
       nchrom = 2L * nS, nsnp = nV)
}
