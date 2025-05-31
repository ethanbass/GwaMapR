#' Extract SNPs from BED file by ID
#'
#' Extracts raw SNP data from a \code{bed} file for SNPs specified by
#' \code{snp_ids}.
#'
#' @importFrom bigsnpr snp_readBed2 snp_attach
#' @param bed Path to BED file.
#' @param snp_ids Vector of SNP IDs.
#' @param fam_cols Which "family" columns to include in the output. Options are
#' \code{FID} and/or \code{"IID"}.
#' @param rename Logical. Whether to simplify the names of the SNPs. Defaults to
#' \code{TRUE}.
#' @author Ethan Bass
#' @export

get_snps <- function(bed, snp_ids, fam_cols = "FID", rename = TRUE){
  match.arg(fam_cols, c("FID", "IID"), several.ok = TRUE)
  map <- data.table::fread(fs::path_ext_set(bed, "bim"),
                           col.names = c("chromosome", "marker.ID", "genetic.dist",
                                         "physical.pos", "allele1", "allele2"))
  fam <- data.table::fread(fs::path_ext_set(bed, "fam"),
                           col.names = c("FID", "IID", "MID", "PID", "SEX", "PHENO"))

  # Create SNP object from VCF
  snp.idx <- match(snp_ids, map$marker.ID)

  tmpfile <- tempfile()
  bigsnpr::snp_readBed2(bed, backingfile = tmpfile,
               ind.col = snp.idx)
  snp_data <-  bigsnpr::snp_attach(paste0(tmpfile, ".rds"))

  # Extract genotypes as a matrix
  geno_matrix <- snp_data$genotypes[]

  # Convert to data frame
  result <- as.data.frame(geno_matrix)

  # Add sample IDs as row names
  rownames(result) <- snp_data$fam$sample.ID

  if (!is.null(names(snp_ids))){
    colnames(result) <- names(snp_ids)
  } else if (rename){
    chr <- stringr::str_extract(snp_data$map$marker.ID, "(\\d+)(?=:)")
    colnames(result) <- make.unique(paste0("qtl", chr))
  } else{
    colnames(result) <- snp_data$map$marker.ID
  }
  return(result)
}
