#' Extract SNPs from BED file by ID
#' @importFrom bigsnpr snp_readBed2 snp_attach
#' @param bed Path to BED file.
#' @param snp_ids Vector of SNP IDs.
#' @author Ethan Bass
#' @export

extract_snps <- function(bed, snp_ids) {
  map <- data.table::fread(fs::path_ext_set(bed, "bim"),
                           col.names = c("chromosome", "marker.ID", "genetic.dist",
                                         "physical.pos", "allele1", "allele2"))

  # Create SNP object from VCF
  snp.idx <- match(snp_ids,map$marker.ID)

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

  if (is.null(names(snp_ids))){
    colnames(result) <- snp_data$map$marker.ID
  } else {
    colnames(result) <- names(snp_ids)
  }

  return(result)
}
