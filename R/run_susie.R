#' Function to run Sum of Single Effects (SuSiE) Regression
#'
#' Wrapper for Sum of Single Effects (SuSiE) regression in the
#' \code{\link[susieR:susie]{susieR}} package.
#'
#' @importFrom susieR susie susie_plot
#' @importFrom  bigsnpr snp_fastImputeSimple snp_readBed2 snp_attach
#' @param Y Vector of observed responses.
#' @param bed Path to BED annotation file.
#' @param genes Path to GFF3 annotations (or a \code{data.table} of annotations).
#' @param chr The chromosome that the SNP is on.
#' @param locus The position of the SNP on chromosome \code{chr}.
#' @param window Half-width window around \code{locus} for inclusion of SNPs.
#' @param coverage Coverage, argument to \code{\link[susieR:susie]{susie}}.
#' @param plot_it Logical. Whether to plot the results or not. Defaults to
#' \code{TRUE}.
#' @param verbose Logical. Whether to print messages to console.
#' @param ... Additional arguments to \code{\link[susieR]{susieR}}.
#' @author Ethan Bass
#' @export

run_susie <- function(Y, bed, genes, chr, locus, window = 1e6,
                          coverage = 0.95, plot_it = TRUE,
                          verbose = getOption("verbose"), ...){
  genes <- parse_annotations(genes)
  if (!inherits(genes, "data.frame"))
    stop("Genes object does not exist", immediate. = TRUE)
  map <- data.table::fread(fs::path_ext_set(bed, "bim"),
               col.names = c("chromosome", "marker.ID", "genetic.dist",
                             "physical.pos", "allele1", "allele2"))
  snp_idx <- with(map, which(chromosome == chr &
                               physical.pos >= (locus - window) &
                               physical.pos <= (locus + window)))

  samples_idx <- if(any(is.na(Y))){
    which(!is.na(Y))
  } else{
    seq_along(Y)
  }

  tmpfile <- tempfile()
  bigsnpr::snp_readBed2(bed, backingfile = tmpfile,
               ind.row = samples_idx,
               ind.col = snp_idx)
  obj <-  bigsnpr::snp_attach(paste0(tmpfile, ".rds"))

  G   <- obj$genotypes
  Gfi <- snp_fastImputeSimple(G)
  if (any(is.na(Y))){
    Y <- Y[samples_idx]
  }

  fitted <- susieR::susie(as.matrix(Gfi[]), Y, coverage = coverage,
                  verbose = verbose, ...)

  if (plot_it){
    susieR::susie_plot(fitted, y = "PIP", add_legend = TRUE, add_bar = TRUE)
  }

  map_sel <- map[snp_idx]
  cs <- fitted$sets$cs

  annotations <- lapply(cs, function(L){
    lapply(L, function(i){
      get_genes(genes, map_sel[[i, "chromosome"]],
                       loc = map_sel[[i, "physical.pos"]])
    })
  })

  snps <- lapply(cs, function(L){
    map_sel[L]
  })

  list(susie = fitted, snps = snps, annotations = annotations)
}

parse_annotations <- function(genes){
  if (inherits(genes, "character")){
    exists <- fs::file_exists(genes)
    if (!exists){
      stop(sprintf("File %s could not be found", basename(genes)))
    } else{
      genes <- data.table::fread(genes)
    }
  }
}
