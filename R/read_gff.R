#' Read GFF3
#' @import data.table
#' @param path Path to GFF3 file
#' @param what What features to include (\code{gene} or \code{mRNA})
#' @author Ethan Bass
#' @export

read_gff <- function(path, what = c("gene", "mRNA")){
  type <- NULL # due to NSE notes in R CMD check
  what <- match.arg(what, c("gene", "mRNA"))
  df <- data.table::fread(path, skip = 1, sep = "\t", sep2 = ";",
              col.names = c("seqid", "source", "type", "start", "end",
                            "score", "strand", "phase", "attributes"))
  genes <- df[type == what]
  genes
}
