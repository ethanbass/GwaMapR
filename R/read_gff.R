#' Read GFF3
#' @import data.table
#' @param path Path to GFF3 file.
#' @param what What features to include (\code{gene} or \code{mRNA}).
#' @author Ethan Bass
#' @export

read_gff <- function(path, what = c("gene", "mRNA", "CDS", "exon", "all")){
  type <- NULL # due to NSE notes in R CMD check
  comment_lines <- count_comment_lines(path)
  what <- match.arg(what, c("gene", "mRNA", "CDS", "exon", "all"))
  genes <- data.table::fread(path, skip = comment_lines, sep = "\t", sep2 = ";",
              col.names = c("seqid", "source", "type", "start", "end",
                            "score", "strand", "phase", "attributes"))
  if (what != "all"){
    genes <- genes[type == what]
  }
  genes
}


#' Count number of comment lines
#' @noRd
count_comment_lines <- function(file, pattern = "^#") {
  con <- file(file, "r")
  on.exit(close(con))

  count <- 0
  while(TRUE) {
    line <- readLines(con, n = 1)
    if(length(line) == 0 || !grepl(pattern, line)) {
      break
    }
    count <- count + 1
  }
  return(count)
}
