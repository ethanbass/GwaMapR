#' Read GFF3
#' @import data.table
#' @param path Path to GFF3 file.
#' @param what What features to include (\code{gene} or \code{mRNA}).
#' @param expand_attributes Logical. Whether to expand attributes or not.
#' Defaults to \code{FALSE}.
#' @author Ethan Bass
#' @export

read_gff <- function(path, what = c("gene", "mRNA", "CDS", "exon", "all"),
                     expand_attributes = FALSE){
  type <- NULL # due to NSE notes in R CMD check
  comment_lines <- count_comment_lines(path)
  what <- match.arg(what, c("gene", "mRNA", "CDS", "exon", "all"))
  gff <- data.table::fread(path, skip = comment_lines, sep = "\t",
                           col.names = c("seqid", "source", "type", "start",
                                         "end", "score", "strand", "phase",
                                         "attributes"), fill=TRUE)
  if (what != "all"){
    gff <- gff[type == what]
  }
  if (expand_attributes){
    gff <- parse_attributes(gff)
  }
  gff
}

#' Parse attributes
#' @noRd
parse_attributes <- function(dt){
  # add row index to rejoin later
  dt[, .row := .I]

  # split attributes into long format
  attrs_long <- dt[, list(attr = unlist(strsplit(attributes, ";"))), by = .row]

  # remove empty strings from trailing semicolons
  attrs_long <- attrs_long[trimws(attr) != ""]

  # split into key and value
  attrs_long[, c("key", "value") := tstrsplit(trimws(attr), "=",
                                              fixed = TRUE, keep = 1:2)]
  attrs_long[, attr := NULL]

  # pivot wide
  attrs_wide <- dcast(attrs_long, .row ~ key, value.var = "value",
                      fun.aggregate = paste, collapse = ",")

  # join back and clean up
  dt <- dt[attrs_wide, on = ".row"]
  dt[, c(".row", "attributes") := NULL]
  dt[]
  return(dt)
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
