utils::globalVariables(c(
  "LTM"
))


#' Extract gene ranges
#'
#' Summarizes per-gene genomic ranges from a data.table of features, returning
#' the min/max coordinates, strand, and fragment count for each gene.
#'
#' @param x A data frame or \code{data.table} with columns \code{gene},
#'   \code{start}, \code{end}, \code{molecule}, \code{strand}, and \code{sseqid}.
#' @param extend Integer. Number of bases to extend the range in both directions.
#'   Defaults to \code{0}. Cannot result in negative start coordinates.
#'
#' @return A \code{data.table} with one row per gene and columns \code{gene},
#'   \code{start}, \code{end}, \code{molecule}, \code{strand}, \code{sseqid},
#'   and \code{n_fragments}.
#'
#' @export
extract_gene_ranges <- function(x, extend=0){
  if (!inherits(x, "data.table")){
    x <- as.data.table(x)
  }
  x[, list(
    start = max(0, min(start) - extend),
    end = max(end) + extend,
    molecule = first(molecule),
    strand = first(strand),
    sseqid = first(sseqid),
    n_fragments = .N
  ), by = gene]
}

#' Extract CDS ranges
#' @importFrom IRanges IRanges start end
#' @export
extract_cds_ranges <- function(x){
  x |>
    dplyr::group_by(gene_id, sseqid) |>
    dplyr::summarize(
      ranges = list(IRanges(start, end)),
      strand = ifelse(mean(send - sstart) > 0, "+", "-"),  # crude strand guess
      .groups = "drop"
    )
}

#' Extract CDS sequences
#' @importFrom Biostrings DNAStringSet xscat
#' @export
extract_cds_sequence <- function(cds_ranges, genome) {
  cds_parts <- DNAStringSet()
  for (i in seq_len(nrow(cds_ranges))) {
    chr <- cds_ranges$seqid[i]
    gene <- cds_ranges$ID[i]
    start <- cds_ranges$start[i]
    end <- cds_ranges$end[i]

    part <- genome[[grep(chr, names(genome))]][start:end]
    cds_parts <- c(cds_parts, DNAStringSet(part))
  }
  # Concatenate all CDS parts (like splicing)
  full_seq <- do.call(xscat, as.list(cds_parts))

  # Reverse complement if the strand is negative
  if (cds_ranges$strand[i] == "-") {
    full_seq <- reverseComplement(full_seq)
  }
  full_seq
}

#' Extract gene sequences
#' @importFrom Biostrings DNAStringSet reverseComplement subseq
extract_gene_sequences <- function(gene_ranges, genome) {
  sequences <- list()
  for (i in seq_len(nrow(gene_ranges))) {
    chr <- gene_ranges$sseqid[i]
    gene <- gene_ranges$gene[i]
    # irs <- IRanges(gene_ranges$start[i], gene_ranges$end[i])

    # Extract each segment
    cds_parts <- DNAStringSet()  # start empty

    full_seq <- subseq(genome[[chr]], gene_ranges$start[i], gene_ranges$end[i])

    # if (what == "CDS"){
    # for (j in seq_len(nrow(gene_ranges))) {
    #   ir <- irs[j]
    #   part <- genome[[chr]][start(ir):end(ir)]
    #   cds_parts <- c(cds_parts, DNAStringSet(part))
    # }

    # Concatenate all CDS parts (like splicing)
    # full_seq <- do.call(xscat, as.list(cds_parts))
    # }

    # Reverse complement if the strand is negative
    if (gene_ranges$strand[i] == "-") {
      full_seq <- reverseComplement(full_seq)
    }
    sequences[[gene]] <- full_seq
  }
  DNAStringSet(sequences)
}


#' Splice CDS
#'
#' Extracts and splices CDS (or other) sequences from a genome FASTA file given
#' a data.table of genomic features. Handles strand-aware sorting so minus-strand
#' features are concatenated in the correct order.

#' @param dt A \code{data.table} (or coercible object) containing genomic features
#'   with at minimum columns for seqname, start, end, strand, and attributes.
#' @param genome_path Character. Path to an indexed FASTA file.
#' @param ID Character or \code{NULL}. Name of the column in \code{dt} to use as
#'   the feature ID for grouping. If \code{NULL}, the ID is parsed from the
#'   \code{attributes} column.
#' @param seqnames_field Character or \code{NULL}. Name of the column in \code{dt}
#'   containing sequence names. If \code{NULL}, defaults to the column matching
#'   \code{"seqid"}.
#' @param splice Logical. If \code{TRUE} (default), returns spliced sequences as a
#'   \code{DNAStringSet} with one entry per feature ID. If \code{FALSE}, returns
#'   unspliced per-range sequences.
#'
#' @importFrom S4Vectors mcols endoapply
#' @importFrom Biostrings getSeq DNAStringSet
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom Rsamtools FaFile
#'
#' @return A \code{DNAStringSet} of spliced (or unspliced) sequences, named by
#'   feature ID.
#'
#' @export
splice_cds <- function(dt, genome_path, ID=NULL, seqnames_field = NULL, splice=TRUE) {
  dt_proc <- as.data.table(dt)
  if (is.null(seqnames_field)){
    seqnames_field <- grep("seqid", colnames(dt), value = TRUE)
  }
  if (is.null(ID)){
    dt_proc[, ID := sub(".*ID=([^;]+);.*", "\\1", attributes)]
  } else{
    dt_proc[,ID:=dt_proc[[ID]]]
  }
  gr <- makeGRangesFromDataFrame(
    dt_proc,
    keep.extra.columns = TRUE,
    seqnames.field = seqnames_field
  )

  gr_list <- split(gr, mcols(gr)$ID)

  gr_list <- endoapply(gr_list, function(g) {
    if (all(as.character(strand(g)) == "-")) {
      sort(g, decreasing = TRUE)
    } else {
      sort(g)
    }
  })

  genome_ref <- Rsamtools::FaFile(genome_path)

  flat_gr <- unlist(gr_list)

  flat_seqs <- getSeq(genome_ref, flat_gr)

  cds_seqs_list <- split(flat_seqs, names(flat_gr))

  final_cds <- DNAStringSet(lapply(cds_seqs_list, unlist))
  if (splice){
    return(final_cds)
  } else{
    return(flat_seqs)
  }

}
