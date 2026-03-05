utils::globalVariables(c(
  "sstart", "send", "gene_id", "molecule", "gene",
  "start", "end", "strand", "pident", "evalue",
  "bitscore", "sseqid", "qstart", "qend", "qoverlap", "gap",
  "sstart_trim", "cluster_size", "center","start_scaled","end_scaled"
))

#' Read blast results
#'
#' Read blast results from tabular (format 7) files.
#'
#' @param path Path to blast (format 7) output file.
#' @importFrom data.table fread
#' @export
read_blast <- function(path){
  data.table::fread(path, col.names = c("qseqid", "sseqid", "pident", "length",
                            "mismatch", "gapopen", "qstart", "qend",
                            "sstart", "send", "evalue", "bitscore"))
}

#' Cluster blast
#'
#' Cluster blast results
#'
#' This is useful for grouping together exons that belong to the same gene.
#'
#' @param df Output from \code{\link{read_blast}}.
#' @param max_gap Maximum gap between hits for clustering. Defaults to
#' \code{5000}.
#' @param min_elements The mimimum number of elements in a cluster. Clusters
#' with fewer elements will be filtered out. Defaults to \code{1}.
cluster_blast <- function(df, max_gap = 5000, min_elements = 1) {
  df |>
    dplyr::mutate(
      start = pmin(sstart, send),
      end   = pmax(sstart, send),
      strand = ifelse(sstart < send, "+", "-")
    ) |>
    dplyr::arrange(sseqid, start) |>  # sort by genomic start for clustering
    dplyr::group_by(sseqid) |>
    dplyr::mutate(
      gap = start - dplyr::lag(end, default = -1L),
      gene = cumsum(gap > max_gap | dplyr::row_number() == 1)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(gene_id = paste0(sseqid, "_", gene)) |>
    # this second part is necessary for trimming
    dplyr::group_by(gene_id) |>
    dplyr::arrange(qstart, .by_group = TRUE) |>  # sort by query start for overlap trimming
    dplyr::mutate(
      qoverlap = pmax(0L, dplyr::lag(qend, default = -1L) - qstart + 1L),
      qstart_trim = qstart + qoverlap,
      sstart_trim = dplyr::case_when(
        strand == "+" ~ sstart + qoverlap,
        strand == "-" ~ sstart - qoverlap,
        TRUE ~ sstart
      )
    ) |>
    dplyr::mutate(cluster_size = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::filter(cluster_size >= min_elements) |>
    dplyr::select(-cluster_size) |>
    dplyr::mutate(
      start = pmin(sstart_trim, send),
      end   = pmax(sstart_trim, send),
      strand = ifelse(sstart_trim < send, "+", "-")
    ) |> dplyr::arrange(sseqid, start)
}

#' Prepare for gggenes
#' @noRd
prepare_for_gggenes <- function(blast_clustered, filter = FALSE) {
  df <- blast_clustered |>
    dplyr::mutate(
      start = pmin(sstart, send),
      end = pmax(sstart, send),
      strand = ifelse(sstart < send, 1, -1),
      # Use gene_id for coloring by cluster
      gene = gene_id
    ) |> dplyr::rename(molecule = gene_id)
  if (filter){
    df <- dplyr::select(df, molecule , gene, start, end, strand,
                  pident, evalue, bitscore, sseqid)
  }
  df
}

#' Plot genes
plot_all_genes <- function(gggenes_data, title = NULL, fill = "gene") {
  # Calculate max span for any gene cluster
  gene_spans <- gggenes_data |>
    dplyr::group_by(molecule) |>
    dplyr::summarize(span = max(end) - min(start), .groups = "drop")
  max_span <- max(gene_spans$span)

  # Center each gene cluster and give it the same width
  gggenes_centered <- gggenes_data |>
    dplyr::group_by(molecule) |>
    dplyr::mutate(center = (min(start) + max(end)) / 2,
           start_scaled = start - center,
           end_scaled = end - center) |>
    dplyr::ungroup()

  ggplot(gggenes_centered,
         aes(xmin = start_scaled, xmax = end_scaled, y = molecule,
             fill = .data[[fill]], forward = strand > 0)) +
    gggenes::geom_gene_arrow() +
    ggplot2::facet_wrap(~molecule, scales = "free_y", ncol = 1) +
    ggplot2::xlim(-max_span/2, max_span/2) +
    gggenes::theme_genes() +
    ggplot2::labs(title = title,
         x = "Position relative to center (bp)",
         y = "Gene Cluster")
}

#' Blast primers
#'
#' Blast primers against reference sequence
#'
#' @param primer Path to fasta file with primer sequences.
#' @param ref Path to fasta file with reference.
#' @param path_out Path to directory to export files.
#' @param format Blast format. Defaults to \code{7}.
blast_primers <- function(primer, ref, path_out = NULL, format = 7){
  primer_base <- fs::path_file(fs::path_ext_remove(primer))
  if (is.null(path_out)){
    path_out <- fs::path_dir(primer)
  }
  ref_base <- substr(basename(ref),1,3)
  file_out <- fs::path(path_out, paste(primer_base, ref_base,
                                       sep = "_"), ext = "out")
  command <- sprintf('blastn -task blastn-short -query "%s" -db "%s" -outfmt %d -dust no | sed \'/^#/d\' > "%s"',
                     primer, ref, format, file_out)
  system(command)
  return(invisible(file_out))
}

#' Blast gene
#' @export
blast_gene <- function(seq, ref, path_out = NULL, format = 7){
  seq_base <- fs::path_file(fs::path_ext_remove(seq))
  if (is.null(path_out)){
    path_out <- fs::path_dir(seq)
  }
  path_out <- fs::path_expand(fs::path_abs(path_out))
  if (!fs::dir_exists(path_out))
      stop("Directory not found. Please check path and try again.")
  ref_base <- substr(basename(ref), 1, 3)
  file_out <- fs::path(path_out, paste(seq_base, ref_base, sep="_"),
                       ext = sprintf("out%d", format))
  command <- sprintf('blastn -query %s -db %s -out %s -outfmt %d',
                     seq, ref, file_out, format)
  system(command)
  return(invisible(file_out))
}

#' Blast sequence
#' @export
blast_sequence <- function(seq, ref, what = c("gene", "primer"), format = 7,
                           path_out = NULL, seq_name = NULL, ref_name = NULL){
  what <- match.arg(what, c("gene", "primer"))
  if (is.null(path_out)){
    path_out <- fs::path_dir(seq)
  }
  if (is.null(seq_name)){
    seq_name <- fs::path_file(fs::path_ext_remove(seq))
  }
  if (is.null(ref_name)){
    ref_name <- substr(basename(ref), 1, 3)
  }
  file_out <- fs::path(path_out, paste(seq_name, ref_name, sep = "_"),
                       ext = sprintf("out%d",format))
  if (what == "gene"){
    command <- sprintf('blastn -query "%s" -db "%s" -outfmt %d | sed \'/^#/d\' > "%s"',
                       seq, ref, format, file_out)
  } else if (what == "primer"){
    command <- sprintf('blastn -task blastn-short -query "%s" -db "%s" -outfmt %d -dust no | sed \'/^#/d\' > "%s"',
                       seq, ref, format, file_out)
  }
  system(command)
  return(invisible(file_out))
}

#' Plot primer conservation
#' @importFrom Biostrings DNAString reverseComplement subseq start end
#' @importFrom ggplot2 ggplot aes geom_area geom_line theme_minimal labs ylim ggtitle
plot_conservation_primers <- function(primerF, primerR, seqs, conservation){
  f_binding <- find_seq(reverseComplement(DNAString(primerF)), seqs)
  r_binding <- find_seq(primerR, seqs)

  subseqs <- subseq(seqs, start(r_binding), end(r_binding))

  logo_R <- ggmsa::seqlogo(subseqs) + ggtitle("Primer R")
  cons_df <- data.frame(position = start(r_binding):end(r_binding),
                        conservation = conservation[start(r_binding):end(r_binding)])
  p_cons_R <- ggplot(cons_df, aes(x = .data$position, y = .data$conservation)) +
    geom_area(fill = "darkblue", alpha = 0.5) +
    geom_line(color = "darkblue", size = 1) +
    theme_minimal() +
    labs(x = "Position", y = "Conservation") + ylim(c(0, 1))

  subseqs <- subseq(seqs, start(f_binding), end(f_binding))

  logo_F <- ggmsa::seqlogo(subseqs) + ggtitle("Primer F")

  cons_df <- data.frame(position = start(f_binding):end(f_binding),
                        conservation = conservation[start(f_binding):end(f_binding)])

  p_cons_F <- ggplot(cons_df, aes(x = .data$position, y = .data$conservation)) +
    geom_area(fill = "darkblue", alpha = 0.5) +
    geom_line(color = "darkblue", size = 1) +
    theme_minimal() +
    labs(x = "Position", y = "Conservation") + ylim(c(0,1))
  ggpubr::ggarrange((logo_F/p_cons_F), (logo_R/p_cons_R), ncol = 2)
}

#' Plot conservation
#' @importFrom Biostrings DNAString reverseComplement
#' @importFrom ggplot2 ggplot aes geom_area geom_line theme_minimal labs ylim ggtitle
plot_conservation <- function(seq, seqs, conservation, rc = FALSE, title = NULL){
  seq_pos <- find_seq(seq, seqs, rc = rc)
  subseqs <- subseq(conservation$seqs, start(seq_pos), end(seq_pos))
  p_logo <- ggmsa::seqlogo(subseqs) + ggtitle(title)
  seq_range <- start(seq_pos):end(seq_pos)
  cons_df <- data.frame(position = seq_range,
                        conservation = conservation$conservation[seq_range])
  p_conservation <- ggplot(cons_df, aes(x = .data$position, y = .data$conservation)) +
    geom_area(fill = "darkblue", alpha = 0.5) +
    geom_line(color = "darkblue", size = 1) +
    theme_minimal() +
    labs(x = "Position", y = "Conservation") + ylim(c(0,1))
  patchwork::wrap_plots(p_logo, p_conservation, nrow=2)
}

#' Find sequence
#' @importFrom Biostrings DNAString reverseComplement vmatchPattern
find_seq <- function(seq, seqs, rc = FALSE){
  if (inherits(seq, "character")){
    seq <- DNAString(seq)
  }
  if (rc){
    seq <- reverseComplement(seq)
  }
  match_found <- NULL
  matched_seq_index <- NULL

  for (j in 1:length(seqs)) {
    matches <- vmatchPattern(seq,
                             seqs[j]
    )

    if (length(matches) > 0 && length(matches[[1]]) > 0) {
      match_found <- matches
      matched_seq_index <- j
      break  # Stop searching once we find a match
    }
  }
  matches[[1]]
}

#' Plot blast results
#' @param path Path to blast results.
#' @param max_gap Max gap between blast hits for clustering. Defaults to
#' \code{100000}.
plot_blast <- function(path, max_gap = 100000){
  df <- read_blast(path)
  blast_clustered <- cluster_blast(df, max_gap = max_gap)
  gggenes_data <- prepare_for_gggenes(blast_clustered)
  plot_all_genes(gggenes_data)
}
