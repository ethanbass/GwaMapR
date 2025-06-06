#' Get genes
#'
#' Returns a list of genes associated with a specified SNP based on the GFF3
#' annotation file.
#'
#' @param genes GFF annotation file as \code{\link[data.table:data.table]{data.table}}.
#' @param chr String denoting chromosome
#' @param loc String denoting position on chromosome
#' @param half_width Half width. Defaults to \code{50000}.
#' @author Ethan Bass
#' @export

get_genes <- function(genes, chr, loc, half_width = 50000){
  # Create GRanges object for the single peak position
  peak_point <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = loc, end = loc)
  )

  # Convert gene locations to GRanges
  gene_ranges <- GenomicRanges::GRanges(
    seqnames = genes$seqid,
    ranges = IRanges::IRanges(
      start = genes$start,
      end = genes$end
    ),
    Name = genes$attributes
  )

  # Find distances to peak
  distances <- GenomicRanges::distanceToNearest(gene_ranges[GenomicRanges::seqnames(gene_ranges) == chr],
                                                peak_point)

  # Add distances to data frame and order
  df.local <- genes |>
    dplyr::filter(.data$seqid == chr) |>
    dplyr::mutate(
      distance_to_peak = GenomicRanges::mcols(distances)$distance,
      # Calculate distance from peak to gene center
      distance_to_center = abs((.data$start + .data$end)/2 - loc)
    ) |>
    dplyr::arrange(.data$distance_to_peak) |>
    dplyr::filter(.data$distance_to_peak <= half_width)

  df.local
  # If you want to keep track of whether genes are upstream/downstream
  df.local <- df.local |>
    dplyr::mutate(
      position = dplyr::case_when(
        .data$strand == "+" & loc < .data$start ~ "upstream",
        .data$strand == "+" & loc > .data$end ~ "downstream",
        .data$strand == "-" & loc > .data$end ~ "upstream",
        .data$strand == "-" & loc < .data$start ~ "downstream",
        TRUE ~ "overlapping"
      )
    )
  df.local
}

#' Get genes by index
#' @param genes GFF annotation file as \code{data.table}.
#' @param half_width Half width. Defaults to \code{50000}.
#' @noRd
get_genes_by_idx <- function(df, genes, idx = 1, half_width = 50000){
  df_sel <- head(df[order(df$p_lrt),])
  chr <- df_sel[[idx, "chr"]]
  loc <- df_sel[[idx, "ps"]]

  # Create GRanges object for the single peak position
  peak_point <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = loc, end = loc)
  )

  # Convert gene locations to GRanges
  gene_ranges <- GenomicRanges::GRanges(
    seqnames = genes$seqid,
    ranges = IRanges::IRanges(
      start = genes$start,
      end = genes$end
    ),
    Name = genes$attributes
  )

  # Find distances to peak
  distances <- GenomicRanges::distanceToNearest(gene_ranges[GenomicRanges::seqnames(gene_ranges) == chr],
                                                peak_point)

  # Add distances to data frame and order
  df.local <- genes |>
    dplyr::filter(.data$seqid == chr) |>
    dplyr::mutate(
      distance_to_peak = GenomicRanges::mcols(distances)$distance,
      # Calculate distance from peak to gene center
      distance_to_center = abs((.data$start + .data$end)/2 - loc)
    ) |>
    dplyr::arrange(.data$distance_to_peak) |>
    dplyr::filter(.data$distance_to_peak <= half_width)

  df.local
  # If you want to keep track of whether genes are upstream/downstream
  df.local <- df.local |>
    dplyr::mutate(
      position = dplyr::case_when(
        # For genes on + strand
        .data$strand == "+" & .data$start > loc ~ "downstream",
        .data$strand == "+" & .data$end < loc ~ "upstream",
        # For genes on - strand
        .data$strand == "-" & .data$start > loc ~ "upstream",
        .data$strand == "-" & .data$end < loc ~ "downstream",
        TRUE ~ "overlapping"
      )
    )
  df.local
}

#' Get peaks
#'
#' Returns a list of genes associated with the phenotype and the number of times
#' the gene appears.
#'
#' @importFrom utils head
#' @param df A \code{data.table} of GWAS results from \code{gemma}.
#' @param genes GFF annotation file as \code{data.table}.
#' @param n Number of genes to return.
#' @param t Significance threshold
#' @param half_width Half width. Defaults to \code{1000}.
#' @author Ethan Bass
#' @export

get_peaks <- function(df, genes, n = NULL, t = NULL, half_width = 1000){
  p_lrt <- NULL # due to NSE notes in R CMD check
  peaks_df <- df[order(df$p_lrt)]
  if (!is.null(n)){
    peaks_df <- head(peaks_df, n = n)
  } else if (!is.null(t)){
    peaks_df <- peaks_df[p_lrt < t]
  } else{
    stop("Please specify `n` or `t`.")
  }

  # First create GRanges for your GWAS peaks
  # Assuming you have a data frame of peaks with columns: chr, position, pvalue/LOD
  peaks_gr <- GenomicRanges::GRanges(
    seqnames = peaks_df$chr,
    ranges = IRanges::IRanges(
      start = peaks_df$ps - half_width,
      end = peaks_df$ps + half_width
    ),
    peak_pos = peaks_df$ps,
    score = peaks_df$p_lrt  # or LOD score
  )

  # Create GRanges for genes
  genes_gr <- GenomicRanges::GRanges(
    seqnames = genes$seqid,
    ranges = IRanges::IRanges(
      start = genes$start,
      end = genes$end
    ),
    gene_id = genes$attributes,
    strand = genes$strand
  )

  # Find overlaps between peaks and genes
  overlaps <- GenomicRanges::findOverlaps(peaks_gr, genes_gr)

  # Create a summary of which peaks map to which genes
  peak_gene_mapping <- data.frame(
    peak_position = peaks_gr$peak_pos[S4Vectors::queryHits(overlaps)],
    chr = GenomicRanges::seqnames(peaks_gr)[S4Vectors::queryHits(overlaps)],
    gene_id = genes_gr$gene_id[S4Vectors::subjectHits(overlaps)],
    score = peaks_gr$score[S4Vectors::queryHits(overlaps)]
  ) |>
    # Group by gene to see which ones have multiple peaks
    dplyr::group_by(.data$gene_id) |>
    dplyr::summarize(
      n_peaks = n(),
      peak_positions = paste(.data$peak_position, collapse = ","),
      max_score = max(.data$score),
      chr = .data$chr[[1]]
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_peaks), dplyr::desc(.data$max_score))
  peak_gene_mapping
}
