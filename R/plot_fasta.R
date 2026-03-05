#' Plot fasta
#'
#' Plot gene arrow diagram from fasta headers
#'
#' Plots gene arrow diagrams for the locations specified in fasta header
#'
#' @importFrom utils head
#' @param x A \code{data.table} of GWAS hits.
#' @param gff GFF annotation file as \code{\link[data.table:data.table]{data.table}}.
#' @param bed Path to BED file.
#' @param G SNPs in \code{bigsnpr} format or path to BED file.
#' @param threshold Significance threshold for SNPs to include (-log10(p)).
#' @param max_hits Number of SNPs to plot.
#' @param nrow Number of genes to plot per page.
#' @param n_genes Maximum number of genes to plot per SNP.
#' @param title Whether to include titles.
#' @param ... Additional arguments to \code{plot_gwas_single}.
#' @author Ethan Bass
#' @export

plot_gwas <- function(x, gff, bed, G, threshold = 7.5, max_hits = 21,
                      nrow = 7, n_genes = 10, title = TRUE, ...){
  headers <- grep(">", readLines(fasta), value = TRUE)
  headers <- gsub(">LTM","", headers)
  info <- stringr::str_split_fixed(headers, ":", 2)
  info <- data.frame(chr=info[,1], apply(stringr::str_split_fixed(info[,2], "-", 2),2,as.numeric))
  colnames(info)[2:3] <- c("start", "end")
  info$fragment <- seq_len(nrow(info))
  info$seqid <- 1
  plot_gwas_single(genes, "v2.2_chr3", 6321824, n = 5) +
    geom_subgene_arrow(data = example_subgenes,
                       aes(xmin = start, xmax = end, y = molecule, fill = gene,
                           xsubmin = from, xsubmax = to), color="black", alpha=.7)
  # gggenes::geom_gene_arrow(data=info[c(1:2),], aes(xmin = .data$start, xmax = .data$end,
  #                 fill = .data$fragment, y = .data$fragment), inherit.aes = FALSE)
  p_lrt <- NULL # due to NSE notes in R CMD check
  G <- get_G_from_bed(G)
  if (!fs::file_exists(bed)){
    stop()
  }
  x_sel <- clump_snps(x = x, G = G, threshold = threshold)
  x_sel <- head(x_sel[order(p_lrt)], n = max_hits)
  plot_list <- apply(x_sel, 1, function(x){
    p <- plot_gwas_single(gff = gff, chr = x[["chr"]], pos = as.numeric(x[["ps"]]),
                          n = n_genes, ...)
    if (title){
      tb <- "                   "
      p <- p + ggtitle(label = "",
                       subtitle = paste0(x[["rs"]],
                                         tb, "-log10(p) = ",
                                         round(-log10(as.numeric(x[["p_lrt"]])), 2),
                                         "\n", tb, tb, tb, "         ",
                                         "MAF = ", calculate_maf(bed = bed, rs = x[["rs"]])))
    }
  })
  ggpubr::ggarrange(plotlist = plot_list, nrow = nrow, align = "v")
}
