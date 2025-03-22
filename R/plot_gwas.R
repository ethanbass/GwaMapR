#' Plot gene arrow diagram from GWAS results
#' @importFrom gggenes geom_gene_arrow geom_feature theme_genes
#' @importFrom ggplot2 ggplot aes guides ylab guide_legend theme unit ggtitle
#' ylab .data scale_fill_brewer
#' @importFrom utils head
#' @param gff GFF annotation file as \code{\link[data.table:data.table]{data.table}}.
#' @param chr Chromosome or contig to plot.
#' @param pos Position in chromosome or contig to plot.
#' @param n Number of genes to include. This will be the top \code{n} genes
#' closest to the provided SNP.
#' @param legend_col Number of columns to include in legend.
#' @param legend_size Size of legend (in cm). Defaults to \code{0.5}.
#' @param feature_width Line width for feature annotations. Defaults to \code{1.2}.
#' @author Ethan Bass
#' @export

plot_gwas <- function(gff, chr, pos, n = 10, legend_col=2,
                      legend_size = 0.5, feature_width = 1.2){
  gns <- head(get_genes_by_pos(gff, chr, pos), n = n)
  gns$gene <- stringr::str_split_fixed(gns$attributes, ",", 2)[, 2]
  gns$gene <- ifelse(gns$gene == "", "Unknown", gns$gene)
  gns$orientation <- as.numeric(gns$strand == "+")
  gns$gene <- factor(gns$gene, levels=unique(gns$gene[order(gns$start)]), )
  ggplot(gns, aes(xmin = .data$start, xmax = .data$end, y = .data$seqid,
                  fill = .data$gene, forward = .data$orientation)) +
    gggff::geom_gene_arrow() +
    scale_fill_brewer(palette = "Set3") + gggenes::theme_genes() +
    gggenes::geom_feature(
      aes(x = pos, y = chr), linewidth = feature_width
    ) + ylab("") +
    guides(fill = guide_legend(ncol = legend_col)) +
    theme(legend.key.size = unit(legend_size, 'cm'))
}

#' Convert GFF to 'gggenes' format
#' @importFrom utils head
#' @param gff GFF annotation file as \code{\link[data.table:data.table]{data.table}}.
#' @param chr Chromosome or contig to convert.
#' @param pos Position in chromosome or contig to convert.
#' @noRd
gff_2gggenes <- function(genes, chr, pos){
  gns <- head(get_genes_by_pos(genes, chr, pos), n = 5)
  gns$gene <- stringr::str_split_fixed(gns$attributes, ",", 2)[, 2]
  gns$gene <- ifelse(gns$gene == "", "Unknown", gns$gene)
  gns$orientation <- as.numeric(gns$strand == "+")
  gns
}
