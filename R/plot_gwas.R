#' Plot gene arrows
#' @importFrom gggenes geom_gene_arrow geom_feature
#' @importFrom ggplot2 ggplot aes guides ylab guide_legend theme unit ggtitle
#' ylab .data scale_fill_brewer
#' @importFrom utils head
#' @param genes GFF annotation file as data.table.
#' @param chr String denoting chromosome.
#' @param pos String denoting position on chromosome.
#' @param legend_size Legend size in cm. Defaults to \code{0.5}.
#' @param feature_size Defaults to \code{1.2}.

plot_gwas <- function(genes, chr, pos, n = 10, legend_col = 2,
                      legend_size = 0.5, feature_size = 1.2){
  gns <- head(get_genes_by_pos(genes, chr, pos), n = n)
  gns$gene <- stringr::str_split_fixed(gns$attributes, ",", 2)[,2]
  gns$gene <- ifelse(gns$gene == "", "Unknown", gns$gene)
  gns$orientation <- as.numeric(gns$strand == "+")
  gns$gene <- factor(gns$gene, levels = unique(gns$gene[order(gns$start)]))
  ggplot2::ggplot(gns, aes(xmin = .data$start, xmax = .data$end, y = .data$seqid,
                  fill = .data$gene, forward = .data$orientation)) +
    gggenes::geom_gene_arrow() +
    scale_fill_brewer(palette = "Set3") + gggenes::theme_genes() +
    gggenes::geom_feature(
      aes(x = pos, y = chr), linewidth = 2
    ) + ylab("") +
    guides(fill = guide_legend(ncol = legend_col)) +
    theme(legend.key.size = unit(legend_size, 'cm'))
}

#' Convert GFF (general feature format) to 'gggenes' format
#' @importFrom utils head
#' @param genes GFF annotation file as data.table
#' @param chr Chromosome
#' @param pos Position on chromosome
gff_2gggenes <- function(genes, chr, pos){
  gns <- head(get_genes_by_pos(genes, chr, pos), n=5)
  gns$gene <- stringr::str_split_fixed(gns$attributes, ",", 2)[,2]
  gns$gene <- ifelse(gns$gene == "", "Unknown", gns$gene)
  gns$orientation <- as.numeric(gns$strand == "+")
  gns
}
