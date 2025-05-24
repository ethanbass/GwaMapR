#' Plot gene arrow diagram from GWAS results
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

#' Calculate MAF
#'
#' Calculate minor allele frequency
#'
#' Uses \code{bigsnpr} to pull the SNP specified by \code{rs} and calculates the
#' minor allele frequency by dividing the number of individuals with the minor
#' minor allele by the total number of individuals.
#'
#' @param bed Path to bed file.
#' @param rs String specifying the target SNP (from column RS in the output from gemma).
calculate_maf <- function(bed, rs){
  snp <- get_snps(bed, rs)
  table(snp)[[2]]/sum(table(snp))
}

#' Get G from BED
#' @noRd
get_G_from_bed <- function(G){
  if (inherits(G, "character")){
    tmpfile <- tempfile()
    bigsnpr::snp_readBed2(G, backingfile = tmpfile)
    snp_data <-  bigsnpr::snp_attach(paste0(tmpfile, ".rds"))
    G <- snp_data$genotypes
  } else if (!inherits(G,  "FBM.code256")){
    stop()
  }
  G
}

#' Clump SNPs
#' @noRd
clump_snps <- function(x, G, threshold = 7.5, thr.r2 = 0.2, size = 100/thr.r2,
                       ncores = 2){
  chr <- as.numeric(gsub(".*([0-9])$", "\\1", x$chr))
  clmp <- bigsnpr::snp_clumping(G = G, infos.chr = chr,
                                S = -log10(x$p_lrt),
                                thr.r2 = thr.r2,
                                size = size,
                                infos.pos = x$ps,
                                exclude = which(-log10(x$p_lrt) < threshold),
                                ncores = ncores)
  x[clmp]
}

#' Read BED
#' @noRd
read_bed <- function(bed){
  tmpfile <- tempfile()
  bigsnpr::snp_readBed2(bed, backingfile = tmpfile)
  snp_data <-  bigsnpr::snp_attach(paste0(tmpfile, ".rds"))
  snp_data$genotypes
}

#' Plot gene arrow diagram for single GWAS hit
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

plot_gwas_single <- function(gff, chr, pos, n = 10, legend_col=2,
                      legend_size = 0.5, feature_width = 1.2){
  gns <- head(get_genes(gff, chr, pos), n = n)
  gns$gene <- stringr::str_split_fixed(gns$attributes, ",", 2)[, 2]
  gns$gene <- ifelse(gns$gene == "", "Unknown", gns$gene)
  gns$orientation <- as.numeric(gns$strand == "+")
  gns$gene <- factor(gns$gene, levels=unique(gns$gene[order(gns$start)]), )
  ggplot(gns, aes(xmin = .data$start, xmax = .data$end, y = .data$seqid,
                  fill = .data$gene, forward = .data$orientation)) +
    gggenes::geom_gene_arrow() +
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
  gns <- head(get_genes(genes, chr, pos), n = 5)
  gns$gene <- stringr::str_split_fixed(gns$attributes, ",", 2)[, 2]
  gns$gene <- ifelse(gns$gene == "", "Unknown", gns$gene)
  gns$orientation <- as.numeric(gns$strand == "+")
  gns
}
