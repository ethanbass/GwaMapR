utils::globalVariables(c(
  "type"
))

#' Plot genes
#' Plot genes from GFF annotations
#' @export
plot_genes <- function(gff, what = c('exons', 'genes'),
                       legend_col = 2, legend_size = 0.5,
                       feature_width = 1.2){
  match.arg(what, c('exons', 'genes'))
  gff$gene <- stringr::str_split_fixed(gff$attributes, ",", 2)[, 2]
  gff$gene <- ifelse(gff$gene == "", "Unknown", gff$gene)
  gff$orientation <- as.numeric(gff$strand == "+")
  gff$gene <- factor(gff$gene, levels=unique(gff$gene[order(gff$start)]), )
  ggplot(gff, aes(xmin = .data$start, xmax = .data$end, y = .data$seqid,
                  fill = .data$gene, forward = .data$orientation)) +
    gggenes::geom_gene_arrow() +
    scale_fill_brewer(palette = "Set3") + gggenes::theme_genes() +
    ylab("") +
    guides(fill = guide_legend(ncol = legend_col)) +
    theme(legend.key.size = unit(legend_size, 'cm'))
}

#' Split attributes
#' @keywords internal
#' @noRd
split_attributes <- function(x){
  stringr::str_split_fixed(x$attributes, ";", n=5)
}

#' Get gene ID
#' @noRd
get_id <- function(x){
  ID <- grep("^ID=", split_attributes(x), value = TRUE)
  gsub("ID=", "", ID)
}

#' Get mRNA
#' Get mRNA from gene parent
#' @export
get_mRNA <- function(x, gff){
  if (all(x$type != "gene")){
    stop("Input must be a gene")
  }
  ID <- get_id(x)
  mRNA <- gff[type == "mRNA"]
  parents <- gsub("Parent=", "", split_attributes(mRNA)[,2])
  idx <- which(parents %in% ID)
  mRNA[idx]
}

#' Get exons from parent
#' @export
get_children <- function(x, gff, what = c("CDS","exon")){
  x_cp <- as.data.table(x)
  what <- match.arg(what, c("CDS","exon"))
  if (all(x_cp$type == "gene")){
    x <- get_mRNA(x_cp, gff)
  }
  ID <- get_id(x)
  cld <- gff[type == what]
  parents <- gsub("Parent=","", split_attributes(cld)[,2])
  idx <- which(parents %in% ID)
  cld[idx]
}

#' Plot exons
#' @export
plot_exons <- function(x, gff){
  # mrna <- get_mRNA(x,gff)
  exons <- get_exons(x, gff)
  plot_genes(exons)
}

#' Extract CDS
#' @export
extract_CDS <- function(x, gff, genome_path){
  cds <- get_children(x, gff, what = "CDS")
  splice_cds(exons, genome_path = genome_path)
}
