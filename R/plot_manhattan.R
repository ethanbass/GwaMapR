#' Manhattan plot
#'
#' Wrapper around the \code{fastman_gg} function in the
#' \code{\link[fastman:fastman]{fastman}} package.
#'
#' @importFrom ggplot2 ggtitle ylab theme
#' @param x A \code{data.frame} containing the GWAS results.
#' @param chr A string denoting the name of the chromosome column in \code{x}.
#' Defaults to \code{"chr"} corresponding to the output from gemma. Argument to
#' \code{fastman_gg}.
#' @param bp #' @param chr A string denoting the column name for the chromosomal
#' position in \code{x}. Defaults to \code{"ps"} corresponding to the output
#' from gemma. Argument to \code{fastman_gg}.
#' @param p #' @param chr A string denoting the column name for the p-value in
#' \code{x}. Defaults to \code{"p_lrt"} corresponding to the output from gemma.
#' Argument to \code{fastman_gg}.
#' @param n The number of SNPs for calculation of Bonferroni threshold. Defaults to
#' \code{nrow(x)}.
#' @param title A string denoting the title to be printed at the top of the plot.
#' @param ... Additional arguments to \code{\link[fastman:fastman]{fastman::fastman_gg}}.
plot_manhattan <- function(x, chr = "chr", bp = "ps",
                           p = "p_lrt", n = NULL, title = NULL, ...){
  check_packages(c("fastman","mdthemes"))
  if (is.null(n)){
    n <- nrow(x)
  }
  fastman::fastman_gg(x, chr = chr, bp = bp, p = p,
                      genomewideline = -log10( .05 /n),
                      suggestiveline = -log10(1/n), ...) +
    ggtitle(ifelse(!is.null(title), "BC-ratio", "")) +
    mdthemes::md_theme_classic() +
    ylab("-log<sub>10</sub>(p)") + theme(legend.position = "none")
}
