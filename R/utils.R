#' Check for suggested package
#'
#' This function checks for a suggested package and returns an error if the
#' package is not installed (if \code{return_boolean} is FALSE. Otherwise, it
#' returns a boolean value.
#'
#' @noRd

check_for_pkg <- function(pkg, return_boolean = FALSE,
                          source = c("cran", "bioconductor")){
  source <- match.arg(source, c("cran", "bioconductor"))
  pkg_exists <- requireNamespace(pkg, quietly = TRUE)
  if (return_boolean){
    return(pkg_exists)
  } else if (!pkg_exists) {
    suggestion <- switch(source, "cran" = sprintf("`install.packages('%s')`", pkg),
                         "bioconductor" = sprintf("`BiocManager::install('%s')`", pkg))
    stop(sprintf("Package `%s` must be installed to perform this action:
                 try `%s`", pkg, suggestion))
    stop(paste(
      "Package", sQuote(pkg), "must be installed to perform this action:
          try", paste0("`install.packages('", pkg, "')`.")),
      call. = FALSE
    )
  }
  invisible(pkg_exists)
}

#' Check packages
#' @noRd
check_packages <- function(pkg){
  invisible(sapply(pkg, check_for_pkg))
}
