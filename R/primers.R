
parse_primer3 <- function(file_path) {

  # Read the file
  lines <- readLines(file_path)

  # Split into individual primer pairs (separated by "=")
  pair_indices <- which(lines == "=")

  # Initialize list to store results
  results_list <- list()

  # Parse each primer pair
  start_idx <- 1
  for (i in seq_along(pair_indices)) {
    end_idx <- pair_indices[i] - 1
    pair_lines <- lines[start_idx:end_idx]

    # Parse key-value pairs
    pair_data <- list()
    for (line in pair_lines) {
      if (grepl("=", line)) {
        parts <- strsplit(line, "=", fixed = TRUE)[[1]]
        key <- parts[1]
        value <- ifelse(length(parts) > 1, parts[2], "")
        pair_data[[key]] <- value
      }
    }

    # Add to results if we have data
    if (length(pair_data) > 0) {
      results_list[[i]] <- pair_data
    }

    # Update start index for next pair
    start_idx <- pair_indices[i] + 1
  }

  # Convert to data frame
  df <- dplyr::bind_rows(results_list)

  # Convert numeric columns
  numeric_cols <- grep("TM|GC_PERCENT|PENALTY|TH|STABILITY|SIZE",
                       names(df), value = TRUE)
  df <- df |>
    dplyr::mutate(dplyr::across(dplyr::all_of(numeric_cols), as.numeric))

  return(df)
}

write_primer3_batch_file <- function(pairs, id = NULL, seqF = NULL, seqR = NULL,
                                     path_out = NULL){
  if (is.null(id)){
    pairs$ID <- seq_len(nrow(pairs))
  }
  if (is.null(path_out)){
    path_out <- tempfile()
  }
  con <- file("path_out", open = "w")
  for (i in seq_len(nrow(pairs))){
    writeLines(sprintf("SEQUENCE_ID=pair%d", pairs[i,"ID"]), con)
    writeLines("PRIMER_TASK=check_primers", con)
    writeLines(sprintf("SEQUENCE_PRIMER=%s", pairs[i,"sequence_F"]), con)
    writeLines(sprintf("SEQUENCE_PRIMER_REVCOMP=%s", pairs[i,"sequence_R"]), con)
    writeLines("=", con)
  }
  close(con)
  return(invisible(path_out))
}
