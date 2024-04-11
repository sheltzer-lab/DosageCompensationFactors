report_list2lines <- function(report) {
  lines <- NULL
  for (item in names(report)) {
    if (is.list(report[[item]])) {
      lines <- c(lines,
                 sprintf("\n--- %s ---\n", item),
                 report_list2lines(report[[item]]),
                 "------\n")
    }
    else {
      lines <- c(lines, sprintf("%s: %s\n", item, report[[item]]))
    }
  }
  return(lines)
}

write_report <- function(report, report_file) {
  file.remove(report_file)
  lines <- c("===== Report =====\n",
             format(Sys.time(), "%a %b %d %X %Y"),
             "\n\n",
             report_list2lines(report))
  writeLines(lines, con = report_file, sep = "")
}

df2reportlist <- function(df, collapse_cols = "\t\t") {
  result <- df %>%
    as.list() %>%
    transpose() %>%
    lapply(\(lst) paste(lst, collapse = collapse_cols))
  result <- append(result, paste(colnames(df), collapse = collapse_cols), after = 0)
  names(result) <- c(0, rownames(df))   # Using "Columns" as label of first row may cause row to shift
  return(result)
}