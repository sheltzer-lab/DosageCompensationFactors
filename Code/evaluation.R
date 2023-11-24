write_list <- function(report, report_file) {
  for (item in names(report)) {
    if (is.list(report[[item]])) {
      cat(sprintf("\n--- %s ---\n", item), file = report_file, append = TRUE)
      write_list(report[[item]], report_file)
      cat("------\n", file = report_file, append = TRUE)
    }
    else {
      cat(sprintf("%s: %s\n", item, report[[item]]), file = report_file, append = TRUE)
    }
  }
}

write_report <- function(report, report_file) {
  cat("===== Report =====\n", file = report_file)
  cat(format(Sys.time(), "%a %b %d %X %Y"), file = report_file, append = TRUE)
  cat("\n\n", file = report_file, append = TRUE)
  write_list(report, report_file)
}