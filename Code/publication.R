# Export Citations
## ToDo: Write to file

getPackagesInCurrentSession <- function() {
  session <- sessionInfo()
  return(c('base', names(session$otherPkgs)))
}

createCitationList <- function(packages) {
  packages <- unique(packages)
  citations <- unlist(lapply(X = packages, FUN = citation), recursive = FALSE)
  return(citations)
}

# Create README for supplemental tables exported in the Apache Parquet format
createParquetReadme <- function(table_description, col_descriptions, title = "Supplemental Data - README",
                                readme_path = "README.md", file_path = "filename.parquet", parquet_version = "2.4") {
  header <- c(
    paste("#", title),
    "===========================",
    "",
    paste0("The supplemental data \"", basename(file_path), "\" accompanied by this README is stored in the Apache Parquet format, version ", parquet_version, " (https://parquet.apache.org/)."),
    "Parquet is a column-oriented data file format supporting data type annotations and compression.",
    "",
    "You can open Parquet files in R using the 'arrow' package (https://arrow.apache.org/docs/r/index.html).",
    "With the arrow package installed you can open the supplemental data in R by executing:",
    "",
    "```r",
    "library(arrow)",
    paste0("df <- read_parquet(\"", basename(file_path), "\")"),
    "```",
    "",
    "## Table Description",
    "--------------------",
    table_description,
    "",
    "## Column Descriptions",
    "--------------------"
  )

  descriptions <- paste0("• ", names(col_descriptions), ": ", col_descriptions)
  content <- c(header, descriptions)

  if (file.exists(readme_path)) file.remove(readme_path)
  writeLines(content, con = readme_path, sep = "\n")
}