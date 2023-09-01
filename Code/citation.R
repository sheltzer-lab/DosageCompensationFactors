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