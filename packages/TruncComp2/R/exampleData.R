.truncComp2ExtdataPath <- function(filename) {
  ns_path <- getNamespaceInfo(asNamespace("TruncComp2"), "path")
  candidates <- c(
    file.path(ns_path, "extdata", filename),
    file.path(ns_path, "inst", "extdata", filename)
  )

  path <- candidates[file.exists(candidates)][1]
  if(length(path) == 0 || !nzchar(path)) {
    stop("Could not locate ", filename, ".")
  }

  path
}

loadTruncComp2Example <- function() {
  readRDS(.truncComp2ExtdataPath("TruncComp2Example.rds"))
}

loadTruncComp2AdjustedExample <- function() {
  readRDS(.truncComp2ExtdataPath("TruncComp2AdjustedExample.rds"))
}
