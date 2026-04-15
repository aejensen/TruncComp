load_local_trunccomp <- function(repo_root) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The pkgload package is required to load the local TruncComp package.", call. = FALSE)
  }

  pkgload::load_all(
    file.path(repo_root, "packages", "TruncComp"),
    quiet = TRUE,
    export_all = FALSE,
    helpers = FALSE,
    attach_testthat = FALSE
  )

  invisible(TRUE)
}

load_example_data <- function(repo_root) {
  data_env <- new.env(parent = emptyenv())
  load(file.path(repo_root, "packages", "TruncComp", "data", "TruncCompExample.RData"), envir = data_env)
  data_env$TruncCompExample
}

load_application_data <- function(manuscript_dir) {
  utils::read.csv(file.path(manuscript_dir, "covidData.csv"), sep = ";")
}

load_simulation_results <- function(manuscript_dir) {
  data_env <- new.env(parent = emptyenv())
  load(file.path(manuscript_dir, "simulations", "simulationResults.RData"), envir = data_env)
  as.list(data_env)
}
