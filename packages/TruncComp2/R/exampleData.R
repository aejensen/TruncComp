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

#' Load the packaged TruncComp2 example data
#'
#' Loads the fixed example data used in the package documentation, tests, and
#' examples.
#'
#' @return A data frame with 50 observations and 2 variables: `R`, the binary
#'   treatment indicator coded `0` or `1`, and `Y`, the combined outcome where
#'   the atom value `0` denotes an unobserved or undefined continuous outcome.
#' @details
#' The example is stored as an `.rds` file under
#' `inst/extdata/TruncComp2Example.rds`. This helper provides a stable public
#' loading mechanism without relying on `data()`.
#' @examples
#' d <- loadTruncComp2Example()
#' head(d)
#' @export
loadTruncComp2Example <- function() {
  readRDS(.truncComp2ExtdataPath("TruncComp2Example.rds"))
}

#' Load the packaged adjusted TruncComp2 example data
#'
#' Loads the fixed example data used to illustrate additive covariate adjustment
#' in both the parametric and semi-parametric procedures.
#'
#' @return A data frame with 50 observations and 3 variables: `R`, the binary
#'   treatment indicator; `L`, a three-level categorical baseline covariate with
#'   levels `low`, `mid`, and `high`; and `Y`, the combined outcome where the
#'   atom value `0` denotes an unobserved or undefined continuous outcome.
#' @details
#' The example is stored as an `.rds` file under
#' `inst/extdata/TruncComp2AdjustedExample.rds`. This helper provides a stable
#' public loading mechanism without relying on `data()`.
#' @examples
#' d <- loadTruncComp2AdjustedExample()
#' head(d)
#' @export
loadTruncComp2AdjustedExample <- function() {
  readRDS(.truncComp2ExtdataPath("TruncComp2AdjustedExample.rds"))
}
