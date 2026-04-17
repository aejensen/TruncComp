.load_package_dataset <- function(name) {
  dataset_env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "TruncComp2", envir = dataset_env)

  if(!exists(name, envir = dataset_env, inherits = FALSE)) {
    stop("Could not load dataset ", name, ".", call. = FALSE)
  }

  get(name, envir = dataset_env, inherits = FALSE)
}

#' Deprecated compatibility wrapper for `trunc_comp_example`
#'
#' Returns the packaged example dataset formerly exposed through a dedicated
#' loader function. Use `data("trunc_comp_example", package = "TruncComp2")`
#' instead.
#'
#' @return A data frame identical to the `trunc_comp_example` dataset.
#' @examples
#' d <- load_trunc_comp2_example()
#' head(d)
#' @export
load_trunc_comp2_example <- function() {
  .Deprecated(msg = paste(
    "`load_trunc_comp2_example()` is deprecated;",
    "use `data(\"trunc_comp_example\", package = \"TruncComp2\")` instead."
  ))
  .load_package_dataset("trunc_comp_example")
}

#' Deprecated compatibility wrapper for `trunc_comp_adjusted_example`
#'
#' Returns the packaged adjusted example dataset formerly exposed through a
#' dedicated loader function. Use
#' `data("trunc_comp_adjusted_example", package = "TruncComp2")` instead.
#'
#' @return A data frame identical to the `trunc_comp_adjusted_example` dataset.
#' @examples
#' d <- load_trunc_comp2_adjusted_example()
#' head(d)
#' @export
load_trunc_comp2_adjusted_example <- function() {
  .Deprecated(msg = paste(
    "`load_trunc_comp2_adjusted_example()` is deprecated;",
    "use `data(\"trunc_comp_adjusted_example\", package = \"TruncComp2\")` instead."
  ))
  .load_package_dataset("trunc_comp_adjusted_example")
}

loadTruncComp2Example <- function() {
  .Deprecated("load_trunc_comp2_example")
  load_trunc_comp2_example()
}

loadTruncComp2AdjustedExample <- function() {
  .Deprecated("load_trunc_comp2_adjusted_example")
  load_trunc_comp2_adjusted_example()
}
