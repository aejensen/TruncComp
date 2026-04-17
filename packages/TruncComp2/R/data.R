#' Example Truncated-Outcome Data
#'
#' Fixed unadjusted example data used throughout the package documentation,
#' regression tests, and worked examples.
#'
#' @format A data frame with 50 rows and 2 variables:
#' \describe{
#'   \item{R}{Binary treatment indicator coded `0` or `1`.}
#'   \item{Y}{Combined outcome where the atom value `0` denotes an unobserved or
#'   undefined continuous outcome.}
#' }
#' @usage data("trunc_comp_example")
"trunc_comp_example"

#' Adjusted Example Truncated-Outcome Data
#'
#' Fixed adjusted example data used to illustrate additive covariate adjustment
#' in both the parametric and semi-parametric procedures.
#'
#' @format A data frame with 50 rows and 3 variables:
#' \describe{
#'   \item{R}{Binary treatment indicator coded `0` or `1`.}
#'   \item{L}{Three-level baseline covariate with levels `low`, `mid`, and `high`.}
#'   \item{Y}{Combined outcome where the atom value `0` denotes an unobserved or
#'   undefined continuous outcome.}
#' }
#' @usage data("trunc_comp_adjusted_example")
"trunc_comp_adjusted_example"
