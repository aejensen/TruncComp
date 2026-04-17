#' TruncComp2: Two-sample comparison of truncated continuous outcomes
#'
#' Implements parametric and semi-parametric likelihood-ratio procedures for
#' comparing two groups when a continuous outcome has a distinguished atom value
#' that denotes an unobserved or undefined outcome.
#'
#' @description
#' `TruncComp2` provides a common user interface for two estimation paths:
#'
#' - `method = "LRT"` fits the fully parametric likelihood-ratio procedure.
#' - `method = "SPLRT"` fits the semi-parametric likelihood-ratio procedure with
#'   an internal empirical-likelihood engine for the observed-outcome component.
#'
#' The main entry point is [truncComp()], which returns an object of class
#' `"trunc_comp_fit"` with summaries, confidence intervals, simultaneous confidence
#' regions for unadjusted fits, and packaged example data helpers.
#'
#' @details
#' The package is organized around a standardized internal data representation
#' with columns `Y`, `A`, and `R`, optionally augmented with baseline covariates
#' supplied through `adjust`.
#'
#' Exported user-facing helpers include:
#'
#' - [truncComp()] for fitting the model.
#' - [summary.trunc_comp_fit()] and [print.trunc_comp_fit()] for displaying results.
#' - [confint.trunc_comp_fit()] and [joint_contrast_ci()] for marginal and joint
#'   inference.
#' - [simulate_truncated_data()], [load_trunc_comp2_example()], and
#'   [load_trunc_comp2_adjusted_example()] for reproducible examples.
#'
#' For successful unadjusted fits the package stores the derived combined-outcome
#' contrast `delta` as a point estimate. Confidence intervals for `delta` are
#' computed on demand through [confint.trunc_comp_fit()] using
#' `method = "welch"`, `"profile"`, or `"projected"`.
#'
#' Adjusted fits return conditional treatment effects from the fitted regression
#' components. In that setting, `delta`, `delta` intervals, and joint
#' confidence regions are intentionally unavailable.
#'
#' @references
#' Andreas Kryger Jensen and Theis Lange.
#' *A novel high-power test for continuous outcomes truncated by death*.
#'
#' @docType package
#' @name TruncComp2-package
#' @aliases TruncComp2
#' @keywords package
#' @importFrom stats get_all_vars model.matrix optimize terms
"_PACKAGE"
