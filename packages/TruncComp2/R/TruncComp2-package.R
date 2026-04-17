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
#' `"TruncComp2"` with summaries, confidence intervals, simultaneous confidence
#' regions for unadjusted fits, and packaged example data helpers.
#'
#' @details
#' The package is organized around a standardized internal data representation
#' with columns `Y`, `A`, and `R`, optionally augmented with baseline covariates
#' supplied through `adjust`.
#'
#' Exported user-facing helpers include:
#'
#' - [truncComp()] and [truncComp.default()] for fitting the model.
#' - [summary.TruncComp2()] and [print.TruncComp2()] for displaying results.
#' - [confint.TruncComp2()] and [jointContrastCI()] for marginal and joint
#'   inference.
#' - [simulateTruncatedData()], [loadTruncComp2Example()], and
#'   [loadTruncComp2AdjustedExample()] for reproducible examples.
#'
#' For successful unadjusted fits the package stores the derived combined-outcome
#' contrast `Delta` as a point estimate. Confidence intervals for `Delta` are
#' computed on demand through [confint.TruncComp2()] using
#' `method = "welch"`, `"profile"`, or `"projected"`.
#'
#' Adjusted fits return conditional treatment effects from the fitted regression
#' components. In that setting, `Delta`, `Delta` intervals, and joint
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
