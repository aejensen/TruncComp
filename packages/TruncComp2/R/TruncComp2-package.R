#' TruncComp2: Two-sample comparison of truncated continuous outcomes
#'
#' Implements parametric and semi-parametric likelihood-ratio procedures for
#' comparing two groups when a continuous outcome has a distinguished atom value
#' that denotes an unobserved or undefined outcome.
#'
#' @description
#' `TruncComp2` provides a common user interface for two estimation paths:
#'
#' - `method = "lrt"` fits the fully parametric likelihood-ratio procedure.
#' - `method = "splrt"` fits the semi-parametric likelihood-ratio procedure with
#'   an internal empirical-likelihood engine for the observed-outcome component.
#' - [trunc_comp_bayes()] fits the experimental Bayesian two-part Dirichlet
#'   process mixture model for the no-covariate setting, with either real-line
#'   Gaussian kernels or positive-support Gamma kernels for the non-atom
#'   outcomes.
#'
#' The main entry point is [trunc_comp()], which returns an object of class
#' `"trunc_comp_fit"` with summaries, confidence intervals, simultaneous confidence
#' regions for unadjusted fits, and packaged example datasets.
#'
#' @details
#' The package is organized around a standardized internal data representation
#' with columns `Y`, `A`, and `R`, optionally augmented with baseline covariates
#' supplied through `adjust`.
#'
#' Exported user-facing helpers include:
#'
#' - [trunc_comp()] for fitting the model.
#' - [trunc_comp_bayes()] for the experimental Bayesian pathway.
#' - [posterior_density_plot()] for plotting the arm-specific posterior outcome
#'   densities implied by the Bayesian fit.
#' - [posterior_predictive_check()] for `bayesplot`-based posterior predictive
#'   checks of the Bayesian fit.
#' - [summary.trunc_comp_fit()] and [print.trunc_comp_fit()] for displaying results.
#' - [summary.trunc_comp_bayes_fit()] and [print.trunc_comp_bayes_fit()] for
#'   Bayesian posterior summaries.
#' - [confint.trunc_comp_fit()] and [joint_contrast_surface()] for marginal and joint
#'   inference.
#' - [simulate_truncated_data()], `trunc_comp_example`, and
#'   `trunc_comp_adjusted_example` for reproducible examples.
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
#' @useDynLib TruncComp2, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom stats get_all_vars model.matrix optimize terms
"_PACKAGE"
