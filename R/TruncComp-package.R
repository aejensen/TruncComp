#' TruncComp: Two-component comparisons for outcomes with a substantive atom
#'
#' Implements parametric and semiparametric likelihood ratio procedures for
#' comparing two groups when an endpoint comprises a substantive atom and a
#' continuous non-atom outcome.
#'
#' @description
#' `TruncComp` provides a common user interface for three analysis paths:
#'
#' - `method = "lrt"` fits the parametric likelihood ratio procedure.
#' - `method = "splrt"` fits the semiparametric likelihood ratio procedure with
#'   empirical likelihood for the conditional non-atom mean.
#' - [trunc_comp_bayes()] fits the experimental Bayesian two-part Dirichlet
#'   process mixture model for the no-covariate setting. Supported non-atom
#'   kernels cover real-valued, positive, bounded continuous, and bounded
#'   reported score outcomes.
#'
#' The main entry point is [trunc_comp()], which returns an object of class
#' `"trunc_comp_fit"` with component summaries, confidence intervals, and
#' simultaneous confidence regions. The package also supplies reproducible
#' example datasets.
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
#' - [posterior_predictive_pvalues()] for discrepancy-based posterior
#'   predictive p values from the Bayesian fit.
#' - [summary.trunc_comp_fit()] and [print.trunc_comp_fit()] for displaying results.
#' - [summary.trunc_comp_bayes_fit()] and [print.trunc_comp_bayes_fit()] for
#'   Bayesian posterior summaries.
#' - [confint.trunc_comp_fit()] and [joint_contrast_surface()] for marginal and
#'   joint inference.
#' - [simulate_truncated_data()], `trunc_comp_example`, and
#'   `trunc_comp_adjusted_example` for reproducible examples.
#'
#' For successful unadjusted fits the package stores the coded endpoint mean
#' difference `delta`. Confidence intervals for `delta` are
#' computed on demand through [confint.trunc_comp_fit()] using
#' `method = "welch"`, `"profile"`, or `"projected"`.
#'
#' Adjusted fits return conditional treatment contrasts from the fitted
#' regression components. Their joint component confidence region profiles the
#' adjustment coefficients. Because those coefficients do not define a marginal
#' coded endpoint contrast without a standardization target, `delta` and its
#' intervals are unavailable for adjusted fits.
#'
#' @docType package
#' @name TruncComp-package
#' @aliases TruncComp
#' @keywords package
#' @useDynLib TruncComp, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom stats get_all_vars model.matrix optimize terms
"_PACKAGE"
