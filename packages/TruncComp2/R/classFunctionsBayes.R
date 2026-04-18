bayes_contrast_labels <- function() {
  c(
    delta_atom = "Difference in atom probability",
    mu_delta = "Difference in means among non-atom outcomes",
    alpha_delta = "Odds ratio of being observed",
    delta = "Combined-outcome contrast (delta)"
  )
}

bayes_arm_labels <- function() {
  c(
    rho_0 = "Atom probability, control",
    rho_1 = "Atom probability, treatment",
    pi_0 = "Observed probability, control",
    pi_1 = "Observed probability, treatment",
    mu_0_c = "Mean among non-atom outcomes, control",
    mu_1_c = "Mean among non-atom outcomes, treatment"
  )
}

bayes_ppc_labels <- function() {
  c(
    atom = "Atom model",
    continuous = "Continuous model"
  )
}

bayes_interval_matrix <- function(draws, parameter, conf.level) {
  rows <- lapply(parameter, function(name) {
    bayes_equal_tail_interval(draws[[name]], conf.level)
  })

  mat <- do.call(rbind, rows)
  rownames(mat) <- parameter
  colnames(mat) <- ciColumnLabels(conf.level)
  mat
}

#' Summarize a Bayesian truncated-comparison fit
#'
#' Prints posterior treatment contrasts, arm-specific summaries, and sampler
#' diagnostics for a fitted `"trunc_comp_bayes_fit"` object.
#'
#' @param object A `"trunc_comp_bayes_fit"` object returned by
#'   [trunc_comp_bayes()].
#' @param ... Unused additional arguments.
#' @return Invisibly returns a list of class `"trunc_comp_bayes_summary"` with
#'   elements `contrasts`, `arms`, `ppc`, `diagnostics`, and `settings`.
#' @export
summary.trunc_comp_bayes_fit <- function(object, ...) {
  if(!is.null(object$call)) {
    cat("Call:\n")
    print(object$call)
    cat("\n")
  }

  cat("Method:", bayes_trunc_comp_method_label(), "\n")
  cat("Experimental: TRUE\n")
  cat("Confidence level = ", object$conf.level * 100, "%\n", sep = "")
  if(!is.null(object$atom)) {
    cat("Atom =", object$atom, "\n")
  }
  if(!is.null(object$settings$continuous_support)) {
    cat("Continuous support =", object$settings$continuous_support, "\n")
  }
  cat("\n")

  ppc_table <- NULL
  if(isTRUE(object$success)) {
    ppc_table <- tryCatch(
      posterior_predictive_pvalues(object),
      error = identity
    )

    if(inherits(ppc_table, "error")) {
      ppc_table <- NULL
    }
  }

  summary_object <- structure(
    list(
      contrasts = object$summary_table,
      arms = object$arm_table,
      ppc = ppc_table,
      diagnostics = object$diagnostics,
      settings = object$settings,
      success = object$success,
      error = object$error,
      call = object$call,
      conf.level = object$conf.level
    ),
    class = c("trunc_comp_bayes_summary", "list")
  )

  if(!isTRUE(object$success)) {
    cat("The Bayesian estimation procedure failed.\n")
    cat("The error message was:", object$error, "\n\n")
    return(invisible(summary_object))
  }

  contrast_labels <- bayes_contrast_labels()
  contrast_mat <- as.matrix(object$summary_table[, c("estimate", "conf.low", "conf.high", "posterior_prob")])
  rownames(contrast_mat) <- unname(contrast_labels[rownames(object$summary_table)])
  colnames(contrast_mat) <- c("Estimate", "CrI Lower", "CrI Upper", "Posterior Prob")

  arm_labels <- bayes_arm_labels()
  arm_mat <- as.matrix(object$arm_table[, c("estimate", "conf.low", "conf.high")])
  rownames(arm_mat) <- unname(arm_labels[rownames(object$arm_table)])
  colnames(arm_mat) <- c("Estimate", "CrI Lower", "CrI Upper")

  cat("Treatment contrasts\n")
  print.default(contrast_mat)
  cat("\nPosterior probabilities are against 0, except alpha_delta which is against 1.\n\n")

  cat("Arm-specific summaries\n")
  print.default(arm_mat)
  cat("\n")

  if(!is.null(ppc_table)) {
    ppc_labels <- bayes_ppc_labels()
    ppc_print <- ppc_table[, c("p_value", "statistic", "scale", "ndraws"), drop = FALSE]
    rownames(ppc_print) <- unname(ppc_labels[rownames(ppc_table)])
    colnames(ppc_print) <- c("Posterior Predictive P", "Statistic", "Scale", "Draws")

    cat("Posterior predictive checks\n")
    print(ppc_print)
    cat("\n")
  }

  diagnostics <- object$diagnostics
  cat("Sampler diagnostics\n")
  cat("  Divergences:", diagnostics$divergences, "\n")
  cat("  Max Rhat:", format(diagnostics$max_rhat, digits = 4), "\n")
  cat("  Min bulk ESS:", format(diagnostics$min_bulk_ess, digits = 6), "\n")
  cat("  Min tail ESS:", format(diagnostics$min_tail_ess, digits = 6), "\n")
  cat("  Diagnostic OK:", isTRUE(diagnostics$diagnostic_ok), "\n\n")

  if(!isTRUE(diagnostics$diagnostic_ok)) {
    cat("Sampler diagnostics indicate potential issues. Interpret posterior summaries with caution.\n\n")
  }

  invisible(summary_object)
}

#' Print a Bayesian truncated-comparison fit
#'
#' Prints a concise overview of a fitted `"trunc_comp_bayes_fit"` object.
#'
#' @param x A `"trunc_comp_bayes_fit"` object returned by [trunc_comp_bayes()].
#' @param ... Unused additional arguments.
#' @return The input object, returned invisibly.
#' @export
print.trunc_comp_bayes_fit <- function(x, ...) {
  cat("<trunc_comp_bayes_fit>\n")

  if(!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
  }

  cat("Method:", bayes_trunc_comp_method_label(), "\n")
  cat("Experimental:", TRUE, "\n")
  if(!is.null(x$settings$continuous_support)) {
    cat("Continuous support:", x$settings$continuous_support, "\n")
  }
  cat("Success:", isTRUE(x$success), "\n")

  if(isTRUE(x$success)) {
    estimates <- c(
      delta = x$summary_table["delta", "estimate"],
      mu_delta = x$summary_table["mu_delta", "estimate"],
      delta_atom = x$summary_table["delta_atom", "estimate"],
      alpha_delta = x$summary_table["alpha_delta", "estimate"]
    )
    print.default(estimates)
    cat("\nUse summary() for posterior intervals and diagnostics.\n")
  } else {
    cat("Error:", x$error, "\n")
  }

  invisible(x)
}

#' Extract posterior point summaries from a Bayesian truncated-comparison fit
#'
#' @param object A `"trunc_comp_bayes_fit"` object returned by
#'   [trunc_comp_bayes()].
#' @param ... Unused additional arguments.
#' @return A named numeric vector containing posterior medians for `mu_delta`,
#'   `delta_atom`, `alpha_delta`, and `delta`, where `delta_atom` is the
#'   treatment-minus-control difference in atom probability, `mu_delta` is the
#'   treatment-minus-control difference in mean among non-atom outcomes,
#'   `alpha_delta` is the treatment-to-control odds ratio for being observed
#'   away from the atom, and `delta` is the combined-outcome contrast on the
#'   original outcome scale.
#' @export
coef.trunc_comp_bayes_fit <- function(object, ...) {
  if(!isTRUE(object$success)) {
    return(stats::setNames(rep(NA_real_, 4), c("mu_delta", "delta_atom", "alpha_delta", "delta")))
  }

  stats::setNames(
    c(
      object$summary_table["mu_delta", "estimate"],
      object$summary_table["delta_atom", "estimate"],
      object$summary_table["alpha_delta", "estimate"],
      object$summary_table["delta", "estimate"]
    ),
    c("mu_delta", "delta_atom", "alpha_delta", "delta")
  )
}

#' Credible intervals for a Bayesian truncated-comparison fit
#'
#' Computes equal-tail credible intervals from the stored posterior draws in a
#' fitted `"trunc_comp_bayes_fit"` object.
#'
#' @param object A `"trunc_comp_bayes_fit"` object returned by
#'   [trunc_comp_bayes()].
#' @param parameter Parameter selection. Supported values are `"mu_delta"`,
#'   `"delta_atom"`, `"alpha_delta"`, and `"delta"`, with the same definitions
#'   as in [trunc_comp_bayes()].
#' @param conf.level Credible level for the requested interval.
#' @param method Unsupported for Bayesian fits and must be left missing.
#' @param ... Unsupported additional arguments.
#' @return Invisibly returns a printed matrix of equal-tail credible intervals.
#' @export
confint.trunc_comp_bayes_fit <- function(object,
                                         parameter = bayes_parameter_names("contrast"),
                                         conf.level = object$conf.level,
                                         method,
                                         ...) {
  dots <- list(...)

  if(!missing(method) && !is.null(method)) {
    stop("Bayesian credible intervals do not support method-specific variants.")
  }

  if(length(dots) > 0L) {
    stop("Additional arguments are not supported for Bayesian credible intervals.")
  }

  if(!isTRUE(object$success)) {
    stop("Estimation failed. Cannot display credible intervals.")
  }

  conf.level <- validateConfidenceLevel(conf.level)
  parameter <- bayes_parameter_aliases(parameter)
  parameter <- unique(match.arg(
    parameter,
    choices = bayes_parameter_names("contrast"),
    several.ok = TRUE
  ))

  interval_mat <- bayes_interval_matrix(object$draws, parameter, conf.level)
  print.default(interval_mat)
  invisible(interval_mat)
}
