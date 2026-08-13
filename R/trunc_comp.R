normalizeFormulaAdjust <- function(adjust, outcome_name, treatment_name) {
  if(is.null(adjust)) {
    return(NULL)
  }

  if(!inherits(adjust, "formula") || length(adjust) != 2) {
    stop("adjust must be NULL or a one-sided formula like ~ age + sex.")
  }

  adjust_terms <- stats::terms(adjust)
  labels <- attr(adjust_terms, "term.labels")
  orders <- attr(adjust_terms, "order")
  if("." %in% labels) {
    stop("adjust must not use '.'. Please list adjustment covariates explicitly.")
  }

  if(any(orders > 1)) {
    stop("adjust must be additive and must not contain interaction terms.")
  }

  adjust_vars <- all.vars(adjust)
  if(outcome_name %in% adjust_vars || treatment_name %in% adjust_vars) {
    stop("adjust must not include the outcome or treatment variable from the main formula.")
  }

  reserved_names <- c("Y", "A", "R")
  if(any(adjust_vars %in% reserved_names)) {
    stop("adjust must not use the reserved analysis variable names Y, A, or R.")
  }

  if(length(labels) == 0) {
    return(NULL)
  }

  adjust
}

prepareFormulaAdjustment <- function(data, adjust) {
  if(is.null(adjust)) {
    return(list(data = NULL, formula = NULL))
  }

  adjust_vars <- all.vars(adjust)
  if(length(adjust_vars) == 0) {
    return(list(data = NULL, formula = NULL))
  }

  adjust_data <- stats::model.frame(
    stats::reformulate(adjust_vars),
    data = data,
    na.action = stats::na.fail
  )

  if(anyDuplicated(names(adjust_data))) {
    stop("Adjustment covariate names must be unique.")
  }

  list(data = adjust_data, formula = adjust)
}

prepareDefaultAdjustment <- function(adjust, n) {
  if(is.null(adjust)) {
    return(list(data = NULL, formula = NULL))
  }

  if(is.matrix(adjust)) {
    adjust <- as.data.frame(adjust, stringsAsFactors = FALSE)
    if(is.null(colnames(adjust))) {
      colnames(adjust) <- paste0("L", seq_len(ncol(adjust)))
    }
  }

  if(!is.data.frame(adjust)) {
    stop("adjust must be NULL, a data.frame, or a matrix.")
  }

  if(nrow(adjust) != n) {
    stop("adjust must have the same number of rows as y.")
  }

  if(ncol(adjust) == 0) {
    return(list(data = NULL, formula = NULL))
  }

  if(is.null(names(adjust)) || any(names(adjust) == "")) {
    stop("All adjustment covariates must have non-empty names.")
  }

  if(anyDuplicated(names(adjust))) {
    stop("Adjustment covariate names must be unique.")
  }

  if(any(names(adjust) %in% c("Y", "A", "R"))) {
    stop("adjust must not use the reserved analysis variable names Y, A, or R.")
  }

  adjust_data <- stats::na.fail(adjust)
  adjust_formula <- stats::reformulate(names(adjust_data))

  list(data = adjust_data, formula = adjust_formula)
}

resolveDefaultAtom <- function(y, a, atom = NULL) {
  if(!is.null(atom)) {
    if(!(length(atom) == 1 && is.numeric(atom) && is.finite(atom))) {
      stop("atom must be NULL or a single finite numeric value.")
    }
    return(as.numeric(atom))
  }

  atom_values <- unique(y[a == 0])
  atom_values <- atom_values[is.finite(atom_values)]
  if(length(atom_values) == 1) {
    return(as.numeric(atom_values))
  }

  stop("atom must be supplied for trunc_comp(y, a, r, ...) unless y[a == 0] has exactly one unique finite value.")
}

trunc_comp_core <- function(y, a, r, method, conf.level = 0.95,
                            adjust_data = NULL, adjust_formula = NULL,
                            atom = NULL, call = NULL) {
  method <- normalize_trunc_comp_method(method)
  a <- as.numeric(a)
  r <- as.numeric(r)

  d <- data.frame(Y = y, A = a, R = r)
  if(!is.null(adjust_data)) {
    d <- cbind(d, adjust_data)
  }

  adjust_spec <- adjustmentSpecification(adjust_formula)

  if(!isDataOkay(d)) {
    error <- "Estimation failed due to data error."
    warning(error)
    return(new_failed_trunc_comp_fit(
      error,
      method,
      conf.level,
      data = d,
      adjust = adjust_spec,
      adjust_formula = adjust_formula,
      atom = atom,
      call = call
    ))
  }

  if(method == "lrt") {
    out <- LRT(
      data = d,
      conf.level = conf.level,
      adjust = adjust_formula,
      adjust_spec = adjust_spec,
      atom = atom,
      call = call
    )
  } else {
    out <- SPLRT(
      data = d,
      conf.level = conf.level,
      adjust = adjust_formula,
      adjust_spec = adjust_spec,
      atom = atom,
      call = call
    )
  }

  out$data <- d
  out$atom <- atom
  out$call <- call
  out
}

#' Fit a two-component treatment comparison
#'
#' Compares two groups when an endpoint comprises a substantive atom and a
#' continuous non-atom outcome.
#'
#' @param formula For the formula interface, a formula with the coded endpoint
#'   on the left-hand side and a single binary treatment indicator on the
#'   right-hand side. For the default interface, the first positional argument
#'   is the outcome storage vector.
#' @param atom A single numeric value used to represent the atom. The non-atom
#'   indicator is reconstructed internally as `Y != atom`. For the
#'   default interface, `atom` may be omitted only when `y[a == 0]` has exactly
#'   one unique finite value, in which case it is inferred and stored on the
#'   fitted object.
#' @param data A data frame containing the variables referenced by `formula`.
#' @param method Either `"lrt"` for the parametric likelihood ratio method or
#'   `"splrt"` for the semiparametric likelihood ratio method. Values are
#'   matched without regard to case.
#' @param conf.level Confidence level used for the reported intervals.
#' @param adjust Optional covariate adjustment. For the formula interface this
#'   must be `NULL` or a one-sided additive formula such as `~ age + sex`. For
#'   the default interface this must be `NULL`, a data frame, or a matrix of
#'   baseline covariates. The same additive adjustment is used in both the
#'   conditional non-atom mean and Bernoulli components. Adjusted fits support
#'   component confidence intervals and a joint component confidence surface.
#' @param ... Additional arguments are not supported.
#' @return An S3 object of class `"trunc_comp_fit"` with component estimates,
#'   confidence intervals, the joint likelihood ratio statistic, metadata about
#'   the fitting method, the standardized analysis data, the matched call, and
#'   the fitted atom value. Failed fits return the same class with
#'   `success = FALSE` and an error message.
#' @details
#' For successful unadjusted fits, the package also computes the coded endpoint
#' mean difference
#'
#' `delta = [p1 * mu1 + (1 - p1) * atom] - [p0 * mu0 + (1 - p0) * atom]`.
#'
#' The fitted object stores only the `delta` point estimate. Confidence
#' intervals for `delta` are computed on demand through
#' `confint(fit, parameter = "delta", method = "welch" | "profile" | "projected")`.
#'
#' Adjusted fits return `NA` for `delta` because their treatment contrasts are
#' conditional regression coefficients rather
#' than standardized marginal contrasts. Their joint component confidence
#' surface remains available through `confint(fit, parameter = "joint")`.
#'
#' @seealso [summary.trunc_comp_fit()], [print.trunc_comp_fit()],
#'   [confint.trunc_comp_fit()], [joint_contrast_surface()]
#' @examples
#' library(TruncComp)
#' f0 <- function(n) stats::rnorm(n, 3, 1)
#' f1 <- function(n) stats::rnorm(n, 3.5, 1)
#' d <- simulate_truncated_data(n = 100, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)
#'
#' # Formula interface
#' trunc_comp(Y ~ R, atom = 0, data = d, method = "lrt")
#' trunc_comp(Y ~ R, atom = 0, data = d, method = "splrt")
#'
#' data("trunc_comp_adjusted_example", package = "TruncComp")
#' trunc_comp(Y ~ R, atom = 0, data = trunc_comp_adjusted_example, method = "lrt", adjust = ~ L)
#' trunc_comp(Y ~ R, atom = 0, data = trunc_comp_adjusted_example, method = "splrt", adjust = ~ L)
#'
#' # Default interface
#' trunc_comp(d$Y, d$A, d$R, method = "lrt", atom = 0)
#' @rdname trunc_comp
#' @export
trunc_comp <- function(formula, ...) {
  UseMethod("trunc_comp")
}

#' @rdname trunc_comp
#' @export
trunc_comp.formula <- function(formula, atom, data, method = c("lrt", "splrt"),
                               conf.level = 0.95, adjust = NULL, ...) {
  if(length(list(...)) > 0) {
    stop("Additional arguments are not supported.")
  }

  if(!inherits(formula, "formula")) {
    stop("The formula must be a formula.")
  }

  method <- normalize_trunc_comp_method(method)

  if(length(attr(terms(formula), "term.labels")) != 1) {
    stop("The formula must contain exactly one treatment term.")
  }

  outcome_name <- all.vars(formula[[2]])
  treatment_name <- all.vars(formula[[3]])
  if(length(outcome_name) != 1 || length(treatment_name) != 1) {
    stop("The formula must have a single outcome and a single binary treatment variable.")
  }

  outcome_name <- outcome_name[[1]]
  treatment_name <- treatment_name[[1]]
  adjust <- normalizeFormulaAdjust(adjust, outcome_name, treatment_name)

  variables <- stats::model.frame(formula, data = data, na.action = stats::na.fail)
  outcome <- variables[,1]
  treatment <- variables[,2]
  adjustment <- prepareFormulaAdjustment(data, adjust)

  if(!(is.numeric(outcome) && all(is.finite(outcome)))) {
    stop("The outcome must contain only finite numeric values.")
  }

  treatment_levels <- sort(unique(treatment))
  if(!(length(treatment_levels) == 2L &&
       all(treatment_levels == c(0, 1)))) {
    stop("The treatment indicator must be binary and coded 0/1.")
  }

  if(!(length(atom) == 1L && is.numeric(atom) && is.finite(atom))) {
    stop("atom must be a single finite numeric value.")
  }

  non_atom <- as.numeric(outcome != atom)

  if(all(non_atom[treatment == 0] == 0) | all(non_atom[treatment == 1] == 0)) {
    stop("One treatment group has no non-atom outcomes, so estimation is unavailable.")
  }

  if(all(non_atom == 1)) {
    stop("The data contain no atom observations, so a two-component analysis is not applicable.")
  }


  trunc_comp_core(
    outcome,
    non_atom,
    treatment,
    method,
    conf.level,
    adjust_data = adjustment$data,
    adjust_formula = adjustment$formula,
    atom = atom,
    call = match.call()
  )
}

#' Default interface for [trunc_comp()]
#'
#' @param a Binary indicator for whether the endpoint is outside the atom.
#' @param r Binary treatment indicator for the default interface.
#' @rdname trunc_comp
#' @exportS3Method trunc_comp default
trunc_comp.default <- function(formula, a, r, method = c("lrt", "splrt"),
                               conf.level = 0.95, adjust = NULL, atom = NULL, ...) {
  if(length(list(...)) > 0) {
    stop("Additional arguments are not supported.")
  }

  y <- formula
  method <- normalize_trunc_comp_method(method)

  if(!(is.numeric(y) && length(y) > 0L && all(is.finite(y)))) {
    stop("y must contain only finite numeric values.")
  }
  if(length(a) != length(y) || length(r) != length(y)) {
    stop("y, a, and r must have the same length.")
  }
  if(!(is.numeric(a) || is.logical(a)) ||
     any(is.na(a)) || !all(a %in% c(0, 1))) {
    stop("a must be binary and coded 0/1.")
  }
  if(!(is.numeric(r) || is.logical(r)) ||
     any(is.na(r)) || !all(r %in% c(0, 1)) ||
     !setequal(unique(r), c(0, 1))) {
    stop("r must contain both treatment groups and be coded 0/1.")
  }

  adjustment <- prepareDefaultAdjustment(adjust, length(y))
  atom <- resolveDefaultAtom(y, a, atom = atom)

  if(all(a == 1)) {
    stop("The data contain no atom observations, so a two-component analysis is not applicable.")
  }

  trunc_comp_core(
    y,
    a,
    r,
    method,
    conf.level,
    adjust_data = adjustment$data,
    adjust_formula = adjustment$formula,
    atom = atom,
    call = match.call()
  )
}
