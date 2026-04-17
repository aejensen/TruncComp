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

  stop("atom must be supplied for truncComp(y, a, r, ...) unless y[a == 0] has exactly one unique finite value.")
}

truncComp_core <- function(y, a, r, method, conf.level = 0.95, init = NULL,
                           adjust_data = NULL, adjust_formula = NULL,
                           atom = NULL) {
  if(!(method == "LRT" || method == "SPLRT")) {
    stop("Only LRT or SPLRT supported as methods.")
  }

  d <- data.frame(Y = y, A = a, R = r)
  if(!is.null(adjust_data)) {
    d <- cbind(d, adjust_data)
  }

  adjust_spec <- adjustmentSpecification(adjust_formula)

  if(!isDataOkay(d)) {
    error <- "Estimation failed due to data error."
    warning(error)
    return(returnErrorData(error,
                           method,
                           conf.level,
                           init = init,
                           data = d,
                           adjust = adjust_spec,
                           atom = atom))
  }

  if(method == "LRT") {
    out <- LRT(d, init, conf.level, adjust = adjust_formula, adjust_spec = adjust_spec, atom = atom)
  } else if(method == "SPLRT") {
    out <- SPLRT(d, conf.level, adjust = adjust_formula, adjust_spec = adjust_spec, atom = atom)
  }

  out$data <- d
  out$atom <- atom
  out
}

#' Fit the TruncComp2 two-sample comparison model
#'
#' Compares two groups when a continuous outcome has a distinguished atom value
#' representing an unobserved or undefined outcome.
#'
#' @param formula A formula with a continuous outcome on the left-hand side and a
#'   single binary treatment indicator on the right-hand side.
#' @param atom A single numeric value used for the special atom outcome. The
#'   observation indicator is reconstructed internally as `Y != atom`. For the
#'   default interface, `atom` may be omitted only when `y[a == 0]` has exactly
#'   one unique finite value, in which case it is inferred and stored on the
#'   fitted object.
#' @param data A data frame containing the variables referenced by `formula`.
#' @param method Either `"LRT"` for the parametric likelihood-ratio method or
#'   `"SPLRT"` for the semi-parametric likelihood-ratio method.
#' @param conf_level Confidence level used for the reported intervals.
#' @param init Optional compatibility argument retained from older parametric
#'   implementations. The current model-backed `LRT` path stores it on the
#'   returned object but does not use it for estimation.
#' @param adjust Optional covariate adjustment. For the formula interface this
#'   must be `NULL` or a one-sided additive formula such as `~ age + sex`. For
#'   the default interface this must be `NULL`, a data frame, or a matrix of
#'   baseline covariates. The same additive adjustment is used in both the
#'   observed-outcome and observation components. Adjusted fits support only
#'   component confidence intervals.
#' @param ... Unused additional arguments.
#' @return An S3 object of class `"trunc_comp_fit"` with component estimates,
#'   confidence intervals, the joint likelihood-ratio statistic, metadata about
#'   the fitting method, the standardized analysis data, and the fitted atom
#'   value. Failed fits return the same class with `success = FALSE` and an
#'   error message.
#' @details
#' For successful unadjusted fits, the package also computes a derived
#' combined-outcome contrast
#'
#' `delta = [p1 * mu1 + (1 - p1) * atom] - [p0 * mu0 + (1 - p0) * atom]`.
#'
#' The fitted object stores only the `delta` point estimate. Confidence
#' intervals for `delta` are computed on demand through
#' `confint(fit, parameter = "delta", method = "welch" | "profile" | "projected")`.
#'
#' Adjusted fits deliberately return `NA` for `delta` because the implemented
#' adjusted treatment effects are conditional regression coefficients rather
#' than standardized marginal contrasts.
#'
#' @seealso [summary.trunc_comp_fit()], [print.trunc_comp_fit()],
#'   [confint.trunc_comp_fit()], [joint_contrast_ci()]
#' @examples
#' library(TruncComp2)
#' f0 <- function(n) stats::rnorm(n, 3, 1)
#' f1 <- function(n) stats::rnorm(n, 3.5, 1)
#' d <- simulate_truncated_data(n = 100, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)
#'
#' # Formula interface
#' truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
#' truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
#'
#' d_adjusted <- load_trunc_comp2_adjusted_example()
#' truncComp(Y ~ R, atom = 0, data = d_adjusted, method = "LRT", adjust = ~ L)
#' truncComp(Y ~ R, atom = 0, data = d_adjusted, method = "SPLRT", adjust = ~ L)
#'
#' # Default interface
#' truncComp(d$Y, d$A, d$R, method = "LRT", atom = 0)
#' @rdname truncComp
#' @export
truncComp <- function(formula, ...) {
  UseMethod("truncComp")
}

#' @rdname truncComp
#' @export
truncComp.formula <- function(formula, atom, data, method = c("LRT", "SPLRT"),
                              conf_level = 0.95, init = NULL, adjust = NULL, ...) {
  if(!inherits(formula, "formula")) {
    stop("The formula must be a formula.")
  }

  method <- match.arg(method)

  if(length(attr(terms(formula), "term.labels")) != 1) {
    stop("The current implementation must have one covariate in the formula.")
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

  if(!all(sort(unique(treatment)) == c(0,1))) {
    stop("The covariate must be binary (0/1) indicating the two treatments.")
  }

  if(!(length(atom) == 1 & is.numeric(atom))) {
    stop("atom must be a single numeric value.")
  }

  alive <- as.numeric(outcome != atom)

  if(all(alive[treatment == 0] == 0) | all(alive[treatment == 1] == 0)) {
    stop("Nothing has been observed in one of the groups. Cannot do estimation.")
  }

  if(all(alive == 1)) {
    stop("Everything seems to have been observed. You should use a different method.")
  }


  truncComp_core(outcome,
                 alive,
                 treatment,
                 method,
                 conf_level,
                 init,
                 adjust_data = adjustment$data,
                 adjust_formula = adjustment$formula,
                 atom = atom)
}

#' Default interface for [truncComp()]
#'
#' @param y Outcome vector for the default interface.
#' @param a Binary indicator for whether the continuous outcome is observed.
#' @param r Binary treatment indicator for the default interface.
#' @rdname truncComp
#' @exportS3Method truncComp default
truncComp.default <- function(y, a, r, method = c("LRT", "SPLRT"),
                              conf_level = 0.95, init = NULL,
                              adjust = NULL, atom = NULL, ...) {
  method <- match.arg(method)
  adjustment <- prepareDefaultAdjustment(adjust, length(y))
  atom <- resolveDefaultAtom(y, a, atom = atom)

  truncComp_core(y,
                 a,
                 r,
                 method,
                 conf_level,
                 init,
                 adjust_data = adjustment$data,
                 adjust_formula = adjustment$formula,
                 atom = atom)
}
