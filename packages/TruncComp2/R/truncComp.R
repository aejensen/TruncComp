normalizeFormulaAdjust <- function(adjust, outcome_name, treatment_name) {
  if(is.null(adjust)) {
    return(NULL)
  }

  if(!inherits(adjust, "formula") || length(adjust) != 2) {
    stop("adjust must be NULL or a one-sided formula like ~ age + sex.")
  }

  adjust_terms <- stats::terms(adjust)
  labels <- attr(adjust_terms, "term.labels")
  if("." %in% labels) {
    stop("adjust must not use '.'. Please list adjustment covariates explicitly.")
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

truncComp_core <- function(y, a, r, method, conf.level = 0.95, init = NULL,
                           adjust_data = NULL, adjust_formula = NULL) {
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
                           adjust = adjust_spec))
  }

  if(!is.null(adjust_formula) && method == "SPLRT") {
    stop('Covariate adjustment is currently only implemented for method = "LRT".')
  }

  if(method == "LRT") {
    out <- LRT(d, init, conf.level, adjust = adjust_formula, adjust_spec = adjust_spec)
  } else if(method == "SPLRT") {
    out <- SPLRT(d, conf.level)
  }

  out$data <- d
  out
}

truncComp <- function(formula, atom, data, method, conf.level = 0.95, init = NULL,
                      adjust = NULL) {
  if(!inherits(formula, "formula")) {
    stop("The formula must be a formula.")
  }

  if(!(method == "LRT" | method == "SPLRT")) {
    stop("Only LRT or SPLRT supported as methods.")
  }

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
                 conf.level,
                 init,
                 adjust_data = adjustment$data,
                 adjust_formula = adjustment$formula)
}


truncComp.default <- function(y, a, r, method, conf.level = 0.95, init = NULL,
                              adjust = NULL) {
  adjustment <- prepareDefaultAdjustment(adjust, length(y))

  truncComp_core(y,
                 a,
                 r,
                 method,
                 conf.level,
                 init,
                 adjust_data = adjustment$data,
                 adjust_formula = adjustment$formula)
}
