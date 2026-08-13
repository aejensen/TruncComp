isDataOkay <- function(d) {
  non_atom_data <- subset(d, d$A == 1)

  y0 <- non_atom_data$Y[non_atom_data$R == 0]
  y1 <- non_atom_data$Y[non_atom_data$R == 1]

  if(length(y0) < 2 || length(y1) < 2) {
    return(FALSE)
  }

  TRUE
}

adjustmentSpecification <- function(adjust) {
  if(is.null(adjust)) {
    return(NULL)
  }

  labels <- attr(stats::terms(adjust), "term.labels")
  if(length(labels) == 0) {
    return(NULL)
  }

  paste(labels, collapse = " + ")
}

is_valid_trunc_comp_fit <- function(object) {
  isTRUE(object$success)
}

normalize_trunc_comp_method <- function(method) {
  match.arg(tolower(method), c("lrt", "splrt"))
}

trunc_comp_method_label <- function(method) {
  method <- normalize_trunc_comp_method(method)

  switch(
    method,
    lrt = "Parametric likelihood ratio test",
    splrt = "Semiparametric likelihood ratio test"
  )
}

numeric_or_na <- function(x) {
  if(is.null(x) || length(x) == 0) {
    return(NA_real_)
  }

  as.numeric(x[[1]])
}

new_trunc_comp_fit <- function(mu_delta = NULL, mu_delta_ci = NULL,
                               alpha_delta = NULL, alpha_delta_ci = NULL,
                               delta = NULL,
                               statistic = NULL, p.value = NULL,
                               method, conf.level, success,
                               error = "", data = NULL,
                               adjust = NULL, adjust_formula = NULL,
                               atom = NULL, call = NULL) {
  method <- normalize_trunc_comp_method(method)

  out <- list(
    mu_delta = mu_delta,
    mu_delta_ci = mu_delta_ci,
    alpha_delta = alpha_delta,
    alpha_delta_ci = alpha_delta_ci,
    delta = delta,
    statistic = statistic,
    p.value = p.value,
    method = method,
    conf.level = conf.level,
    success = success,
    error = error,
    data = data,
    adjust = adjust,
    adjust_formula = adjust_formula,
    atom = atom,
    call = call
  )
  class(out) <- c("trunc_comp_fit", "list")
  out
}

new_failed_trunc_comp_fit <- function(error, method, conf.level, data = NULL,
                                      adjust = NULL,
                                      adjust_formula = NULL, atom = NULL,
                                      call = NULL) {
  new_trunc_comp_fit(
    method = method,
    conf.level = conf.level,
    success = FALSE,
    error = error,
    data = data,
    adjust = adjust,
    adjust_formula = adjust_formula,
    atom = atom,
    call = call
  )
}
