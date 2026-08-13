adjusted_SPLRT <- function(data, conf.level = 0.95, adjust = NULL, adjust_spec = NULL,
                           atom = NULL, call = NULL) {
  if(is.null(adjust_spec)) {
    adjust_spec <- adjustmentSpecification(adjust)
  }

  fits <- parametric_fit_models(data, adjust = adjust)
  glm_regular <- parametric_glm_is_regular(fits$bernoulli_null, tol = 1e-8) &&
    parametric_glm_is_regular(fits$bernoulli_alt, term = "R", tol = 1e-8)

  if(!glm_regular) {
    return(new_failed_trunc_comp_fit(
      "The adjusted semiparametric likelihood ratio test is not estimable under the supplied covariates.",
      method = "splrt",
      conf.level = conf.level,
      data = data,
      adjust = adjust_spec,
      adjust_formula = adjust,
      atom = atom,
      call = call
    ))
  }

  ll_null <- parametric_loglik_value(fits$bernoulli_null)
  ll_alt <- parametric_loglik_value(fits$bernoulli_alt)
  alphaDelta <- exp(parametric_term_estimate(fits$bernoulli_alt, "R"))
  alphaDeltaCI <- parametric_glm_profile_interval(
    fits$bernoulli_alt,
    "R",
    conf.level,
    transform = exp
  )

  if(!all(is.finite(c(ll_null, ll_alt, alphaDelta, alphaDeltaCI)))) {
    return(new_failed_trunc_comp_fit(
      "The adjusted semiparametric likelihood ratio test is not estimable under the supplied covariates.",
      method = "splrt",
      conf.level = conf.level,
      data = data,
      adjust = adjust_spec,
      adjust_formula = adjust,
      atom = atom,
      call = call
    ))
  }

  el_fit <- el_regression_fit(
    data = fits$normal_data,
    formula = fits$formulas$normal_alt,
    term = "R",
    mu = 0,
    conf.level = conf.level
  )

  if(!isTRUE(el_fit$success)) {
    return(new_failed_trunc_comp_fit(
      el_fit$error,
      method = "splrt",
      conf.level = conf.level,
      data = data,
      adjust = adjust_spec,
      adjust_formula = adjust,
      atom = atom,
      call = call
    ))
  }

  alphaW <- parametric_clamp_statistic(2 * (ll_alt - ll_null))
  W <- parametric_clamp_statistic(el_fit$statistic + alphaW)

  new_trunc_comp_fit(
    mu_delta = as.numeric(el_fit$estimate),
    mu_delta_ci = as.numeric(el_fit$conf.int),
    alpha_delta = as.numeric(alphaDelta),
    alpha_delta_ci = as.numeric(alphaDeltaCI),
    delta = NA_real_,
    statistic = W,
    p.value = stats::pchisq(W, 2, lower.tail = FALSE),
    method = "splrt",
    conf.level = conf.level,
    success = TRUE,
    data = data,
    adjust = adjust_spec,
    adjust_formula = adjust,
    atom = atom,
    call = call
  )
}

SPLRT <- function(data, conf.level = 0.95, adjust = NULL, adjust_spec = NULL,
                  atom = NULL, call = NULL) {
  if(!is.null(adjust)) {
    return(adjusted_SPLRT(
      data,
      conf.level = conf.level,
      adjust = adjust,
      adjust_spec = adjust_spec,
      atom = atom,
      call = call
    ))
  }

  y0 <- data[data$R == 0 & data$A == 1, "Y"]
  y1 <- data[data$R == 1 & data$A == 1, "Y"]

  ELRT <- el_mean_diff_fit(y1, y0, conf.level = conf.level)

  muDelta <- as.numeric(ELRT$estimate)
  muDeltaCI <- as.numeric(ELRT$conf.int)
  muW <- as.numeric(ELRT$statistic)

  data$R <- as.numeric(data$R)

  m0 <- stats::glm(A ~ 1, family = stats::binomial(), data = data)
  m1 <- stats::glm(A ~ R, family = stats::binomial(), data = data)

  alphaDelta <- exp(stats::coef(m1)["R"])
  alphaDeltaCI <- parametric_glm_profile_interval(
    m1,
    "R",
    conf.level,
    transform = exp
  )

  binomTest <- stats::anova(m0, m1, test = "LRT")
  alphaW <- as.numeric(binomTest$Deviance[2])

  W <- as.numeric(muW + alphaW)
  p.value <- stats::pchisq(W, 2, lower.tail = FALSE)

  out <- new_trunc_comp_fit(
    mu_delta = muDelta,
    mu_delta_ci = muDeltaCI,
    alpha_delta = alphaDelta,
    alpha_delta_ci = alphaDeltaCI,
    delta = NA_real_,
    statistic = W,
    p.value = p.value,
    method = "splrt",
    conf.level = conf.level,
    success = TRUE,
    data = data,
    atom = atom,
    call = call
  )

  augmentDeltaInference(out)
}
