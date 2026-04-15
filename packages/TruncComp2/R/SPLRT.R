adjusted_SPLRT <- function(data, conf.level = 0.95, adjust = NULL, adjust_spec = NULL,
                           atom = NULL) {
  if(is.null(adjust_spec)) {
    adjust_spec <- adjustmentSpecification(adjust)
  }

  fits <- parametric_fit_models(data, adjust = adjust)
  glm_regular <- parametric_glm_is_regular(fits$bernoulli_null, tol = 1e-8) &&
    parametric_glm_is_regular(fits$bernoulli_alt, term = "R", tol = 1e-8)

  if(!glm_regular) {
    return(returnErrorData(
      "Adjusted semi-parametric SPLRT is not estimable under the supplied covariates.",
      method = "SPLRT",
      conf.level = conf.level,
      data = data,
      adjust = adjust_spec,
      atom = atom
    ))
  }

  ll_null <- parametric_loglik_value(fits$bernoulli_null)
  ll_alt <- parametric_loglik_value(fits$bernoulli_alt)
  alphaDelta <- exp(parametric_term_estimate(fits$bernoulli_alt, "R"))
  alphaDeltaCI <- parametric_term_interval(fits$bernoulli_alt, "R", conf.level, transform = exp)

  if(!all(is.finite(c(ll_null, ll_alt, alphaDelta, alphaDeltaCI)))) {
    return(returnErrorData(
      "Adjusted semi-parametric SPLRT is not estimable under the supplied covariates.",
      method = "SPLRT",
      conf.level = conf.level,
      data = data,
      adjust = adjust_spec,
      atom = atom
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
    return(returnErrorData(
      el_fit$error,
      method = "SPLRT",
      conf.level = conf.level,
      data = data,
      adjust = adjust_spec,
      atom = atom
    ))
  }

  alphaW <- parametric_clamp_statistic(2 * (ll_alt - ll_null))
  W <- parametric_clamp_statistic(el_fit$statistic + alphaW)

  newTruncComp2(
    muDelta = as.numeric(el_fit$estimate),
    muDeltaCI = as.numeric(el_fit$conf.int),
    alphaDelta = as.numeric(alphaDelta),
    alphaDeltaCI = as.numeric(alphaDeltaCI),
    Delta = NA_real_,
    DeltaCI = delta_na_interval(),
    DeltaMarginalCI = delta_na_interval(),
    DeltaProjectedCI = delta_na_interval(),
    DeltaProfileCI = delta_na_interval(),
    W = W,
    p = stats::pchisq(W, 2, lower.tail = FALSE),
    method = "SPLRT",
    conf.level = conf.level,
    success = TRUE,
    init = NULL,
    data = data,
    adjust = adjust_spec,
    atom = atom
  )
}

SPLRT <- function(data, conf.level = 0.95, adjust = NULL, adjust_spec = NULL, atom = NULL) {
  if(!is.null(adjust)) {
    return(adjusted_SPLRT(data, conf.level = conf.level, adjust = adjust, adjust_spec = adjust_spec,
                          atom = atom))
  }

  yAlive1 <- data[data$R == 0 & data$A == 1, "Y"]
  yAlive2 <- data[data$R == 1 & data$A == 1, "Y"]

  #Empirical likelihood ratio test
  ELRT <- el_mean_diff_fit(yAlive2, yAlive1, conf.level = conf.level)

  muDelta <- as.numeric(ELRT$estimate)
  muDeltaCI <- as.numeric(ELRT$conf.int)
  muW <- as.numeric(ELRT$statistic)
  #muP <- as.numeric(ELRT$p.value)

  data$R <- as.numeric(data$R) #beware
  
  #Fit logistic models
  m0 <- stats::glm(A ~ 1, family=stats::binomial(), data=data)
  m1 <- stats::glm(A ~ R, family=stats::binomial(), data=data)

  binomConfint <- suppressMessages(stats::confint(m1, level = conf.level))
  alphaDelta <- exp(stats::coef(m1)["R"])
  alphaDeltaCI <- as.numeric(exp(binomConfint["R",]))

  binomTest <- stats::anova(m0, m1, test="LRT")
  alphaW <- as.numeric(binomTest$Deviance[2])
  #alphaP <- binomTest$"Pr(>Chi)"[2]
  #W.binom <- 2 * (logLik(m1) - logLik(m0))

  #Joint likelihood ratio test
  W <- as.numeric(muW + alphaW) #Joint test statistic
  p <- 1 - stats::pchisq(W, 2)  #Joint p-value

  out <- newTruncComp2(muDelta = muDelta,
                       muDeltaCI = muDeltaCI,
                       alphaDelta = alphaDelta,
                       alphaDeltaCI = alphaDeltaCI,
                       Delta = NA_real_,
                       DeltaCI = delta_na_interval(),
                       DeltaMarginalCI = delta_na_interval(),
                       DeltaProjectedCI = delta_na_interval(),
                       DeltaProfileCI = delta_na_interval(),
                       W = W,
                       p = p,
                       method = "SPLRT",
                       conf.level = conf.level,
                       success = TRUE,
                       init = NULL,
                       data = data,
                       atom = atom)

  augmentDeltaInference(out)
}
