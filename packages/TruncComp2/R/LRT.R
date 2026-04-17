parametric_safe_xlogy <- function(x, y) {
  out <- numeric(length(x))
  nonzero <- x != 0

  if(any(nonzero)) {
    out[nonzero] <- x[nonzero] * log(y[nonzero])
  }

  out
}

parametric_odds_ratio <- function(pi1, pi0) {
  if((pi1 == 0 && pi0 == 0) || (pi1 == 1 && pi0 == 1)) {
    return(1)
  }

  if(pi1 == 1 || pi0 == 0) {
    return(Inf)
  }

  if(pi1 == 0 || pi0 == 1) {
    return(0)
  }

  (pi1 / (1 - pi1)) / (pi0 / (1 - pi0))
}

parametric_wald_interval <- function(estimate, se, conf.level) {
  if(!(length(se) == 1 && is.finite(se) && se >= 0)) {
    return(c(NA_real_, NA_real_))
  }

  z <- stats::qnorm((1 + conf.level) / 2)
  estimate + c(-1, 1) * z * se
}

parametric_loglik_bernoulli <- function(k, n, p) {
  sum(parametric_safe_xlogy(c(k, n - k), c(p, 1 - p)))
}

parametric_model_formulas <- function(adjust = NULL) {
  if(is.null(adjust)) {
    return(list(
      bernoulli_null = A ~ 1,
      bernoulli_alt = A ~ R,
      normal_null = Y ~ 1,
      normal_alt = Y ~ R
    ))
  }

  list(
    bernoulli_null = stats::update(adjust, A ~ .),
    bernoulli_alt = stats::update(adjust, A ~ R + .),
    normal_null = stats::update(adjust, Y ~ .),
    normal_alt = stats::update(adjust, Y ~ R + .)
  )
}

parametric_fit_glm <- function(formula, data) {
  suppressWarnings(tryCatch(
    stats::glm(formula, family = stats::binomial(), data = data),
    error = function(e) NULL
  ))
}

parametric_fit_lm <- function(formula, data) {
  suppressWarnings(tryCatch(
    stats::lm(formula, data = data),
    error = function(e) NULL
  ))
}

parametric_fit_models <- function(data, adjust = NULL) {
  formulas <- parametric_model_formulas(adjust)
  observed_variables <- c("Y", "R")
  if(!is.null(adjust)) {
    observed_variables <- unique(c(observed_variables, all.vars(adjust)))
  }
  observed_data <- droplevels(data[data$A == 1, observed_variables, drop = FALSE])

  list(
    formulas = formulas,
    bernoulli_null = parametric_fit_glm(formulas$bernoulli_null, data),
    bernoulli_alt = parametric_fit_glm(formulas$bernoulli_alt, data),
    normal_data = observed_data,
    normal_null = parametric_fit_lm(formulas$normal_null, observed_data),
    normal_alt = parametric_fit_lm(formulas$normal_alt, observed_data)
  )
}

parametric_loglik_value <- function(fit, reml = NULL) {
  if(is.null(fit)) {
    return(NA_real_)
  }

  ll <- tryCatch(
    if(is.null(reml)) {
      stats::logLik(fit)
    } else {
      stats::logLik(fit, REML = reml)
    },
    error = function(e) NA_real_
  )

  as.numeric(ll)
}

parametric_clamp_statistic <- function(statistic, tol = 1e-12) {
  statistic <- as.numeric(statistic)

  if(is.finite(statistic) && statistic < tol) {
    return(0)
  }

  statistic
}

parametric_term_variance <- function(fit, term) {
  covariance <- tryCatch(stats::vcov(fit), error = function(e) NULL)

  if(is.null(covariance) || !(term %in% rownames(covariance))) {
    return(NA_real_)
  }

  as.numeric(covariance[term, term])
}

parametric_term_estimate <- function(fit, term) {
  coefficients <- tryCatch(stats::coef(fit), error = function(e) NULL)

  if(is.null(coefficients) || !(term %in% names(coefficients))) {
    return(NA_real_)
  }

  as.numeric(coefficients[[term]])
}

parametric_term_interval <- function(fit, term, conf.level, transform = identity) {
  estimate <- parametric_term_estimate(fit, term)
  variance <- parametric_term_variance(fit, term)

  if(!is.finite(estimate) || !is.finite(variance) || variance < 0) {
    return(c(NA_real_, NA_real_))
  }

  interval <- transform(parametric_wald_interval(estimate, sqrt(variance), conf.level))

  if(!all(is.finite(interval))) {
    return(c(NA_real_, NA_real_))
  }

  as.numeric(interval)
}

parametric_glm_is_regular <- function(fit, term = NULL, tol = 1e-8) {
  if(is.null(fit) || !isTRUE(fit$converged)) {
    return(FALSE)
  }

  design <- tryCatch(stats::model.matrix(fit), error = function(e) NULL)
  if(is.null(design) || fit$rank < ncol(design)) {
    return(FALSE)
  }

  ll <- parametric_loglik_value(fit)
  if(!is.finite(ll)) {
    return(FALSE)
  }

  fitted_values <- tryCatch(stats::fitted(fit), error = function(e) NULL)
  if(is.null(fitted_values) ||
     any(!is.finite(fitted_values)) ||
     any(fitted_values <= tol | fitted_values >= 1 - tol)) {
    return(FALSE)
  }

  if(is.null(term)) {
    return(TRUE)
  }

  estimate <- parametric_term_estimate(fit, term)
  variance <- parametric_term_variance(fit, term)

  is.finite(estimate) && is.finite(variance) && variance > tol
}

parametric_lm_is_regular <- function(fit, term = NULL, tol = 1e-8) {
  if(is.null(fit)) {
    return(FALSE)
  }

  design <- tryCatch(stats::model.matrix(fit), error = function(e) NULL)
  if(is.null(design) || fit$rank < ncol(design)) {
    return(FALSE)
  }

  ll <- parametric_loglik_value(fit, reml = FALSE)
  sigma <- tryCatch(summary(fit)$sigma, error = function(e) NA_real_)
  if(!is.finite(ll) || !is.finite(sigma) || sigma <= tol) {
    return(FALSE)
  }

  if(is.null(term)) {
    return(TRUE)
  }

  estimate <- parametric_term_estimate(fit, term)
  variance <- parametric_term_variance(fit, term)

  is.finite(estimate) && is.finite(variance) && variance > tol
}

parametric_bernoulli_summary <- function(a, r, fits, conf.level, tol = 1e-12) {
  n0 <- sum(r == 0)
  n1 <- sum(r == 1)
  k0 <- sum(a[r == 0])
  k1 <- sum(a[r == 1])

  pi0 <- k0 / n0
  pi1 <- k1 / n1
  pi <- (k0 + k1) / (n0 + n1)

  ell_h0_fallback <- parametric_loglik_bernoulli(k0 + k1, n0 + n1, pi)
  ell_ha_fallback <- parametric_loglik_bernoulli(k0, n0, pi0) +
    parametric_loglik_bernoulli(k1, n1, pi1)

  ll_h0 <- parametric_loglik_value(fits$bernoulli_null)
  ll_ha <- parametric_loglik_value(fits$bernoulli_alt)

  use_model_loglik <- !is.null(fits$bernoulli_null) &&
    !is.null(fits$bernoulli_alt) &&
    isTRUE(fits$bernoulli_null$converged) &&
    isTRUE(fits$bernoulli_alt$converged) &&
    all(is.finite(c(ll_h0, ll_ha)))

  W <- if(use_model_loglik) {
    2 * (ll_ha - ll_h0)
  } else {
    2 * (ell_ha_fallback - ell_h0_fallback)
  }
  W <- parametric_clamp_statistic(W, tol = tol)

  regular_cells <- all(c(k0, n0 - k0, k1, n1 - k1) > 0)
  coefficient <- parametric_term_estimate(fits$bernoulli_alt, "R")
  alpha_model <- exp(coefficient)
  alpha_exact <- parametric_odds_ratio(pi1, pi0)

  if(regular_cells && is.finite(alpha_model)) {
    alphaDelta <- alpha_model
    alphaDeltaCI <- parametric_term_interval(
      fits$bernoulli_alt,
      "R",
      conf.level,
      transform = exp
    )
  } else {
    alphaDelta <- alpha_exact
    alphaDeltaCI <- c(NA_real_, NA_real_)
  }

  list(
    n0 = n0,
    n1 = n1,
    k0 = k0,
    k1 = k1,
    pi0 = pi0,
    pi1 = pi1,
    pi = pi,
    ell_h0 = if(use_model_loglik) ll_h0 else ell_h0_fallback,
    ell_ha = if(use_model_loglik) ll_ha else ell_ha_fallback,
    W = W,
    alphaDelta = alphaDelta,
    alphaDeltaCI = alphaDeltaCI,
    fit_null = fits$bernoulli_null,
    fit_alt = fits$bernoulli_alt
  )
}

parametric_normal_summary <- function(y0, y1, fits, conf.level, tol = 1e-12) {
  m0 <- length(y0)
  m1 <- length(y1)
  m <- m0 + m1

  mu0 <- mean(y0)
  mu1 <- mean(y1)
  mu <- mean(c(y0, y1))

  sse1 <- sum((y0 - mu0)^2) + sum((y1 - mu1)^2)
  sse0 <- sum((c(y0, y1) - mu)^2)

  sse1 <- max(as.numeric(sse1), 0)
  sse0 <- max(as.numeric(sse0), 0)

  if(sse0 < sse1 && (sse1 - sse0) <= tol * max(1, sse1)) {
    sse0 <- sse1
  }

  ll_h0 <- parametric_loglik_value(fits$normal_null, reml = FALSE)
  ll_ha <- parametric_loglik_value(fits$normal_alt, reml = FALSE)

  if(sse1 <= tol) {
    W <- if(sse0 <= tol) 0 else Inf
    muDeltaCI <- c(NA_real_, NA_real_)
  } else if(all(is.finite(c(ll_h0, ll_ha)))) {
    W <- 2 * (ll_ha - ll_h0)
    W <- parametric_clamp_statistic(W, tol = tol)
    muDeltaCI <- parametric_term_interval(fits$normal_alt, "R", conf.level)
  } else {
    W <- parametric_clamp_statistic(m * log(sse0 / sse1), tol = tol)
    muDeltaCI <- c(NA_real_, NA_real_)
  }

  coefficient <- parametric_term_estimate(fits$normal_alt, "R")
  muDelta <- if(is.finite(coefficient)) coefficient else mu1 - mu0

  if(sse1 <= tol || !all(is.finite(muDeltaCI))) {
    muDeltaCI <- c(NA_real_, NA_real_)
  }

  list(
    m0 = m0,
    m1 = m1,
    m = m,
    mu0 = mu0,
    mu1 = mu1,
    mu = mu,
    sse0 = sse0,
    sse1 = sse1,
    sigma2_hat = if(sse1 <= tol) 0 else sse1 / m,
    ell_h0 = ll_h0,
    ell_ha = ll_ha,
    W = W,
    muDelta = muDelta,
    muDeltaCI = muDeltaCI,
    fit_null = fits$normal_null,
    fit_alt = fits$normal_alt
  )
}

parametric_adjusted_lrt_fit <- function(data, adjust, conf.level = 0.95,
                                        tol = 1e-12, regularity_tol = 1e-8) {
  fits <- parametric_fit_models(data, adjust = adjust)

  glm_regular <- parametric_glm_is_regular(fits$bernoulli_null, tol = regularity_tol) &&
    parametric_glm_is_regular(fits$bernoulli_alt, term = "R", tol = regularity_tol)
  lm_regular <- parametric_lm_is_regular(fits$normal_null, tol = regularity_tol) &&
    parametric_lm_is_regular(fits$normal_alt, term = "R", tol = regularity_tol)

  if(!(glm_regular && lm_regular)) {
    return(list(
      success = FALSE,
      error = "Adjusted parametric LRT is not estimable under the supplied covariates."
    ))
  }

  ll_glm_null <- parametric_loglik_value(fits$bernoulli_null)
  ll_glm_alt <- parametric_loglik_value(fits$bernoulli_alt)
  ll_lm_null <- parametric_loglik_value(fits$normal_null, reml = FALSE)
  ll_lm_alt <- parametric_loglik_value(fits$normal_alt, reml = FALSE)

  if(!all(is.finite(c(ll_glm_null, ll_glm_alt, ll_lm_null, ll_lm_alt)))) {
    return(list(
      success = FALSE,
      error = "Adjusted parametric LRT is not estimable under the supplied covariates."
    ))
  }

  W_A <- parametric_clamp_statistic(2 * (ll_glm_alt - ll_glm_null), tol = tol)
  W_Y <- parametric_clamp_statistic(2 * (ll_lm_alt - ll_lm_null), tol = tol)
  W <- parametric_clamp_statistic(W_A + W_Y, tol = tol)

  log_alpha_delta <- parametric_term_estimate(fits$bernoulli_alt, "R")
  mu_delta <- parametric_term_estimate(fits$normal_alt, "R")
  alpha_delta_ci <- parametric_term_interval(fits$bernoulli_alt, "R", conf.level, transform = exp)
  mu_delta_ci <- parametric_term_interval(fits$normal_alt, "R", conf.level)

  if(!is.finite(log_alpha_delta) ||
     !is.finite(mu_delta) ||
     !all(is.finite(alpha_delta_ci)) ||
     !all(is.finite(mu_delta_ci))) {
    return(list(
      success = FALSE,
      error = "Adjusted parametric LRT is not estimable under the supplied covariates."
    ))
  }

  list(
    success = TRUE,
    mu_delta = mu_delta,
    mu_delta_ci = mu_delta_ci,
    alpha_delta = exp(log_alpha_delta),
    alpha_delta_ci = alpha_delta_ci,
    delta = NA_real_,
    W_A = W_A,
    W_Y = W_Y,
    W = W,
    statistic_A = W_A,
    statistic_Y = W_Y,
    statistic = W,
    p.value = stats::pchisq(W, df = 2, lower.tail = FALSE),
    bernoulli = list(
      fit_null = fits$bernoulli_null,
      fit_alt = fits$bernoulli_alt,
      ell_h0 = ll_glm_null,
      ell_ha = ll_glm_alt
    ),
    normal = list(
      fit_null = fits$normal_null,
      fit_alt = fits$normal_alt,
      ell_h0 = ll_lm_null,
      ell_ha = ll_lm_alt
    )
  )
}

parametric_lrt_fit <- function(data, conf.level = 0.95, adjust = NULL,
                               tol = 1e-12, regularity_tol = 1e-8) {
  if(!is.null(adjust)) {
    return(parametric_adjusted_lrt_fit(
      data,
      adjust = adjust,
      conf.level = conf.level,
      tol = tol,
      regularity_tol = regularity_tol
    ))
  }

  y0 <- data$Y[data$R == 0 & data$A == 1]
  y1 <- data$Y[data$R == 1 & data$A == 1]
  fits <- parametric_fit_models(data)

  bernoulli <- parametric_bernoulli_summary(data$A, data$R, fits, conf.level, tol = tol)
  normal <- parametric_normal_summary(y0, y1, fits, conf.level, tol = tol)

  W <- parametric_clamp_statistic(bernoulli$W + normal$W, tol = tol)

  list(
    success = TRUE,
    mu_delta = normal$muDelta,
    mu_delta_ci = normal$muDeltaCI,
    alpha_delta = bernoulli$alphaDelta,
    alpha_delta_ci = bernoulli$alphaDeltaCI,
    delta = NA_real_,
    W_A = bernoulli$W,
    W_Y = normal$W,
    W = W,
    statistic_A = bernoulli$W,
    statistic_Y = normal$W,
    statistic = W,
    p.value = stats::pchisq(W, df = 2, lower.tail = FALSE),
    bernoulli = bernoulli,
    normal = normal
  )
}

LRT <- function(data, init = NULL, conf.level = 0.95, adjust = NULL, adjust_spec = NULL,
                atom = NULL, call = NULL) {
  if(is.null(adjust_spec)) {
    adjust_spec <- adjustmentSpecification(adjust)
  }

  fit <- parametric_lrt_fit(data, conf.level = conf.level, adjust = adjust)

  if(!isTRUE(fit$success)) {
    return(new_failed_trunc_comp_fit(
      fit$error,
      method = "lrt",
      conf.level = conf.level,
      init = init,
      data = data,
      adjust = adjust_spec,
      atom = atom,
      call = call
    ))
  }

  out <- new_trunc_comp_fit(
    mu_delta = fit$mu_delta,
    mu_delta_ci = fit$mu_delta_ci,
    alpha_delta = fit$alpha_delta,
    alpha_delta_ci = fit$alpha_delta_ci,
    delta = fit$delta,
    statistic = fit$statistic,
    p.value = fit$p.value,
    method = "lrt",
    conf.level = conf.level,
    success = TRUE,
    init = init,
    data = data,
    adjust = adjust_spec,
    atom = atom,
    call = call
  )

  augmentDeltaInference(out)
}
