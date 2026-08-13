delta_na_interval <- function() {
  c(NA_real_, NA_real_)
}

delta_from_components <- function(atom, mu0, mu1, pi0, pi1) {
  combined1 <- pi1 * mu1 + (1 - pi1) * atom
  combined0 <- pi0 * mu0 + (1 - pi0) * atom
  as.numeric(combined1 - combined0)
}

delta_welch_interval <- function(y, a, r, atom, conf.level) {
  y <- ifelse(a == 1, y, atom)
  y0 <- y[r == 0]
  y1 <- y[r == 1]

  delta <- mean(y1) - mean(y0)
  v0 <- stats::var(y0)
  v1 <- stats::var(y1)
  n0 <- length(y0)
  n1 <- length(y1)

  if(!is.finite(v0)) {
    v0 <- 0
  }
  if(!is.finite(v1)) {
    v1 <- 0
  }

  se2 <- v1 / n1 + v0 / n0
  if(!is.finite(se2) || se2 <= 0) {
    return(rep(as.numeric(delta), 2))
  }

  numerator <- se2 ^ 2
  denominator <- 0
  if(v1 > 0) {
    denominator <- denominator + (v1 / n1) ^ 2 / (n1 - 1)
  }
  if(v0 > 0) {
    denominator <- denominator + (v0 / n0) ^ 2 / (n0 - 1)
  }

  df <- if(denominator > 0) numerator / denominator else Inf
  crit <- if(is.finite(df) && df > 0) {
    stats::qt((1 + conf.level) / 2, df = df)
  } else {
    stats::qnorm((1 + conf.level) / 2)
  }

  as.numeric(delta + c(-1, 1) * crit * sqrt(se2))
}

delta_unadjusted_point_estimate <- function(object, tol = 1e-8) {
  y0 <- object$data$Y[object$data$R == 0 & object$data$A == 1]
  y1 <- object$data$Y[object$data$R == 1 & object$data$A == 1]
  pi0 <- mean(object$data$A[object$data$R == 0])
  pi1 <- mean(object$data$A[object$data$R == 1])

  if(identical(object$method, "splrt")) {
    theta_fit <- el_mean_diff_theta(object$mu_delta, y1, y0, tol = tol)
    if(isTRUE(theta_fit$feasible) && is.finite(theta_fit$theta)) {
      mu1 <- theta_fit$theta
      mu0 <- theta_fit$theta - object$mu_delta
      return(delta_from_components(object$atom, mu0, mu1, pi0, pi1))
    }
  }

  delta_from_components(object$atom, mean(y0), mean(y1), pi0, pi1)
}

prepareParametricJointReference <- function(data, adjust = NULL) {
  fits <- parametric_fit_models(data, adjust = adjust)
  observed_data <- fits$normal_data
  glm_alt <- fits$bernoulli_alt
  lm_alt <- fits$normal_alt

  list(
    data = data,
    adjust = adjust,
    adjusted = !is.null(adjust),
    observed_data = observed_data,
    glm_alt = glm_alt,
    lm_alt = lm_alt,
    ll_glm_alt = parametric_loglik_value(glm_alt),
    ll_lm_alt = parametric_loglik_value(lm_alt, reml = FALSE),
    observed_mean = mean(observed_data$Y),
    observed_treatment_rate = mean(observed_data$R),
    muDelta_hat = parametric_term_estimate(lm_alt, "R")
  )
}

parametricJointCandidate <- function(parametricReference, muDelta, logORdelta, tol = 1e-12) {
  W_A <- parametric_glm_profile_statistic(
    parametricReference$glm_alt,
    "R",
    logORdelta,
    ll_unrestricted = parametricReference$ll_glm_alt,
    tol = tol
  )
  lm_profile <- parametric_lm_profile_fit(
    parametricReference$lm_alt,
    "R",
    muDelta,
    ll_unrestricted = parametricReference$ll_lm_alt,
    tol = tol
  )
  W_Y <- lm_profile$statistic
  lm_constrained <- lm_profile$fit

  beta0 <- NA_real_
  mu0 <- parametric_term_estimate(lm_constrained, "(Intercept)")
  p0 <- NA_real_
  p1 <- NA_real_
  if(!isTRUE(parametricReference$adjusted)) {
    glm_terms <- attr(stats::terms(parametricReference$glm_alt), "term.labels")
    glm_formula <- stats::reformulate(
      setdiff(glm_terms, "R"),
      response = "A",
      intercept = attr(stats::terms(parametricReference$glm_alt), "intercept")
    )
    glm_constrained <- suppressWarnings(tryCatch(
      stats::glm(
        glm_formula,
        family = stats::binomial(),
        data = parametricReference$glm_alt$model,
        offset = logORdelta * parametricReference$glm_alt$model$R
      ),
      error = function(e) NULL
    ))
    beta0 <- parametric_term_estimate(glm_constrained, "(Intercept)")
    p0 <- if(is.finite(beta0)) stats::plogis(beta0) else NA_real_
    p1 <- if(is.finite(beta0)) stats::plogis(beta0 + logORdelta) else NA_real_
  } else {
    glm_constrained <- NULL
  }
  mu1 <- if(is.finite(mu0)) mu0 + muDelta else NA_real_
  total <- W_A + W_Y
  if(is.finite(total)) {
    total <- parametric_clamp_statistic(total, tol = tol)
  }

  list(
    statistic = total,
    W_A = W_A,
    W_Y = W_Y,
    beta0 = beta0,
    mu0 = mu0,
    mu1 = mu1,
    p0 = p0,
    p1 = p1,
    glm = glm_constrained,
    lm = lm_constrained
  )
}

prepareSPLRTJointReference <- function(data, adjust = NULL) {
  if(is.null(adjust)) {
    return(list(
      adjusted = FALSE,
      yAlive0 = data$Y[data$R == 0 & data$A == 1],
      yAlive1 = data$Y[data$R == 1 & data$A == 1],
      logitReference = logit.prepare(data)
    ))
  }

  fits <- parametric_fit_models(data, adjust = adjust)
  list(
    adjusted = TRUE,
    glm_alt = fits$bernoulli_alt,
    ll_glm_alt = parametric_loglik_value(fits$bernoulli_alt),
    el_design = el_regression_design(
      fits$formulas$normal_alt,
      fits$normal_data,
      term = "R"
    )
  )
}

splrtContinuousCandidate <- function(splrtReference, muDelta, tol = 1e-8) {
  theta_fit <- el_mean_diff_theta(muDelta, splrtReference$yAlive1, splrtReference$yAlive0, tol = tol)

  if(!isTRUE(theta_fit$feasible) || !is.finite(theta_fit$theta)) {
    return(list(
      feasible = FALSE,
      statistic = Inf,
      theta = NA_real_,
      mu0 = NA_real_,
      mu1 = NA_real_
    ))
  }

  mu1 <- theta_fit$theta
  mu0 <- theta_fit$theta - muDelta
  list(
    feasible = TRUE,
    statistic = as.numeric(theta_fit$statistic),
    theta = theta_fit$theta,
    mu0 = mu0,
    mu1 = mu1
  )
}

splrtJointCandidate <- function(splrtReference, muDelta, logORdelta, tol = 1e-8) {
  if(isTRUE(splrtReference$adjusted)) {
    W_A <- parametric_glm_profile_statistic(
      splrtReference$glm_alt,
      "R",
      logORdelta,
      ll_unrestricted = splrtReference$ll_glm_alt,
      tol = tol
    )
    continuous <- el_regression_profile_fit(
      splrtReference$el_design,
      muDelta,
      tol = tol
    )
    W_Y <- if(isTRUE(continuous$success)) {
      as.numeric(continuous$statistic)
    } else {
      Inf
    }
    total <- W_A + W_Y
    if(is.finite(total)) {
      total <- parametric_clamp_statistic(total, tol = tol)
    }

    return(list(
      statistic = as.numeric(total),
      W_A = as.numeric(W_A),
      W_Y = as.numeric(W_Y),
      beta0 = NA_real_,
      p0 = NA_real_,
      p1 = NA_real_,
      theta = NA_real_,
      mu0 = NA_real_,
      mu1 = NA_real_,
      beta = if(isTRUE(continuous$success)) continuous$beta else NULL
    ))
  }

  logit_fit <- logit_profile_fit.prepared(splrtReference$logitReference, logORdelta)
  continuous <- splrtContinuousCandidate(
    splrtReference,
    muDelta,
    tol = tol
  )

  total <- as.numeric(logit_fit$statistic + continuous$statistic)
  if(is.finite(total) && total < tol) {
    total <- 0
  }

  list(
    statistic = total,
    W_A = as.numeric(logit_fit$statistic),
    W_Y = continuous$statistic,
    beta0 = logit_fit$beta0,
    p0 = logit_fit$p0,
    p1 = logit_fit$p1,
    theta = continuous$theta,
    mu0 = continuous$mu0,
    mu1 = continuous$mu1
  )
}

delta_full_reference <- function(object, tol = 1e-10) {
  if(!isTRUE(object$success) || !is.null(object$adjust) || is.null(object$atom) ||
     !(object$method %in% c("lrt", "splrt"))) {
    return(NULL)
  }

  data <- object$data
  y0 <- data$Y[data$R == 0 & data$A == 1]
  y1 <- data$Y[data$R == 1 & data$A == 1]
  n0 <- sum(data$R == 0)
  n1 <- sum(data$R == 1)
  k0 <- sum(data$A[data$R == 0])
  k1 <- sum(data$A[data$R == 1])

  if(length(y0) < 2 || length(y1) < 2 || n0 <= 0 || n1 <= 0 ||
     any(c(k0, n0 - k0, k1, n1 - k1) <= 0)) {
    return(NULL)
  }

  p0_hat <- k0 / n0
  p1_hat <- k1 / n1
  eta0_hat <- stats::qlogis(p0_hat)
  eta1_hat <- stats::qlogis(p1_hat)
  mu0_hat <- mean(y0)
  mu1_hat <- mean(y1)
  sse_hat <- sum((y0 - mu0_hat) ^ 2) + sum((y1 - mu1_hat) ^ 2)

  if(identical(object$method, "lrt") && (!is.finite(sse_hat) || sse_hat <= tol)) {
    return(NULL)
  }
  if(identical(object$method, "splrt") &&
     (stats::var(y0) <= tol || stats::var(y1) <= tol)) {
    return(NULL)
  }

  bernoulli_nll <- function(eta0, eta1) {
    n0 * logit_log1pexp(eta0) - k0 * eta0 +
      n1 * logit_log1pexp(eta1) - k1 * eta1
  }

  center <- c(mu0_hat, mu1_hat, eta0_hat, eta1_hat)
  list(
    method = object$method,
    atom = object$atom,
    y0 = y0,
    y1 = y1,
    n0 = n0,
    n1 = n1,
    k0 = k0,
    k1 = k1,
    p0_hat = p0_hat,
    p1_hat = p1_hat,
    eta0_hat = eta0_hat,
    eta1_hat = eta1_hat,
    mu0_hat = mu0_hat,
    mu1_hat = mu1_hat,
    sse_hat = sse_hat,
    n_observed = length(y0) + length(y1),
    bernoulli_nll = bernoulli_nll,
    bernoulli_nll_hat = bernoulli_nll(eta0_hat, eta1_hat),
    center = center,
    delta_hat = delta_from_components(object$atom, mu0_hat, mu1_hat, p0_hat, p1_hat)
  )
}

delta_full_candidate <- function(reference, par, tol = 1e-10) {
  par <- as.numeric(par)
  if(length(par) != 4 || any(!is.finite(par))) {
    return(list(feasible = FALSE, statistic = Inf, Delta = NA_real_))
  }

  mu0 <- par[1]
  mu1 <- par[2]
  eta0 <- par[3]
  eta1 <- par[4]
  p0 <- stats::plogis(eta0)
  p1 <- stats::plogis(eta1)

  W_A <- 2 * (reference$bernoulli_nll(eta0, eta1) - reference$bernoulli_nll_hat)
  W_A <- parametric_clamp_statistic(W_A, tol = tol)

  if(identical(reference$method, "lrt")) {
    sse <- sum((reference$y0 - mu0) ^ 2) + sum((reference$y1 - mu1) ^ 2)
    W_Y <- if(is.finite(sse) && sse > 0 && reference$sse_hat > 0) {
      reference$n_observed * log(sse / reference$sse_hat)
    } else {
      Inf
    }
    W_Y <- parametric_clamp_statistic(W_Y, tol = tol)
    feasible <- is.finite(W_Y)
  } else {
    fit0 <- el_one_sample_lambda(reference$y0 - mu0)
    fit1 <- el_one_sample_lambda(reference$y1 - mu1)
    feasible <- isTRUE(fit0$feasible) && isTRUE(fit1$feasible)
    W_Y <- if(feasible) fit0$statistic + fit1$statistic else Inf
  }

  total <- W_A + W_Y
  if(is.finite(total)) {
    total <- parametric_clamp_statistic(total, tol = tol)
  }

  list(
    feasible = feasible && is.finite(W_A) && is.finite(total),
    statistic = as.numeric(total),
    W_A = as.numeric(W_A),
    W_Y = as.numeric(W_Y),
    mu0 = mu0,
    mu1 = mu1,
    muDelta = mu1 - mu0,
    eta0 = eta0,
    eta1 = eta1,
    logORdelta = eta1 - eta0,
    p0 = p0,
    p1 = p1,
    Delta = delta_from_components(reference$atom, mu0, mu1, p0, p1)
  )
}

delta_bernoulli_arm_deviance <- function(reference, arm, eta) {
  n <- if(arm == 0L) reference$n0 else reference$n1
  k <- if(arm == 0L) reference$k0 else reference$k1
  eta_hat <- if(arm == 0L) reference$eta0_hat else reference$eta1_hat
  nll <- n * logit_log1pexp(eta) - k * eta
  nll_hat <- n * logit_log1pexp(eta_hat) - k * eta_hat
  as.numeric(2 * (nll - nll_hat))
}

delta_logit_sublevel_bounds <- function(reference, arm, threshold,
                                        tol = 1e-10) {
  n <- if(arm == 0L) reference$n0 else reference$n1
  k <- if(arm == 0L) reference$k0 else reference$k1
  eta_hat <- if(arm == 0L) reference$eta0_hat else reference$eta1_hat
  nll_hat <- n * logit_log1pexp(eta_hat) - k * eta_hat
  lower_bracket <- -(nll_hat + threshold / 2) / k
  upper_bracket <- (nll_hat + threshold / 2) / (n - k)
  objective <- function(eta) {
    delta_bernoulli_arm_deviance(reference, arm, eta) - threshold
  }

  lower <- suppressWarnings(tryCatch(
    stats::uniroot(objective, c(lower_bracket, eta_hat), tol = tol)$root,
    error = function(e) NA_real_
  ))
  upper <- suppressWarnings(tryCatch(
    stats::uniroot(objective, c(eta_hat, upper_bracket), tol = tol)$root,
    error = function(e) NA_real_
  ))
  c(lower, upper)
}

delta_parametric_eta_deviance <- function(reference, eta) {
  delta_bernoulli_arm_deviance(reference, 0L, eta[1]) +
    delta_bernoulli_arm_deviance(reference, 1L, eta[2])
}

delta_parametric_eta_candidate <- function(reference, eta, threshold,
                                           direction = c("lower", "upper"),
                                           tol = 1e-10) {
  direction <- match.arg(direction)
  eta <- as.numeric(eta)
  W_A <- delta_parametric_eta_deviance(reference, eta)
  remaining <- threshold - W_A
  if(length(eta) != 2 || any(!is.finite(eta)) || !is.finite(remaining) ||
     remaining < -tol) {
    return(list(feasible = FALSE, Delta = NA_real_, par = rep(NA_real_, 4)))
  }
  remaining <- max(remaining, 0)

  p0 <- stats::plogis(eta[1])
  p1 <- stats::plogis(eta[2])
  m0 <- length(reference$y0)
  m1 <- length(reference$y1)
  mean_radius_squared <- reference$sse_hat *
    expm1(remaining / reference$n_observed)
  coefficient_norm <- sqrt(p0 ^ 2 / m0 + p1 ^ 2 / m1)
  if(!is.finite(mean_radius_squared) || mean_radius_squared < 0 ||
     !is.finite(coefficient_norm) || coefficient_norm <= 0) {
    return(list(feasible = FALSE, Delta = NA_real_, par = rep(NA_real_, 4)))
  }

  sign_direction <- if(identical(direction, "upper")) 1 else -1
  root_radius <- sqrt(mean_radius_squared)
  mu0 <- reference$mu0_hat - sign_direction *
    root_radius * p0 / (m0 * coefficient_norm)
  mu1 <- reference$mu1_hat + sign_direction *
    root_radius * p1 / (m1 * coefficient_norm)
  par <- c(mu0, mu1, eta)

  list(
    feasible = TRUE,
    Delta = delta_from_components(reference$atom, mu0, mu1, p0, p1),
    par = par,
    W_A = W_A,
    W_Y = remaining
  )
}

delta_parametric_radial_limit <- function(reference, threshold, phi,
                                          eta_bounds, tol = 1e-10) {
  phi <- phi %% (2 * pi)
  unit <- c(cos(phi), sin(phi))
  center <- c(reference$eta0_hat, reference$eta1_hat)
  radial_limits <- rep(Inf, 2)
  for(j in seq_len(2)) {
    if(unit[j] > tol) {
      radial_limits[j] <- (eta_bounds$upper[j] - center[j]) / unit[j]
    } else if(unit[j] < -tol) {
      radial_limits[j] <- (eta_bounds$lower[j] - center[j]) / unit[j]
    }
  }
  outer <- min(radial_limits)
  if(!is.finite(outer) || outer <= 0) return(NA_real_)

  objective <- function(radius) {
    delta_parametric_eta_deviance(reference, center + radius * unit) - threshold
  }
  outer_value <- objective(outer)
  if(!is.finite(outer_value) || outer_value < -sqrt(tol)) return(NA_real_)
  if(abs(outer_value) <= tol) return(as.numeric(outer))

  suppressWarnings(tryCatch(
    stats::uniroot(objective, c(0, outer), tol = tol)$root,
    error = function(e) NA_real_
  ))
}

delta_parametric_radial_extreme <- function(reference, threshold, phi,
                                            eta_bounds,
                                            direction = c("lower", "upper"),
                                            radial_resolution = 65L,
                                            tol = 1e-10) {
  direction <- match.arg(direction)
  phi <- phi %% (2 * pi)
  unit <- c(cos(phi), sin(phi))
  center <- c(reference$eta0_hat, reference$eta1_hat)
  radial_limit <- delta_parametric_radial_limit(
    reference, threshold, phi, eta_bounds, tol = tol
  )
  if(!is.finite(radial_limit)) return(NULL)

  evaluate <- function(radius) {
    delta_parametric_eta_candidate(
      reference,
      center + radius * unit,
      threshold,
      direction = direction,
      tol = tol
    )
  }
  radii <- seq(0, radial_limit, length.out = radial_resolution)
  values <- vapply(radii, function(radius) evaluate(radius)$Delta, numeric(1))
  if(any(!is.finite(values))) return(NULL)
  best_index <- if(identical(direction, "lower")) which.min(values) else which.max(values)
  candidate_radii <- radii[c(1, best_index, length(radii))]

  if(best_index > 1 && best_index < length(radii)) {
    interval <- radii[c(best_index - 1, best_index + 1)]
    objective <- function(radius) evaluate(radius)$Delta
    optimum <- suppressWarnings(tryCatch(
      stats::optimize(
        objective,
        interval = interval,
        maximum = identical(direction, "upper"),
        tol = tol
      ),
      error = function(e) NULL
    ))
    if(!is.null(optimum)) {
      optimum_radius <- if(identical(direction, "upper")) {
        optimum$maximum
      } else {
        optimum$minimum
      }
      candidate_radii <- c(candidate_radii, optimum_radius)
    }
  }

  candidates <- lapply(candidate_radii, evaluate)
  candidate_values <- vapply(candidates, `[[`, numeric(1), "Delta")
  best <- if(identical(direction, "lower")) {
    which.min(candidate_values)
  } else {
    which.max(candidate_values)
  }
  candidates[[best]]
}

delta_parametric_sublevel_endpoint <- function(reference, threshold,
                                               direction = c("lower", "upper"),
                                               angle_resolution = 1440L,
                                               radial_resolution = 65L,
                                               tol = 1e-9) {
  direction <- match.arg(direction)
  eta0_bounds <- delta_logit_sublevel_bounds(reference, 0L, threshold, tol = tol)
  eta1_bounds <- delta_logit_sublevel_bounds(reference, 1L, threshold, tol = tol)
  eta_bounds <- list(
    lower = c(eta0_bounds[1], eta1_bounds[1]),
    upper = c(eta0_bounds[2], eta1_bounds[2])
  )
  if(any(!is.finite(c(eta_bounds$lower, eta_bounds$upper)))) return(NULL)

  angle_step <- 2 * pi / angle_resolution
  angles <- seq(0, 2 * pi - angle_step, length.out = angle_resolution)
  candidates <- lapply(angles, function(phi) {
    delta_parametric_radial_extreme(
      reference, threshold, phi, eta_bounds,
      direction = direction,
      radial_resolution = radial_resolution,
      tol = tol
    )
  })
  values <- vapply(candidates, function(candidate) {
    if(is.null(candidate)) NA_real_ else candidate$Delta
  }, numeric(1))
  if(!any(is.finite(values))) return(NULL)
  best_index <- if(identical(direction, "lower")) {
    which.min(replace(values, !is.finite(values), Inf))
  } else {
    which.max(replace(values, !is.finite(values), -Inf))
  }

  evaluate_angle <- function(phi) {
    delta_parametric_radial_extreme(
      reference, threshold, phi, eta_bounds,
      direction = direction,
      radial_resolution = radial_resolution,
      tol = tol
    )
  }
  angle_center <- angles[best_index]
  angle_objective <- function(phi) {
    candidate <- evaluate_angle(phi)
    if(is.null(candidate)) {
      return(if(identical(direction, "lower")) Inf else -Inf)
    }
    candidate$Delta
  }
  angle_optimum <- suppressWarnings(tryCatch(
    stats::optimize(
      angle_objective,
      interval = angle_center + c(-angle_step, angle_step),
      maximum = identical(direction, "upper"),
      tol = tol
    ),
    error = function(e) NULL
  ))
  finalists <- list(candidates[[best_index]])
  if(!is.null(angle_optimum)) {
    optimum_angle <- if(identical(direction, "upper")) {
      angle_optimum$maximum
    } else {
      angle_optimum$minimum
    }
    finalists[[2]] <- evaluate_angle(optimum_angle)
  }
  finalist_values <- vapply(finalists, function(candidate) candidate$Delta, numeric(1))
  best <- if(identical(direction, "lower")) {
    which.min(finalist_values)
  } else {
    which.max(finalist_values)
  }

  full_candidate <- delta_full_candidate(reference, finalists[[best]]$par, tol = tol)
  statistic_tolerance <- max(1e-6, 100 * tol)
  if(!isTRUE(full_candidate$feasible) ||
     abs(full_candidate$statistic - threshold) > statistic_tolerance) {
    return(NULL)
  }
  list(par = finalists[[best]]$par, candidate = full_candidate)
}

delta_parametric_sublevel_interval <- function(reference, threshold, tol = 1e-9) {
  lower <- delta_parametric_sublevel_endpoint(
    reference, threshold, direction = "lower", tol = tol
  )
  upper <- delta_parametric_sublevel_endpoint(
    reference, threshold, direction = "upper", tol = tol
  )
  if(is.null(lower) || is.null(upper)) return(delta_na_interval())
  as.numeric(c(lower$candidate$Delta, upper$candidate$Delta))
}

delta_splrt_sublevel_interval <- function(reference, threshold, tol = 1e-8) {
  coded0 <- c(rep(reference$atom, reference$n0 - reference$k0), reference$y0)
  coded1 <- c(rep(reference$atom, reference$n1 - reference$k1), reference$y1)
  estimate <- mean(coded1) - mean(coded0)
  delta_range <- el_mean_diff_delta_range(coded1, coded0)
  lower <- el_find_mean_diff_bound(
    coded1, coded0, estimate, delta_range[1], threshold, "lower", tol = tol
  )
  upper <- el_find_mean_diff_bound(
    coded1, coded0, estimate, delta_range[2], threshold, "upper", tol = tol
  )
  as.numeric(c(lower, upper))
}

delta_sublevel_range <- function(object, conf.level = object$conf.level,
                                 df, tol = 1e-7) {
  conf.level <- validateConfidenceLevel(conf.level)
  reference <- delta_full_reference(object, tol = tol)
  if(is.null(reference)) return(delta_na_interval())

  threshold <- stats::qchisq(conf.level, df)
  interval <- if(identical(reference$method, "lrt")) {
    delta_parametric_sublevel_interval(reference, threshold, tol = tol)
  } else {
    delta_splrt_sublevel_interval(reference, threshold, tol = tol)
  }
  if(any(!is.finite(interval))) {
    warning("Coded scale likelihood optimization did not converge.")
    return(delta_na_interval())
  }

  if(interval[1] > interval[2] ||
     interval[1] > reference$delta_hat + tol ||
     interval[2] < reference$delta_hat - tol) {
    warning("Coded scale likelihood optimization failed its feasibility checks.")
    return(delta_na_interval())
  }
  as.numeric(interval)
}

delta_projected_interval.optimize <- function(object,
                                              conf.level = object$conf.level,
                                              tol = 1e-7) {
  delta_sublevel_range(object, conf.level = conf.level, df = 4, tol = tol)
}

delta_profile_interval.optimize <- function(object,
                                            conf.level = object$conf.level,
                                            tol = 1e-7) {
  delta_sublevel_range(object, conf.level = conf.level, df = 1, tol = tol)
}

delta_projected_interval <- function(object, conf.level = object$conf.level,
                                     tol = 1e-8) {
  delta_projected_interval.optimize(
    object,
    conf.level = conf.level,
    tol = tol
  )
}

delta_profile_interval <- function(object, conf.level = object$conf.level,
                                   tol = 1e-8) {
  delta_profile_interval.optimize(
    object,
    conf.level = conf.level,
    tol = tol
  )
}

augmentDeltaInference <- function(object) {
  if(!isTRUE(object$success) || !is.null(object$adjust) || is.null(object$atom)) {
    object$delta <- if(is.null(object$delta)) NA_real_ else object$delta
    return(object)
  }

  object$delta <- delta_unadjusted_point_estimate(object)
  object
}
