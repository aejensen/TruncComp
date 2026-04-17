delta_na_interval <- function() {
  c(NA_real_, NA_real_)
}

delta_from_components <- function(atom, mu0, mu1, pi0, pi1) {
  combined1 <- pi1 * mu1 + (1 - pi1) * atom
  combined0 <- pi0 * mu0 + (1 - pi0) * atom
  as.numeric(combined1 - combined0)
}

delta_welch_interval <- function(y, r, conf.level) {
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

  if(identical(object$method, "Semi-empirical Likelihood Ratio Test")) {
    theta_fit <- el_mean_diff_theta(object$muDelta, y1, y0, tol = tol)
    if(isTRUE(theta_fit$feasible) && is.finite(theta_fit$theta)) {
      mu1 <- theta_fit$theta
      mu0 <- theta_fit$theta - object$muDelta
      return(delta_from_components(object$atom, mu0, mu1, pi0, pi1))
    }
  }

  delta_from_components(object$atom, mean(y0), mean(y1), pi0, pi1)
}

prepareParametricJointReference <- function(data, atom) {
  observed_data <- droplevels(data[data$A == 1, c("Y", "R"), drop = FALSE])
  glm_alt <- suppressWarnings(tryCatch(
    stats::glm(A ~ R, family = stats::binomial(), data = data),
    error = function(e) NULL
  ))
  lm_alt <- suppressWarnings(tryCatch(
    stats::lm(Y ~ R, data = observed_data),
    error = function(e) NULL
  ))

  list(
    data = data,
    observed_data = observed_data,
    glm_alt = glm_alt,
    lm_alt = lm_alt,
    logitReference = logit.prepare(data),
    ll_glm_alt = parametric_loglik_value(glm_alt),
    ll_lm_alt = parametric_loglik_value(lm_alt, reml = FALSE),
    observed_mean = mean(observed_data$Y),
    observed_treatment_rate = mean(observed_data$R),
    muDelta_hat = parametric_term_estimate(lm_alt, "R"),
    atom = atom
  )
}

parametricJointCandidate <- function(parametricReference, muDelta, logORdelta, tol = 1e-12) {
  glm_constrained <- suppressWarnings(tryCatch(
    stats::glm(
      A ~ 1,
      family = stats::binomial(),
      data = parametricReference$data,
      offset = logORdelta * parametricReference$data$R
    ),
    error = function(e) NULL
  ))

  lm_constrained <- suppressWarnings(tryCatch(
    stats::lm(
      Y ~ 1,
      data = parametricReference$observed_data,
      offset = muDelta * parametricReference$observed_data$R
    ),
    error = function(e) NULL
  ))

  ll_glm_constrained <- parametric_loglik_value(glm_constrained)
  ll_lm_constrained <- parametric_loglik_value(lm_constrained, reml = FALSE)

  W_A <- if(is.finite(ll_glm_constrained) && is.finite(parametricReference$ll_glm_alt)) {
    parametric_clamp_statistic(2 * (parametricReference$ll_glm_alt - ll_glm_constrained), tol = tol)
  } else {
    Inf
  }

  W_Y <- if(is.finite(ll_lm_constrained) && is.finite(parametricReference$ll_lm_alt)) {
    parametric_clamp_statistic(2 * (parametricReference$ll_lm_alt - ll_lm_constrained), tol = tol)
  } else {
    Inf
  }

  beta0 <- parametric_term_estimate(glm_constrained, "(Intercept)")
  mu0 <- parametric_term_estimate(lm_constrained, "(Intercept)")
  p0 <- if(is.finite(beta0)) stats::plogis(beta0) else NA_real_
  p1 <- if(is.finite(beta0)) stats::plogis(beta0 + logORdelta) else NA_real_
  mu1 <- if(is.finite(mu0)) mu0 + muDelta else NA_real_
  Delta <- if(all(is.finite(c(mu0, mu1, p0, p1)))) {
    delta_from_components(parametricReference$atom, mu0, mu1, p0, p1)
  } else {
    NA_real_
  }

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
    Delta = Delta,
    glm = glm_constrained,
    lm = lm_constrained
  )
}

prepareSPLRTJointReference <- function(data, atom) {
  list(
    yAlive0 = data$Y[data$R == 0 & data$A == 1],
    yAlive1 = data$Y[data$R == 1 & data$A == 1],
    logitReference = logit.prepare(data),
    atom = atom
  )
}

splrtContinuousCandidate <- function(splrtReference, muDelta, p0, p1, tol = 1e-8) {
  theta_fit <- el_mean_diff_theta(muDelta, splrtReference$yAlive1, splrtReference$yAlive0, tol = tol)

  if(!isTRUE(theta_fit$feasible) || !is.finite(theta_fit$theta)) {
    return(list(
      feasible = FALSE,
      statistic = Inf,
      theta = NA_real_,
      mu0 = NA_real_,
      mu1 = NA_real_,
      Delta = NA_real_
    ))
  }

  mu1 <- theta_fit$theta
  mu0 <- theta_fit$theta - muDelta
  Delta <- delta_from_components(splrtReference$atom, mu0, mu1, p0, p1)

  list(
    feasible = TRUE,
    statistic = as.numeric(theta_fit$statistic),
    theta = theta_fit$theta,
    mu0 = mu0,
    mu1 = mu1,
    Delta = Delta
  )
}

splrtJointCandidate <- function(splrtReference, muDelta, logORdelta, tol = 1e-8) {
  logit_fit <- logit_profile_fit.prepared(splrtReference$logitReference, logORdelta)
  continuous <- splrtContinuousCandidate(
    splrtReference,
    muDelta,
    p0 = logit_fit$p0,
    p1 = logit_fit$p1,
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
    mu1 = continuous$mu1,
    Delta = continuous$Delta
  )
}

delta_projection_candidate_factory <- function(object) {
  if(identical(object$method, "Parametric Likelihood Ratio Test")) {
    reference <- prepareParametricJointReference(object$data, atom = object$atom)
    return(function(muDelta, logORdelta) {
      parametricJointCandidate(reference, muDelta, logORdelta)
    })
  }

  if(identical(object$method, "Semi-empirical Likelihood Ratio Test")) {
    reference <- prepareSPLRTJointReference(object$data, atom = object$atom)
    return(function(muDelta, logORdelta) {
      splrtJointCandidate(reference, muDelta, logORdelta)
    })
  }

  stop("Projected Delta intervals are only implemented for unadjusted LRT and SPLRT fits.")
}

delta_projection_penalty_scale <- function(object) {
  y <- object$data$Y[is.finite(object$data$Y)]
  y_scale <- if(length(y) > 0) diff(range(y)) else NA_real_
  if(!is.finite(y_scale) || y_scale <= 0) {
    y_scale <- max(abs(y), na.rm = TRUE)
  }
  if(!is.finite(y_scale) || y_scale <= 0) {
    y_scale <- 1
  }

  1e5 * y_scale
}

delta_projection_center <- function(object, bounds) {
  log_alpha <- suppressWarnings(log(as.numeric(object$alphaDelta)))

  c(
    if(is.finite(object$muDelta)) object$muDelta else mean(bounds$muDelta),
    if(is.finite(log_alpha)) log_alpha else mean(bounds$logORdelta)
  )
}

delta_projection_start_points <- function(object, bounds) {
  center <- delta_projection_center(object, bounds)

  starts <- rbind(
    center,
    c(bounds$muDelta[1], bounds$logORdelta[1]),
    c(bounds$muDelta[1], bounds$logORdelta[2]),
    c(bounds$muDelta[2], bounds$logORdelta[1]),
    c(bounds$muDelta[2], bounds$logORdelta[2]),
    c(mean(bounds$muDelta), mean(bounds$logORdelta))
  )

  starts[, 1] <- pmin(pmax(starts[, 1], bounds$muDelta[1]), bounds$muDelta[2])
  starts[, 2] <- pmin(pmax(starts[, 2], bounds$logORdelta[1]), bounds$logORdelta[2])
  starts[!duplicated(signif(starts, digits = 10)), , drop = FALSE]
}

delta_projection_objective <- function(par, candidate_fun, threshold, direction, penalty_scale) {
  candidate <- candidate_fun(par[1], par[2])

  if(!is.finite(candidate$Delta) || !is.finite(candidate$statistic)) {
    return(penalty_scale * 1e6)
  }

  penalty <- max(candidate$statistic - threshold, 0)
  target <- if(identical(direction, "lower")) candidate$Delta else -candidate$Delta
  as.numeric(target + penalty_scale * penalty ^ 2)
}

delta_projection_is_feasible <- function(candidate, threshold, tol = 1e-8) {
  is.finite(candidate$Delta) &&
    is.finite(candidate$statistic) &&
    candidate$statistic <= threshold + tol
}

delta_projection_boundary_candidate <- function(candidate_fun, from, to, threshold,
                                                tol = 1e-8, max_iter = 60) {
  from_candidate <- candidate_fun(from[1], from[2])
  if(!delta_projection_is_feasible(from_candidate, threshold, tol = tol)) {
    return(list(par = from, candidate = from_candidate))
  }

  to_candidate <- candidate_fun(to[1], to[2])
  if(delta_projection_is_feasible(to_candidate, threshold, tol = tol)) {
    return(list(par = to, candidate = to_candidate))
  }

  lower_t <- 0
  upper_t <- 1
  best_par <- from
  best_candidate <- from_candidate

  for(i in seq_len(max_iter)) {
    mid_t <- (lower_t + upper_t) / 2
    mid_par <- from + mid_t * (to - from)
    mid_candidate <- candidate_fun(mid_par[1], mid_par[2])

    if(delta_projection_is_feasible(mid_candidate, threshold, tol = tol)) {
      lower_t <- mid_t
      best_par <- mid_par
      best_candidate <- mid_candidate
    } else {
      upper_t <- mid_t
    }

    if(abs(upper_t - lower_t) <= tol) {
      break
    }
  }

  list(par = best_par, candidate = best_candidate)
}

delta_projection_better <- function(candidate, best, direction, tol = 1e-8) {
  if(is.null(best)) {
    return(TRUE)
  }

  if(identical(direction, "lower")) {
    if(candidate$candidate$Delta < best$candidate$Delta - tol) {
      return(TRUE)
    }
    if(abs(candidate$candidate$Delta - best$candidate$Delta) <= tol &&
       candidate$candidate$statistic < best$candidate$statistic - tol) {
      return(TRUE)
    }
    return(FALSE)
  }

  if(candidate$candidate$Delta > best$candidate$Delta + tol) {
    return(TRUE)
  }
  if(abs(candidate$candidate$Delta - best$candidate$Delta) <= tol &&
     candidate$candidate$statistic < best$candidate$statistic - tol) {
    return(TRUE)
  }

  FALSE
}

delta_projection_in_box <- function(object, candidate_fun, bounds, threshold,
                                    direction = c("lower", "upper"), tol = 1e-8) {
  direction <- match.arg(direction)
  starts <- delta_projection_start_points(object, bounds)
  center <- delta_projection_center(object, bounds)
  lower <- c(bounds$muDelta[1], bounds$logORdelta[1])
  upper <- c(bounds$muDelta[2], bounds$logORdelta[2])
  penalty_scale <- delta_projection_penalty_scale(object)
  best <- NULL

  consider_point <- function(par) {
    boundary <- delta_projection_boundary_candidate(
      candidate_fun,
      from = center,
      to = par,
      threshold = threshold,
      tol = tol
    )
    if(delta_projection_is_feasible(boundary$candidate, threshold, tol = tol) &&
       delta_projection_better(boundary, best, direction, tol = tol)) {
      best <<- boundary
    }
  }

  for(i in seq_len(nrow(starts))) {
    start <- as.numeric(starts[i, ])
    consider_point(start)

    opt <- suppressWarnings(tryCatch(
      stats::nlminb(
        start = start,
        objective = delta_projection_objective,
        candidate_fun = candidate_fun,
        threshold = threshold,
        direction = direction,
        penalty_scale = penalty_scale,
        lower = lower,
        upper = upper,
        control = list(eval.max = 120, iter.max = 120, rel.tol = 1e-9, abs.tol = 1e-9)
      ),
      error = function(e) NULL
    ))

    if(!is.null(opt) && is.numeric(opt$par) && length(opt$par) == 2) {
      consider_point(as.numeric(opt$par))
    }
  }

  best
}

delta_projection_near_boundary <- function(solution, bounds, frac = 0.05, tol = 1e-6) {
  if(is.null(solution) || is.null(solution$par)) {
    return(TRUE)
  }

  mu_tol <- max(tol, frac * diff(bounds$muDelta))
  log_tol <- max(tol, frac * diff(bounds$logORdelta))

  solution$par[1] <= bounds$muDelta[1] + mu_tol ||
    solution$par[1] >= bounds$muDelta[2] - mu_tol ||
    solution$par[2] <= bounds$logORdelta[1] + log_tol ||
    solution$par[2] >= bounds$logORdelta[2] - log_tol
}

delta_projection_expand_bounds <- function(bounds, center, expansion = 1.5, additive = c(0.25, 0.25)) {
  mu_width <- diff(bounds$muDelta) / 2
  log_width <- diff(bounds$logORdelta) / 2

  list(
    muDelta = center[1] + c(-1, 1) * (mu_width * expansion + additive[1]),
    logORdelta = center[2] + c(-1, 1) * (log_width * expansion + additive[2])
  )
}

delta_projected_interval.optimize <- function(object, conf.level = object$conf.level,
                                              offset = NULL, expansion = 1.5,
                                              max_iter = 4, tol = 1e-8) {
  conf.level <- validateConfidenceLevel(conf.level)
  threshold <- stats::qchisq(conf.level, 2)
  candidate_fun <- delta_projection_candidate_factory(object)
  bounds <- jointContrastDefaultBounds(object, offset = offset)
  center <- delta_projection_center(object, bounds)
  previous_interval <- NULL
  best_lower <- NULL
  best_upper <- NULL

  for(iter in seq_len(max_iter)) {
    lower_solution <- delta_projection_in_box(
      object,
      candidate_fun,
      bounds,
      threshold,
      direction = "lower",
      tol = tol
    )
    upper_solution <- delta_projection_in_box(
      object,
      candidate_fun,
      bounds,
      threshold,
      direction = "upper",
      tol = tol
    )

    if(!is.null(lower_solution) && !is.null(upper_solution)) {
      current_interval <- c(lower_solution$candidate$Delta, upper_solution$candidate$Delta)
      if(!is.null(previous_interval)) {
        stable <- max(abs(current_interval - previous_interval)) <= max(1e-4, 1e-3 * max(1, abs(current_interval)))
        interior <- !delta_projection_near_boundary(lower_solution, bounds) &&
          !delta_projection_near_boundary(upper_solution, bounds)

        if(stable && interior) {
          best_lower <- lower_solution
          best_upper <- upper_solution
          break
        }
      }

      previous_interval <- current_interval
      best_lower <- lower_solution
      best_upper <- upper_solution
    }

    bounds <- delta_projection_expand_bounds(bounds, center, expansion = expansion)
  }

  if(is.null(best_lower) || is.null(best_upper)) {
    return(delta_na_interval())
  }

  as.numeric(c(best_lower$candidate$Delta, best_upper$candidate$Delta))
}

parametricContinuousCandidate <- function(parametricReference, muDelta, tol = 1e-12) {
  lm_constrained <- suppressWarnings(tryCatch(
    stats::lm(
      Y ~ 1,
      data = parametricReference$observed_data,
      offset = muDelta * parametricReference$observed_data$R
    ),
    error = function(e) NULL
  ))

  ll_lm_constrained <- parametric_loglik_value(lm_constrained, reml = FALSE)
  W_Y <- if(is.finite(ll_lm_constrained) && is.finite(parametricReference$ll_lm_alt)) {
    parametric_clamp_statistic(2 * (parametricReference$ll_lm_alt - ll_lm_constrained), tol = tol)
  } else {
    Inf
  }

  mu0 <- parametric_term_estimate(lm_constrained, "(Intercept)")
  mu1 <- if(is.finite(mu0)) mu0 + muDelta else NA_real_

  list(
    statistic = W_Y,
    mu0 = mu0,
    mu1 = mu1,
    lm = lm_constrained
  )
}

delta_profile_value_scale <- function(object) {
  delta_welch_ci <- delta_welch_interval(object$data$Y, object$data$R, object$conf.level)
  if(length(delta_welch_ci) >= 2 && all(is.finite(delta_welch_ci[1:2]))) {
    width <- diff(delta_welch_ci[1:2])
    if(is.finite(width) && width > 0) {
      return(as.numeric(width))
    }
  }

  y <- object$data$Y[is.finite(object$data$Y)]
  if(length(y) > 0) {
    scale <- diff(range(y))
    if(is.finite(scale) && scale > 0) {
      return(as.numeric(scale))
    }
  }

  1
}

delta_profile_initial_step <- function(object) {
  scale <- delta_profile_value_scale(object)
  max(scale / 4, 0.25)
}

delta_profile_target_tolerance <- function(target, scale = 1) {
  max(1e-8, sqrt(.Machine$double.eps) * max(1, abs(target), scale))
}

delta_profile_interval_width <- function(object) {
  offsets <- jointContrastDefaultOffsets(object)
  max(as.numeric(offsets[["log_or_delta"]]), 0.5)
}

delta_profile_cache_key <- function(x) {
  formatC(as.numeric(x), digits = 16, format = "fg", flag = "#")
}

delta_profile_lrt_candidate_factory <- function(object, tol = 1e-8) {
  reference <- prepareParametricJointReference(object$data, atom = object$atom)
  ybar <- reference$observed_mean
  rbar <- reference$observed_treatment_rate
  atom <- object$atom
  logit_cache <- new.env(parent = emptyenv())

  get_logit_fit <- function(logORdelta) {
    key <- delta_profile_cache_key(logORdelta)
    if(exists(key, envir = logit_cache, inherits = FALSE)) {
      return(get(key, envir = logit_cache, inherits = FALSE))
    }

    fit <- logit_profile_fit.prepared(reference$logitReference, logORdelta)
    assign(key, fit, envir = logit_cache)
    fit
  }

  function(targetDelta, logORdelta) {
    logit_fit <- get_logit_fit(logORdelta)
    coefficient <- logit_fit$p1 * (1 - rbar) + logit_fit$p0 * rbar
    if(!is.finite(coefficient) || coefficient <= tol) {
      return(list(statistic = Inf, Delta = NA_real_, muDelta = NA_real_, logORdelta = logORdelta))
    }

    base <- (logit_fit$p1 - logit_fit$p0) * (ybar - atom)
    muDelta <- (targetDelta - base) / coefficient
    continuous <- parametricContinuousCandidate(reference, muDelta, tol = tol)
    total <- as.numeric(logit_fit$statistic + continuous$statistic)
    if(is.finite(total)) {
      total <- parametric_clamp_statistic(total, tol = tol)
    }

    Delta <- if(all(is.finite(c(continuous$mu0, continuous$mu1, logit_fit$p0, logit_fit$p1)))) {
      delta_from_components(atom, continuous$mu0, continuous$mu1, logit_fit$p0, logit_fit$p1)
    } else {
      NA_real_
    }

    list(
      statistic = total,
      Delta = Delta,
      muDelta = muDelta,
      logORdelta = logORdelta,
      p0 = logit_fit$p0,
      p1 = logit_fit$p1,
      mu0 = continuous$mu0,
      mu1 = continuous$mu1
    )
  }
}

delta_profile_splrt_mu_solver <- function(reference, p0, p1, targetDelta, tol = 1e-8) {
  mu_range <- el_mean_diff_delta_range(reference$yAlive1, reference$yAlive0)
  local_cache <- new.env(parent = emptyenv())
  delta_tol <- delta_profile_target_tolerance(targetDelta)

  evaluate_mu <- function(muDelta) {
    key <- delta_profile_cache_key(muDelta)
    if(exists(key, envir = local_cache, inherits = FALSE)) {
      return(get(key, envir = local_cache, inherits = FALSE))
    }

    candidate <- splrtContinuousCandidate(reference, muDelta, p0 = p0, p1 = p1, tol = tol)
    assign(key, candidate, envir = local_cache)
    candidate
  }

  lower_candidate <- evaluate_mu(mu_range[1])
  upper_candidate <- evaluate_mu(mu_range[2])

  lower_diff <- lower_candidate$Delta - targetDelta
  upper_diff <- upper_candidate$Delta - targetDelta

  if(is.finite(lower_diff) && abs(lower_diff) <= delta_tol) {
    lower_candidate$muDelta <- mu_range[1]
    return(lower_candidate)
  }
  if(is.finite(upper_diff) && abs(upper_diff) <= delta_tol) {
    upper_candidate$muDelta <- mu_range[2]
    return(upper_candidate)
  }

  if(is.finite(lower_diff) && is.finite(upper_diff) && sign(lower_diff) != sign(upper_diff)) {
    root <- stats::uniroot(
      function(muDelta) evaluate_mu(muDelta)$Delta - targetDelta,
      interval = mu_range,
      tol = tol
    )$root
    candidate <- evaluate_mu(root)
    candidate$muDelta <- root
    return(candidate)
  }

  nearest <- suppressWarnings(tryCatch(
    stats::optimize(
      function(muDelta) {
        candidate <- evaluate_mu(muDelta)
        if(!is.finite(candidate$Delta)) {
          return(1e50)
        }
        (candidate$Delta - targetDelta) ^ 2
      },
      interval = mu_range,
      tol = tol
    ),
    error = function(e) NULL
  ))

  if(is.null(nearest)) {
    return(list(feasible = FALSE, statistic = Inf, Delta = NA_real_, muDelta = NA_real_))
  }

  candidate <- evaluate_mu(nearest$minimum)
  candidate$muDelta <- nearest$minimum
  if(is.finite(candidate$Delta) && abs(candidate$Delta - targetDelta) <= delta_tol) {
    return(candidate)
  }

  list(feasible = FALSE, statistic = Inf, Delta = NA_real_, muDelta = NA_real_)
}

delta_profile_splrt_candidate_factory <- function(object, tol = 1e-8) {
  reference <- prepareSPLRTJointReference(object$data, atom = object$atom)
  logit_cache <- new.env(parent = emptyenv())

  get_logit_fit <- function(logORdelta) {
    key <- delta_profile_cache_key(logORdelta)
    if(exists(key, envir = logit_cache, inherits = FALSE)) {
      return(get(key, envir = logit_cache, inherits = FALSE))
    }

    fit <- logit_profile_fit.prepared(reference$logitReference, logORdelta)
    assign(key, fit, envir = logit_cache)
    fit
  }

  function(targetDelta, logORdelta) {
    logit_fit <- get_logit_fit(logORdelta)
    continuous <- delta_profile_splrt_mu_solver(
      reference,
      p0 = logit_fit$p0,
      p1 = logit_fit$p1,
      targetDelta = targetDelta,
      tol = tol
    )

    if(!isTRUE(continuous$feasible) || !is.finite(continuous$Delta)) {
      return(list(statistic = Inf, Delta = NA_real_, muDelta = NA_real_, logORdelta = logORdelta))
    }

    total <- as.numeric(logit_fit$statistic + continuous$statistic)
    if(is.finite(total) && total < tol) {
      total <- 0
    }

    list(
      statistic = total,
      Delta = continuous$Delta,
      muDelta = continuous$muDelta,
      logORdelta = logORdelta,
      p0 = logit_fit$p0,
      p1 = logit_fit$p1,
      mu0 = continuous$mu0,
      mu1 = continuous$mu1
    )
  }
}

delta_profile_candidate_factory <- function(object, tol = 1e-8) {
  if(identical(object$method, "Parametric Likelihood Ratio Test")) {
    return(delta_profile_lrt_candidate_factory(object, tol = tol))
  }

  if(identical(object$method, "Semi-empirical Likelihood Ratio Test")) {
    return(delta_profile_splrt_candidate_factory(object, tol = tol))
  }

  stop("Profile Delta intervals are only implemented for unadjusted LRT and SPLRT fits.")
}

delta_profile_objective_value <- function(candidate) {
  if(is.null(candidate) || !is.finite(candidate$statistic)) {
    return(1e50)
  }

  as.numeric(candidate$statistic)
}

delta_profile_logor_near_boundary <- function(logORdelta, lower, upper, frac = 0.05, tol = 1e-6) {
  span <- upper - lower
  boundary_tol <- max(tol, frac * span)

  logORdelta <= lower + boundary_tol || logORdelta >= upper - boundary_tol
}

delta_profile_optimize_logor <- function(candidate_fun, center, width,
                                         expansion = 1.5, max_iter = 6,
                                         tol = 1e-6) {
  width <- max(width, 0.5)
  best <- NULL

  consider <- function(candidate) {
    if(is.null(candidate)) {
      return()
    }

    if(is.null(best) || delta_profile_objective_value(candidate) < delta_profile_objective_value(best)) {
      best <<- candidate
    }
  }

  eval_cache <- new.env(parent = emptyenv())
  evaluate <- function(logORdelta) {
    key <- delta_profile_cache_key(logORdelta)
    if(exists(key, envir = eval_cache, inherits = FALSE)) {
      return(get(key, envir = eval_cache, inherits = FALSE))
    }

    candidate <- candidate_fun(logORdelta)
    assign(key, candidate, envir = eval_cache)
    candidate
  }

  for(iter in seq_len(max_iter)) {
    lower <- center - width
    upper <- center + width

    consider(evaluate(lower))
    consider(evaluate(center))
    consider(evaluate(upper))

    objective <- function(logORdelta) {
      delta_profile_objective_value(evaluate(logORdelta))
    }

    opt <- suppressWarnings(tryCatch(
      stats::optimize(objective, interval = c(lower, upper), tol = tol),
      error = function(e) NULL
    ))

    if(!is.null(opt)) {
      consider(evaluate(opt$minimum))
    }

    if(!is.null(best) && is.finite(best$logORdelta)) {
      lower_value <- delta_profile_objective_value(evaluate(lower))
      upper_value <- delta_profile_objective_value(evaluate(upper))
      best_value <- delta_profile_objective_value(best)
      interior <- !delta_profile_logor_near_boundary(best$logORdelta, lower, upper)
      separated <- lower_value > best_value + max(1e-4, 1e-4 * max(1, abs(best_value))) &&
        upper_value > best_value + max(1e-4, 1e-4 * max(1, abs(best_value)))

      if(interior && separated) {
        break
      }
    }

    width <- width * expansion + 0.25
  }

  best
}

delta_profile_factory <- function(object, tol = 1e-8) {
  candidate_builder <- delta_profile_candidate_factory(object, tol = tol)
  center <- suppressWarnings(log(as.numeric(object$alphaDelta)))
  if(!is.finite(center)) {
    center <- 0
  }
  width <- delta_profile_interval_width(object)
  scale <- delta_profile_value_scale(object)
  delta_cache <- new.env(parent = emptyenv())

  function(targetDelta) {
    key <- delta_profile_cache_key(targetDelta)
    if(exists(key, envir = delta_cache, inherits = FALSE)) {
      return(get(key, envir = delta_cache, inherits = FALSE))
    }

    if(abs(targetDelta - object$Delta) <= delta_profile_target_tolerance(object$Delta, scale = scale)) {
      fitted_log_or <- suppressWarnings(log(as.numeric(object$alphaDelta)))
      out <- list(
        statistic = 0,
        Delta = object$Delta,
        muDelta = object$muDelta,
        logORdelta = if(is.finite(fitted_log_or)) fitted_log_or else center
      )
      assign(key, out, envir = delta_cache)
      return(out)
    }

    best <- delta_profile_optimize_logor(
      function(logORdelta) candidate_builder(targetDelta, logORdelta),
      center = center,
      width = width,
      tol = tol
    )

    if(is.null(best) || !is.finite(best$Delta) ||
       abs(best$Delta - targetDelta) > delta_profile_target_tolerance(targetDelta, scale = scale)) {
      best <- list(statistic = Inf, Delta = NA_real_, muDelta = NA_real_, logORdelta = NA_real_)
    }

    assign(key, best, envir = delta_cache)
    best
  }
}

delta_profile_find_bound <- function(profile_fun, estimate, crit, direction = c("lower", "upper"),
                                     initial_step, tol = 1e-6, max_expand = 24) {
  direction <- match.arg(direction)
  sign_dir <- if(identical(direction, "lower")) -1 else 1

  inside <- estimate
  step <- initial_step
  outside <- NULL

  for(i in seq_len(max_expand)) {
    candidate <- estimate + sign_dir * step
    statistic <- profile_fun(candidate)$statistic
    if(!is.finite(statistic) || statistic > crit) {
      outside <- candidate
      break
    }
    inside <- candidate
    step <- step * 2
  }

  if(is.null(outside)) {
    return(NA_real_)
  }

  while(abs(outside - inside) > tol * max(1, abs(estimate), abs(outside), abs(inside))) {
    midpoint <- (inside + outside) / 2
    statistic <- profile_fun(midpoint)$statistic
    if(!is.finite(statistic) || statistic > crit) {
      outside <- midpoint
    } else {
      inside <- midpoint
    }
  }

  (inside + outside) / 2
}

delta_profile_interval.optimize <- function(object, conf.level = object$conf.level,
                                            tol = 1e-6) {
  if(!isTRUE(object$success) || !is.null(object$adjust) || is.null(object$atom)) {
    return(delta_na_interval())
  }

  profile_fun <- delta_profile_factory(object, tol = tol)
  crit <- stats::qchisq(conf.level, 1)
  step <- delta_profile_initial_step(object)

  lower <- delta_profile_find_bound(
    profile_fun,
    estimate = object$Delta,
    crit = crit,
    direction = "lower",
    initial_step = step,
    tol = tol
  )
  upper <- delta_profile_find_bound(
    profile_fun,
    estimate = object$Delta,
    crit = crit,
    direction = "upper",
    initial_step = step,
    tol = tol
  )

  if(!all(is.finite(c(lower, upper)))) {
    return(delta_na_interval())
  }

  as.numeric(c(lower, upper))
}

delta_boundary_touched <- function(mask) {
  if(!any(mask)) {
    return(FALSE)
  }

  any(mask[1, ]) || any(mask[nrow(mask), ]) || any(mask[, 1]) || any(mask[, ncol(mask)])
}

delta_surface_for_inference <- function(object, conf.level = object$conf.level,
                                        resolution = 31, offset = NULL,
                                        expansion = 1.5, max_iter = 6) {
  current_offset <- normalizeJointContrastOffsets(object, offset)
  surface_data <- NULL
  threshold <- stats::qchisq(conf.level, 2)

  for(iter in seq_len(max_iter)) {
    surface_data <- jointContrastSurfaceData(
      object,
      conf.level = conf.level,
      plot = FALSE,
      offset = current_offset,
      resolution = resolution,
      include_delta = TRUE
    )
    accepted <- is.finite(surface_data$surface) &
      is.finite(surface_data$deltaSurface) &
      surface_data$surface <= threshold

    if(any(accepted) && !delta_boundary_touched(accepted)) {
      break
    }

    current_offset <- current_offset * expansion + 0.25
  }

  surface_data
}

delta_interval_from_surface <- function(surface_data, threshold) {
  accepted <- is.finite(surface_data$surface) &
    is.finite(surface_data$deltaSurface) &
    surface_data$surface <= threshold

  if(!any(accepted)) {
    return(delta_na_interval())
  }

  as.numeric(range(surface_data$deltaSurface[accepted]))
}

validateDeltaIntervalAlgorithm <- function(algorithm) {
  algorithm <- match.arg(algorithm, c("grid", "optimize"))
  algorithm
}

delta_projected_interval.grid <- function(object, conf.level = object$conf.level,
                                          offset = NULL, resolution = 35,
                                          expansion = 1.5, max_iter = 6) {
  surface_data <- delta_surface_for_inference(
    object,
    conf.level = conf.level,
    resolution = resolution,
    offset = offset,
    expansion = expansion,
    max_iter = max_iter
  )

  delta_interval_from_surface(surface_data, stats::qchisq(conf.level, 2))
}

delta_projected_interval <- function(object, conf.level = object$conf.level,
                                     offset = NULL, resolution = 35,
                                     algorithm = c("grid", "optimize"),
                                     expansion = 1.5, max_iter = 6, tol = 1e-8) {
  algorithm <- validateDeltaIntervalAlgorithm(algorithm)

  if(identical(algorithm, "grid")) {
    return(delta_projected_interval.grid(
      object,
      conf.level = conf.level,
      offset = offset,
      resolution = resolution,
      expansion = expansion,
      max_iter = max_iter
    ))
  }

  delta_projected_interval.optimize(
    object,
    conf.level = conf.level,
    offset = offset,
    expansion = expansion,
    max_iter = max_iter,
    tol = tol
  )
}

delta_profile_interval.grid <- function(object, conf.level = object$conf.level,
                                        resolution = 31, offset = NULL,
                                        expansion = 1.5, max_iter = 6) {
  surface_data <- delta_surface_for_inference(
    object,
    conf.level = conf.level,
    resolution = resolution,
    offset = offset,
    expansion = expansion,
    max_iter = max_iter
  )

  delta_interval_from_surface(surface_data, stats::qchisq(conf.level, 1))
}

delta_profile_interval <- function(object, conf.level = object$conf.level,
                                   resolution = 31, offset = NULL,
                                   algorithm = c("grid", "optimize"),
                                   expansion = 1.5, max_iter = 6, tol = 1e-6) {
  algorithm <- validateDeltaIntervalAlgorithm(algorithm)

  if(identical(algorithm, "grid")) {
    return(delta_profile_interval.grid(
      object,
      conf.level = conf.level,
      resolution = resolution,
      offset = offset,
      expansion = expansion,
      max_iter = max_iter
    ))
  }

  delta_profile_interval.optimize(
    object,
    conf.level = conf.level,
    tol = tol
  )
}

augmentDeltaInference <- function(object) {
  if(!isTRUE(object$success) || !is.null(object$adjust) || is.null(object$atom)) {
    object$delta <- if(is.null(object$delta)) NA_real_ else object$delta
    return(syncTruncComp2Aliases(object))
  }

  object$delta <- delta_unadjusted_point_estimate(object)
  syncTruncComp2Aliases(object)
}
