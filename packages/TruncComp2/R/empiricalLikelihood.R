el_validate_mean_diff_inputs <- function(x, y, conf.level = NULL) {
  if(!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric.")
  }

  if(any(!is.finite(x)) || any(!is.finite(y))) {
    stop("x and y must contain only finite values.")
  }

  if(length(x) < 2 || length(y) < 2) {
    stop("x and y must each contain at least two observations.")
  }

  if(!is.null(conf.level) && !(length(conf.level) == 1 &&
                               is.finite(conf.level) &&
                               conf.level > 0 &&
                               conf.level < 1)) {
    stop("conf.level must be a single number strictly between 0 and 1.")
  }
}

el_root_epsilon <- function(lower, upper) {
  span <- abs(upper - lower)
  scale <- max(1, abs(lower), abs(upper), span)
  min(sqrt(.Machine$double.eps) * scale, span / 4)
}

el_one_sample_lambda <- function(residuals, tol = 1e-10) {
  residuals <- as.numeric(residuals)

  if(all(abs(residuals) <= tol)) {
    return(list(feasible = TRUE,
                lambda = 0,
                statistic = 0,
                gradient = 0))
  }

  r_min <- min(residuals)
  r_max <- max(residuals)
  if(!(r_min < -tol && r_max > tol)) {
    return(list(feasible = FALSE,
                lambda = NA_real_,
                statistic = Inf,
                gradient = NA_real_))
  }

  lower <- -1 / r_max
  upper <- -1 / r_min
  eps <- el_root_epsilon(lower, upper)
  objective <- function(lambda) {
    sum(residuals / (1 + lambda * residuals))
  }

  lambda <- stats::uniroot(objective, c(lower + eps, upper - eps), tol = tol)$root
  denom <- 1 + lambda * residuals
  statistic <- 2 * sum(log1p(lambda * residuals))
  if(abs(statistic) < tol) {
    statistic <- 0
  }

  list(feasible = TRUE,
       lambda = lambda,
       statistic = statistic,
       gradient = -lambda * sum(1 / denom))
}

el_mean_diff_delta_range <- function(x, y) {
  c(min(x) - max(y), max(x) - min(y))
}

el_mean_diff_theta_range <- function(delta, x, y) {
  c(max(min(x), min(y) + delta),
    min(max(x), max(y) + delta))
}

el_mean_diff_theta_eval <- function(theta, delta, x, y) {
  fit_x <- el_one_sample_lambda(x - theta)
  fit_y <- el_one_sample_lambda(y - theta + delta)

  if(!fit_x$feasible || !fit_y$feasible) {
    return(list(feasible = FALSE,
                theta = theta,
                statistic = Inf,
                gradient = NA_real_))
  }

  list(feasible = TRUE,
       theta = theta,
       statistic = fit_x$statistic + fit_y$statistic,
       gradient = fit_x$gradient + fit_y$gradient)
}

el_mean_diff_theta <- function(delta, x, y, tol = 1e-8) {
  theta_range <- el_mean_diff_theta_range(delta, x, y)
  lower <- theta_range[1]
  upper <- theta_range[2]

  if(lower > upper + tol) {
    return(list(feasible = FALSE,
                theta = NA_real_,
                statistic = Inf,
                gradient = NA_real_))
  }

  if(abs(upper - lower) <= tol) {
    return(el_mean_diff_theta_eval(lower, delta, x, y))
  }

  theta_grad <- function(theta) {
    el_mean_diff_theta_eval(theta, delta, x, y)$gradient
  }

  eps <- el_root_epsilon(lower, upper)
  inner_lower <- lower + eps
  inner_upper <- upper - eps

  if(inner_lower >= inner_upper) {
    return(el_mean_diff_theta_eval((lower + upper) / 2, delta, x, y))
  }

  grad_lower <- suppressWarnings(theta_grad(inner_lower))
  grad_upper <- suppressWarnings(theta_grad(inner_upper))

  if(is.finite(grad_lower) && is.finite(grad_upper)) {
    if(abs(grad_lower) <= tol) {
      return(el_mean_diff_theta_eval(inner_lower, delta, x, y))
    }
    if(abs(grad_upper) <= tol) {
      return(el_mean_diff_theta_eval(inner_upper, delta, x, y))
    }
    if(sign(grad_lower) != sign(grad_upper)) {
      theta <- stats::uniroot(theta_grad,
                              c(inner_lower, inner_upper),
                              tol = tol)$root
      return(el_mean_diff_theta_eval(theta, delta, x, y))
    }
  }

  objective <- function(theta) {
    statistic <- el_mean_diff_theta_eval(theta, delta, x, y)$statistic
    if(is.finite(statistic)) statistic else 1e50
  }

  optimum <- stats::optimize(objective, c(inner_lower, inner_upper))
  el_mean_diff_theta_eval(optimum$minimum, delta, x, y)
}

el_mean_diff_statistic <- function(x, y, delta, tol = 1e-8) {
  el_validate_mean_diff_inputs(x, y)

  estimate <- mean(x) - mean(y)
  if(abs(delta - estimate) <= tol) {
    return(0)
  }

  delta_range <- el_mean_diff_delta_range(x, y)
  if(delta < delta_range[1] - tol || delta > delta_range[2] + tol) {
    return(Inf)
  }

  theta_fit <- el_mean_diff_theta(delta, x, y, tol = tol)
  as.numeric(theta_fit$statistic)
}

el_find_mean_diff_bound <- function(x, y, estimate, boundary, crit, side, tol = 1e-8) {
  if(abs(boundary - estimate) <= tol) {
    return(boundary)
  }

  objective <- function(delta) {
    el_mean_diff_statistic(x, y, delta, tol = tol) - crit
  }

  boundary_value <- suppressWarnings(objective(boundary))
  if(is.finite(boundary_value) && boundary_value <= 0) {
    return(boundary)
  }

  if(side == "lower") {
    outer <- boundary
    inner <- estimate
    for(i in 1:200) {
      mid <- (outer + inner) / 2
      value <- suppressWarnings(objective(mid))
      if(!is.finite(value) || value > 0) {
        outer <- mid
      } else {
        inner <- mid
      }
      if(abs(inner - outer) <= tol) {
        break
      }
    }
  } else {
    inner <- estimate
    outer <- boundary
    for(i in 1:200) {
      mid <- (inner + outer) / 2
      value <- suppressWarnings(objective(mid))
      if(!is.finite(value) || value > 0) {
        outer <- mid
      } else {
        inner <- mid
      }
      if(abs(inner - outer) <= tol) {
        break
      }
    }
  }

  (inner + outer) / 2
}

el_mean_diff_confint <- function(x, y, estimate, conf.level, tol = 1e-8) {
  delta_range <- el_mean_diff_delta_range(x, y)
  crit <- stats::qchisq(conf.level, 1)

  lower <- el_find_mean_diff_bound(x, y, estimate, delta_range[1], crit, "lower", tol = tol)
  upper <- el_find_mean_diff_bound(x, y, estimate, delta_range[2], crit, "upper", tol = tol)
  conf_int <- c(lower, upper)
  attr(conf_int, "conf.level") <- conf.level
  conf_int
}

el_mean_diff_fit <- function(x, y, mu = 0, conf.level = 0.95) {
  call <- match.call()
  el_validate_mean_diff_inputs(x, y, conf.level = conf.level)

  estimate <- mean(x) - mean(y)
  statistic <- el_mean_diff_statistic(x, y, mu)
  names(statistic) <- "-2 * LogLikelihood"
  p_value <- stats::pchisq(statistic, 1, lower.tail = FALSE)
  names(p_value) <- "-2 * LogLikelihood"

  out <- list(estimate = structure(estimate, names = "Mean difference"),
              conf.int = el_mean_diff_confint(x, y, estimate, conf.level),
              p.value = p_value,
              statistic = statistic,
              method = "Empirical likelihood mean difference test",
              null.value = mu,
              data.name = paste(deparse(call$x), " and ", deparse(call$y)))
  class(out) <- "htest"
  out
}

el_regression_design <- function(formula, data, term = "R") {
  model_frame <- stats::model.frame(formula, data = data, na.action = stats::na.fail)
  model_frame <- droplevels(model_frame)
  y <- stats::model.response(model_frame)
  X <- stats::model.matrix(stats::terms(model_frame), model_frame)
  term_index <- match(term, colnames(X))

  if(is.na(term_index)) {
    return(list(success = FALSE,
                error = paste0("The treatment term '", term, "' is not identifiable in the observed-outcome design.")))
  }

  qr_x <- qr(X)
  if(qr_x$rank < ncol(X)) {
    return(list(success = FALSE,
                error = "The observed-outcome design matrix is rank deficient."))
  }

  if(nrow(X) <= ncol(X)) {
    return(list(success = FALSE,
                error = "There are too few observed outcomes relative to the adjusted design."))
  }

  lm_fit <- stats::lm.fit(x = X, y = y)
  coefficients <- as.numeric(lm_fit$coefficients)
  names(coefficients) <- colnames(X)

  list(
    success = TRUE,
    formula = formula,
    y = as.numeric(y),
    X = unname(X),
    X_names = colnames(X),
    term = term,
    term_index = term_index,
    qr = qr_x,
    lm_fit = lm_fit,
    coefficients = coefficients,
    estimate = as.numeric(coefficients[[term]])
  )
}

el_regression_beta <- function(delta, nuisance, term_index, p) {
  beta <- numeric(p)
  nuisance_index <- seq_len(p) != term_index
  beta[term_index] <- delta
  beta[nuisance_index] <- nuisance
  beta
}

el_regression_restricted_ls <- function(design, delta) {
  nuisance_index <- seq_len(ncol(design$X)) != design$term_index
  Z <- design$X[, nuisance_index, drop = FALSE]
  response <- design$y - design$X[, design$term_index] * delta

  if(ncol(Z) == 0) {
    return(numeric(0))
  }

  restricted_fit <- stats::lm.fit(x = Z, y = response)
  coefficients <- as.numeric(restricted_fit$coefficients)

  if(any(!is.finite(coefficients))) {
    stop("Unable to construct restricted least-squares starting values for the empirical-likelihood profile.")
  }

  coefficients
}

el_regression_moments <- function(design, beta) {
  residuals <- as.numeric(design$y - design$X %*% beta)
  design$X * residuals
}

el_dual_state <- function(lambda, G) {
  t_val <- as.numeric(G %*% lambda)
  denom <- 1 + t_val
  feasible <- all(is.finite(denom)) && all(denom > 0)

  if(!feasible) {
    return(list(feasible = FALSE,
                objective = Inf,
                statistic = Inf,
                score = rep(NA_real_, ncol(G)),
                hessian = matrix(NA_real_, ncol(G), ncol(G)),
                denom = denom))
  }

  inv_denom <- 1 / denom
  score <- -colSums(G * inv_denom)
  hessian <- crossprod(G, G * (inv_denom ^ 2))
  objective <- -sum(log(denom))

  list(feasible = TRUE,
       objective = objective,
       statistic = 2 * sum(log(denom)),
       score = score,
       hessian = hessian,
       denom = denom)
}

el_dual_fit_newton <- function(G, tol = 1e-8, maxit = 100) {
  p <- ncol(G)
  lambda <- rep(0, p)
  state <- el_dual_state(lambda, G)

  if(max(abs(colSums(G))) <= tol) {
    return(list(success = TRUE,
                feasible = TRUE,
                lambda = lambda,
                statistic = 0,
                score = colSums(G),
                denom = rep(1, nrow(G)),
                iterations = 0,
                method = "newton"))
  }

  for(iter in seq_len(maxit)) {
    if(max(abs(state$score)) <= tol) {
      return(list(success = TRUE,
                  feasible = TRUE,
                  lambda = lambda,
                  statistic = max(0, state$statistic),
                  score = state$score,
                  denom = state$denom,
                  iterations = iter,
                  method = "newton"))
    }

    step_direction <- tryCatch(
      as.numeric(solve(state$hessian, state$score)),
      error = function(e) NULL
    )

    if(is.null(step_direction) || any(!is.finite(step_direction))) {
      break
    }

    step_size <- 1
    accepted <- FALSE
    for(backtrack in seq_len(60)) {
      candidate <- lambda - step_size * step_direction
      candidate_state <- el_dual_state(candidate, G)
      if(candidate_state$feasible && candidate_state$objective < state$objective) {
        lambda <- candidate
        state <- candidate_state
        accepted <- TRUE
        break
      }
      step_size <- step_size / 2
    }

    if(!accepted) {
      break
    }
  }

  list(success = FALSE,
       feasible = FALSE,
       lambda = lambda,
       statistic = Inf,
       score = state$score,
       denom = state$denom,
       iterations = maxit,
       method = "newton")
}

el_dual_fit_nlminb <- function(G, start = NULL, tol = 1e-8) {
  p <- ncol(G)
  if(is.null(start)) {
    start <- rep(0, p)
  }

  objective <- function(lambda) {
    el_dual_state(lambda, G)$objective
  }

  gradient <- function(lambda) {
    el_dual_state(lambda, G)$score
  }

  fit <- tryCatch(
    suppressWarnings(stats::nlminb(
      start = start,
      objective = objective,
      gradient = gradient,
      control = list(rel.tol = tol, x.tol = tol, eval.max = 500, iter.max = 500)
    )),
    error = function(e) NULL
  )

  if(is.null(fit)) {
    return(list(success = FALSE,
                feasible = FALSE,
                lambda = start,
                statistic = Inf,
                score = rep(NA_real_, p),
                denom = rep(NA_real_, nrow(G)),
                iterations = NA_integer_,
                method = "nlminb"))
  }

  state <- el_dual_state(fit$par, G)
  converged <- fit$convergence == 0 &&
    state$feasible &&
    max(abs(state$score)) <= max(tol * 10, 1e-7)

  list(success = converged,
       feasible = state$feasible,
       lambda = fit$par,
       statistic = if(converged) max(0, state$statistic) else Inf,
       score = state$score,
       denom = state$denom,
       iterations = fit$iterations,
       method = "nlminb")
}

el_dual_fit <- function(G, tol = 1e-8) {
  G <- as.matrix(G)

  if(nrow(G) < 1 || ncol(G) < 1) {
    stop("G must have positive dimensions.")
  }

  if(any(!is.finite(G))) {
    return(list(success = FALSE,
                feasible = FALSE,
                statistic = Inf))
  }

  fit <- el_dual_fit_newton(G, tol = tol)
  if(!isTRUE(fit$success)) {
    fit <- el_dual_fit_nlminb(G, start = fit$lambda, tol = tol)
  }

  if(!isTRUE(fit$success)) {
    return(list(success = FALSE,
                feasible = FALSE,
                statistic = Inf))
  }

  weights <- 1 / (nrow(G) * fit$denom)

  list(
    success = TRUE,
    feasible = TRUE,
    lambda = fit$lambda,
    statistic = max(0, fit$statistic),
    score = fit$score,
    weights = as.numeric(weights),
    weighted_moments = as.numeric(colSums(G * weights)),
    method = fit$method
  )
}

el_regression_profile_fit <- function(design, delta, tol = 1e-8) {
  if(!isTRUE(design$success)) {
    return(list(success = FALSE,
                feasible = FALSE,
                statistic = Inf,
                error = design$error))
  }

  p <- ncol(design$X)
  nuisance_index <- seq_len(p) != design$term_index
  nuisance_hat <- design$coefficients[nuisance_index]

  if(abs(delta - design$estimate) <= tol) {
    beta_hat <- as.numeric(design$coefficients)
    moments_hat <- el_regression_moments(design, beta_hat)
    return(list(
      success = TRUE,
      feasible = TRUE,
      beta = beta_hat,
      nuisance = as.numeric(nuisance_hat),
      statistic = 0,
      dual = list(
        success = TRUE,
        feasible = TRUE,
        lambda = rep(0, p),
        statistic = 0,
        score = colSums(moments_hat),
        weights = rep(1 / nrow(design$X), nrow(design$X)),
        weighted_moments = colSums(moments_hat) / nrow(design$X),
        method = "closed-form"
      )
    ))
  }

  objective <- function(nuisance) {
    beta <- el_regression_beta(delta, nuisance, design$term_index, p)
    moments <- el_regression_moments(design, beta)
    dual <- el_dual_fit(moments, tol = tol)

    if(!isTRUE(dual$success)) {
      return(Inf)
    }

    dual$statistic
  }

  nuisance_start <- tryCatch(
    el_regression_restricted_ls(design, delta),
    error = function(e) NULL
  )

  if(is.null(nuisance_start)) {
    return(list(success = FALSE,
                feasible = FALSE,
                statistic = Inf,
                error = "Unable to profile nuisance coefficients for the empirical-likelihood regression fit."))
  }

  if(length(nuisance_start) == 0) {
    beta <- el_regression_beta(delta, numeric(0), design$term_index, p)
    moments <- el_regression_moments(design, beta)
    dual <- el_dual_fit(moments, tol = tol)

    if(!isTRUE(dual$success)) {
      return(list(success = FALSE,
                  feasible = FALSE,
                  statistic = Inf,
                  error = "Unable to solve the empirical-likelihood dual problem for the constrained treatment effect."))
    }

    return(list(success = TRUE,
                feasible = TRUE,
                beta = beta,
                nuisance = numeric(0),
                statistic = dual$statistic,
                dual = dual))
  }

  optimizer <- tryCatch(
    suppressWarnings(stats::nlminb(
      start = nuisance_start,
      objective = objective,
      control = list(rel.tol = tol, x.tol = tol, eval.max = 500, iter.max = 500)
    )),
    error = function(e) NULL
  )

  if(is.null(optimizer) || optimizer$convergence != 0 || !is.finite(optimizer$objective)) {
    return(list(success = FALSE,
                feasible = FALSE,
                statistic = Inf,
                error = "Unable to profile nuisance coefficients for the empirical-likelihood regression fit."))
  }

  beta <- el_regression_beta(delta, optimizer$par, design$term_index, p)
  moments <- el_regression_moments(design, beta)
  dual <- el_dual_fit(moments, tol = tol)

  if(!isTRUE(dual$success)) {
    return(list(success = FALSE,
                feasible = FALSE,
                statistic = Inf,
                error = "Unable to solve the empirical-likelihood dual problem for the constrained treatment effect."))
  }

  list(
    success = TRUE,
    feasible = TRUE,
    beta = beta,
    nuisance = optimizer$par,
    statistic = dual$statistic,
    dual = dual
  )
}

el_regression_statistic <- function(design, delta, tol = 1e-8) {
  profile <- el_regression_profile_fit(design, delta, tol = tol)

  if(!isTRUE(profile$success)) {
    return(Inf)
  }

  as.numeric(profile$statistic)
}

el_regression_variance <- function(design) {
  n <- nrow(design$X)
  p <- ncol(design$X)
  residuals <- as.numeric(design$y - design$X %*% design$coefficients)
  rss <- sum(residuals ^ 2)

  if(!is.finite(rss) || rss <= 0 || n <= p) {
    return(NA_real_)
  }

  xtx_inverse <- tryCatch(
    solve(crossprod(design$X)),
    error = function(e) NULL
  )

  if(is.null(xtx_inverse)) {
    return(NA_real_)
  }

  sigma2 <- rss / (n - p)
  variance <- sigma2 * xtx_inverse[design$term_index, design$term_index]
  as.numeric(variance)
}

el_regression_ci_bound <- function(design, estimate, crit, direction, se, tol = 1e-8) {
  step <- max(0.5, 2 * se, abs(estimate) * 0.25)
  if(!is.finite(step) || step <= 0) {
    step <- 1
  }

  inner <- estimate

  for(i in seq_len(40)) {
    candidate <- estimate + direction * step
    value <- suppressWarnings(el_regression_statistic(design, candidate, tol = tol) - crit)

    if(!is.finite(value)) {
      outer <- candidate
      for(j in seq_len(80)) {
        midpoint <- (inner + outer) / 2
        midpoint_value <- suppressWarnings(el_regression_statistic(design, midpoint, tol = tol) - crit)

        if(is.finite(midpoint_value) && midpoint_value >= 0) {
          outer <- midpoint
          break
        }

        if(is.finite(midpoint_value)) {
          inner <- midpoint
        } else {
          outer <- midpoint
        }
      }

      left <- min(inner, outer)
      right <- max(inner, outer)
      for(j in seq_len(120)) {
        midpoint <- (left + right) / 2
        midpoint_value <- suppressWarnings(el_regression_statistic(design, midpoint, tol = tol) - crit)

        if(!is.finite(midpoint_value) || midpoint_value >= 0) {
          right <- midpoint
        } else {
          left <- midpoint
        }

        if(abs(right - left) <= tol) {
          break
        }
      }

      return(as.numeric((left + right) / 2))
    }

    if(value >= 0) {
      left <- min(inner, candidate)
      right <- max(inner, candidate)
      for(j in seq_len(120)) {
        midpoint <- (left + right) / 2
        midpoint_value <- suppressWarnings(el_regression_statistic(design, midpoint, tol = tol) - crit)

        if(!is.finite(midpoint_value) || midpoint_value >= 0) {
          right <- midpoint
        } else {
          left <- midpoint
        }

        if(abs(right - left) <= tol) {
          break
        }
      }

      return(as.numeric((left + right) / 2))
    }

    inner <- candidate
    step <- step * 2
  }

  NA_real_
}

el_regression_confint <- function(design, estimate, conf.level, tol = 1e-8) {
  variance <- el_regression_variance(design)
  if(!is.finite(variance) || variance <= 0) {
    return(c(NA_real_, NA_real_))
  }

  se <- sqrt(variance)
  crit <- stats::qchisq(conf.level, 1)
  lower <- el_regression_ci_bound(design, estimate, crit, -1, se, tol = tol)
  upper <- el_regression_ci_bound(design, estimate, crit, 1, se, tol = tol)

  if(!all(is.finite(c(lower, upper)))) {
    return(c(NA_real_, NA_real_))
  }

  conf_int <- c(lower, upper)
  attr(conf_int, "conf.level") <- conf.level
  conf_int
}

el_regression_fit <- function(data, formula, term = "R", mu = 0, conf.level = 0.95, tol = 1e-8) {
  conf.level <- validateConfidenceLevel(conf.level)
  design <- el_regression_design(formula, data, term = term)

  if(!isTRUE(design$success)) {
    return(list(success = FALSE, error = design$error))
  }

  if(identical(design$X_names, c("(Intercept)", "R"))) {
    x <- design$y[design$X[, design$term_index] == 1]
    y <- design$y[design$X[, design$term_index] == 0]
    mean_fit <- el_mean_diff_fit(x, y, mu = mu, conf.level = conf.level)
    profile_at_estimate <- el_regression_profile_fit(design, as.numeric(mean_fit$estimate), tol = tol)

    return(list(
      success = TRUE,
      estimate = as.numeric(mean_fit$estimate),
      conf.int = as.numeric(mean_fit$conf.int),
      statistic = as.numeric(mean_fit$statistic),
      p.value = as.numeric(mean_fit$p.value),
      null.value = mu,
      design = design,
      profile = profile_at_estimate
    ))
  }

  estimate <- as.numeric(design$estimate)
  statistic <- el_regression_statistic(design, mu, tol = tol)

  if(!is.finite(statistic)) {
    return(list(success = FALSE,
                error = "Adjusted semi-parametric EL is not estimable under the supplied covariates."))
  }

  conf_int <- el_regression_confint(design, estimate, conf.level, tol = tol)
  if(!all(is.finite(conf_int))) {
    return(list(success = FALSE,
                error = "Adjusted semi-parametric EL is not estimable under the supplied covariates."))
  }

  p_value <- stats::pchisq(statistic, 1, lower.tail = FALSE)
  profile_at_estimate <- el_regression_profile_fit(design, estimate, tol = tol)

  list(
    success = TRUE,
    estimate = estimate,
    conf.int = conf_int,
    statistic = statistic,
    p.value = p_value,
    null.value = mu,
    design = design,
    profile = profile_at_estimate
  )
}
