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
