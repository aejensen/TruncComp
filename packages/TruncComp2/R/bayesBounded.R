bayes_bounded_supports <- function() {
  c("bounded_continuous", "bounded_score")
}

bayes_is_bounded_support <- function(continuous_support) {
  continuous_support %in% bayes_bounded_supports()
}

bayes_fit_bounded_kernel <- function(object) {
  if(is.null(object$settings$bounded_kernel)) {
    return("beta")
  }

  bayes_bounded_kernel(object$settings$bounded_kernel)
}

bayes_validate_bounded_range <- function(score_min, score_max) {
  if(!(length(score_min) == 1L &&
       is.numeric(score_min) &&
       is.finite(score_min))) {
    stop("score_min must be a single finite numeric value.")
  }

  if(!(length(score_max) == 1L &&
       is.numeric(score_max) &&
       is.finite(score_max))) {
    stop("score_max must be a single finite numeric value.")
  }

  if(!(score_min < score_max)) {
    stop("score_min must be strictly less than score_max.")
  }

  list(
    score_min = as.numeric(score_min),
    score_max = as.numeric(score_max),
    score_range = as.numeric(score_max - score_min)
  )
}

bayes_score_grid_tolerance <- function(score_min, score_max, score_step = 1) {
  max(1e-8, 1e-8 * max(1, abs(score_min), abs(score_max), abs(score_step)))
}

bayes_is_grid_multiple <- function(x, step, tol = bayes_score_grid_tolerance(0, x, step)) {
  position <- x / step
  abs(position - round(position)) <= tol
}

bayes_validate_positive_scalar <- function(x, name) {
  if(!(length(x) == 1L &&
       is.numeric(x) &&
       is.finite(x) &&
       x > 0)) {
    stop(name, " must be a single positive finite numeric value.")
  }

  as.numeric(x)
}

bayes_score_values <- function(score_min, score_max, score_step) {
  score_step <- bayes_validate_positive_scalar(score_step, "score_step")
  tol <- bayes_score_grid_tolerance(score_min, score_max, score_step)
  n_steps <- (score_max - score_min) / score_step

  if(!(is.finite(n_steps) && bayes_is_grid_multiple(score_max - score_min, score_step, tol))) {
    stop("score_max - score_min must be an integer multiple of score_step.")
  }

  n_steps <- as.integer(round(n_steps))
  values <- score_min + score_step * seq.int(0L, n_steps)
  values[length(values)] <- score_max
  as.numeric(values)
}

bayes_score_grid_index <- function(y, score_min, score_max, score_step) {
  score_values <- bayes_score_values(score_min, score_max, score_step)
  tol <- bayes_score_grid_tolerance(score_min, score_max, score_step)
  position <- (y - score_min) / score_step
  index <- as.integer(round(position)) + 1L
  valid <- is.finite(y) &
    index >= 1L &
    index <= length(score_values) &
    abs(y - score_values[pmin(pmax(index, 1L), length(score_values))]) <= tol

  if(any(!valid)) {
    if(isTRUE(all.equal(score_step, 1))) {
      stop(
        "For continuous_support = \"bounded_score\" with score_step = 1, all observed non-atom scores must be integers on the score grid."
      )
    }

    stop("For continuous_support = \"bounded_score\", all observed non-atom scores must lie exactly on the score grid.")
  }

  index
}

bayes_score_cell_counts <- function(arm_obs, y_obs_index, j_count) {
  j_count <- validate_bayes_positive_integer(j_count, "j_count", min_value = 1L)
  arm_obs <- as.integer(arm_obs)
  y_obs_index <- as.integer(y_obs_index)

  if(length(arm_obs) != length(y_obs_index)) {
    stop("arm_obs and y_obs_index must have the same length.")
  }
  if(length(arm_obs) == 0L) {
    stop("At least one observed non-atom score is required.")
  }
  if(any(is.na(arm_obs)) || any(!arm_obs %in% 1:2)) {
    stop("arm_obs must contain arm indices 1 or 2.")
  }
  if(any(is.na(y_obs_index)) || any(y_obs_index < 1L | y_obs_index > j_count)) {
    stop("y_obs_index must contain valid score-grid indices.")
  }

  cell_id <- (arm_obs - 1L) * j_count + y_obs_index
  counts <- tabulate(cell_id, nbins = 2L * j_count)
  nonzero <- which(counts > 0L)

  list(
    N_score_cells = as.integer(length(nonzero)),
    score_cell_arm = as.integer((nonzero - 1L) %/% j_count + 1L),
    score_cell_index = as.integer((nonzero - 1L) %% j_count + 1L),
    score_cell_count = as.integer(counts[nonzero])
  )
}

bayes_normalize_heaping <- function(heaping = "shared") {
  match.arg(heaping, c("shared", "arm_specific"))
}

bayes_validate_heaping_grids <- function(heaping_grids, score_step,
                                         score_min, score_max) {
  if(!(is.numeric(heaping_grids) &&
       length(heaping_grids) >= 1L &&
       all(is.finite(heaping_grids)) &&
       all(heaping_grids > 0))) {
    stop("heaping_grids must contain positive finite numeric values.")
  }

  heaping_grids <- sort(unique(as.numeric(heaping_grids)))
  tol <- bayes_score_grid_tolerance(score_min, score_max, score_step)
  multiples <- vapply(
    heaping_grids,
    function(grid) bayes_is_grid_multiple(grid, score_step, tol),
    logical(1)
  )

  if(!all(multiples)) {
    stop("Every value in heaping_grids must be a positive multiple of score_step.")
  }

  heaping_grids
}

bayes_score_compatibility_matrix <- function(score_values, heaping_grids,
                                             score_min, tol) {
  valid <- outer(
    heaping_grids,
    score_values,
    Vectorize(function(grid, value) {
      bayes_is_grid_multiple(value - score_min, grid, tol)
    })
  )

  matrix(as.integer(valid), nrow = length(heaping_grids), ncol = length(score_values))
}

bayes_score_bin_metadata <- function(score_min, score_max, score_step,
                                     heaping_grids) {
  score_values <- bayes_score_values(score_min, score_max, score_step)
  score_range <- score_max - score_min
  tol <- bayes_score_grid_tolerance(score_min, score_max, score_step)
  bin_valid <- bayes_score_compatibility_matrix(
    score_values = score_values,
    heaping_grids = heaping_grids,
    score_min = score_min,
    tol = tol
  )
  bin_lower <- matrix(0, nrow = length(heaping_grids), ncol = length(score_values))
  bin_upper <- matrix(0, nrow = length(heaping_grids), ncol = length(score_values))

  for(k in seq_along(heaping_grids)) {
    grid <- heaping_grids[[k]]

    for(j in seq_along(score_values)) {
      if(bin_valid[k, j] == 0L) {
        next
      }

      lower <- max(score_min, score_values[[j]] - grid / 2)
      upper <- min(score_max, score_values[[j]] + grid / 2)
      bin_lower[k, j] <- (lower - score_min) / score_range
      bin_upper[k, j] <- (upper - score_min) / score_range
    }
  }

  list(
    score_values = as.numeric(score_values),
    bin_lower = bin_lower,
    bin_upper = bin_upper,
    bin_valid = bin_valid
  )
}

bayes_validate_eta_prior <- function(eta_prior, k) {
  if(!(is.numeric(eta_prior) &&
       length(eta_prior) %in% c(1L, k) &&
       all(is.finite(eta_prior)) &&
       all(eta_prior > 0))) {
    stop("prior$eta_prior must be a positive finite scalar or have one value per heaping grid.")
  }

  if(length(eta_prior) == 1L) {
    eta_prior <- rep(as.numeric(eta_prior), k)
  }

  as.numeric(eta_prior)
}

bayes_support_supplied_defaults <- function(score_min = FALSE,
                                            score_max = FALSE,
                                            score_step = FALSE,
                                            heaping_grids = FALSE,
                                            heaping = FALSE) {
  list(
    score_min = isTRUE(score_min),
    score_max = isTRUE(score_max),
    score_step = isTRUE(score_step),
    heaping_grids = isTRUE(heaping_grids),
    heaping = isTRUE(heaping)
  )
}

bayes_normalize_support_options <- function(continuous_support,
                                            score_min = NULL,
                                            score_max = NULL,
                                            score_step = 1,
                                            heaping_grids = 1,
                                            heaping = "shared",
                                            supplied = bayes_support_supplied_defaults()) {
  continuous_support <- bayes_continuous_support(continuous_support)

  if(!bayes_is_bounded_support(continuous_support)) {
    supplied_any <- any(unlist(supplied, use.names = FALSE))
    if(supplied_any) {
      stop(
        "score_min, score_max, score_step, heaping_grids, and heaping are only supported for bounded Bayesian models."
      )
    }

    return(list(continuous_support = continuous_support))
  }

  range <- bayes_validate_bounded_range(score_min, score_max)

  if(identical(continuous_support, "bounded_continuous")) {
    if(isTRUE(supplied$score_step) ||
       isTRUE(supplied$heaping_grids) ||
       isTRUE(supplied$heaping)) {
      stop(
        "score_step, heaping_grids, and heaping are only used with continuous_support = \"bounded_score\"."
      )
    }

    return(c(
      list(continuous_support = continuous_support),
      range
    ))
  }

  score_step <- bayes_validate_positive_scalar(score_step, "score_step")
  score_values <- bayes_score_values(
    score_min = range$score_min,
    score_max = range$score_max,
    score_step = score_step
  )
  heaping_grids <- bayes_validate_heaping_grids(
    heaping_grids = heaping_grids,
    score_step = score_step,
    score_min = range$score_min,
    score_max = range$score_max
  )
  heaping <- bayes_normalize_heaping(heaping)
  bin_metadata <- bayes_score_bin_metadata(
    score_min = range$score_min,
    score_max = range$score_max,
    score_step = score_step,
    heaping_grids = heaping_grids
  )

  eta_groups <- if(identical(heaping, "shared")) 1L else 2L
  eta_group_by_arm <- if(identical(heaping, "shared")) c(1L, 1L) else c(1L, 2L)

  c(
    list(
      continuous_support = continuous_support,
      score_step = score_step,
      score_values = score_values,
      heaping_grids = heaping_grids,
      heaping = heaping,
      eta_groups = eta_groups,
      eta_group_by_arm = eta_group_by_arm
    ),
    range,
    bin_metadata[c("bin_lower", "bin_upper", "bin_valid")]
  )
}

bayes_validate_bounded_data <- function(data, atom, support_options) {
  observed <- data$Y[data$A == 1L]
  score_min <- support_options$score_min
  score_max <- support_options$score_max

  if(score_min <= atom && atom <= score_max) {
    stop(
      "For bounded Bayesian models, atom must not lie on the non-atom outcome support."
    )
  }

  if(any(observed < score_min | observed > score_max)) {
    stop(
      "For bounded Bayesian models, all observed non-atom outcomes must lie in [score_min, score_max]."
    )
  }

  if(identical(support_options$continuous_support, "bounded_continuous")) {
    if(any(observed <= score_min | observed >= score_max)) {
      stop(
        "For continuous_support = \"bounded_continuous\", observed non-atom outcomes must be strictly inside (score_min, score_max). Boundary values require the bounded_score model in this MVP."
      )
    }

    return(invisible(TRUE))
  }

  bayes_score_grid_index(
    observed,
    score_min = score_min,
    score_max = score_max,
    score_step = support_options$score_step
  )

  compatible <- bayes_score_compatibility_matrix(
    score_values = observed,
    heaping_grids = support_options$heaping_grids,
    score_min = score_min,
    tol = bayes_score_grid_tolerance(score_min, score_max, support_options$score_step)
  )

  if(any(colSums(compatible) == 0L)) {
    stop("Every observed non-atom score must be compatible with at least one heaping grid.")
  }

  invisible(TRUE)
}

bayes_log_diff_exp <- function(log_x, log_y) {
  if(log_y > log_x) {
    return(NaN)
  }

  log_x + log1p(-exp(log_y - log_x))
}

bayes_beta_interval_prob <- function(lower, upper, shape1, shape2) {
  if(lower <= 0 && upper >= 1) {
    return(1)
  }

  if(lower <= 0) {
    return(stats::pbeta(upper, shape1, shape2))
  }

  if(upper >= 1) {
    return(stats::pbeta(lower, shape1, shape2, lower.tail = FALSE))
  }

  out <- exp(bayes_log_diff_exp(
    stats::pbeta(upper, shape1, shape2, log.p = TRUE),
    stats::pbeta(lower, shape1, shape2, log.p = TRUE)
  ))

  max(0, min(1, out))
}

bayes_logitnormal_mean_quad <- local({
  cache <- new.env(parent = emptyenv())

  function(n = 20L) {
    n <- validate_bayes_positive_integer(n, "n", min_value = 2L)
    key <- as.character(n)

    if(exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }

    index <- seq_len(n - 1L)
    jacobi <- matrix(0, nrow = n, ncol = n)
    off_diagonal <- sqrt(index / 2)
    jacobi[cbind(index, index + 1L)] <- off_diagonal
    jacobi[cbind(index + 1L, index)] <- off_diagonal

    eig <- eigen(jacobi, symmetric = TRUE)
    order <- order(eig$values)
    nodes <- eig$values[order]
    weights <- eig$vectors[1L, order]^2
    weights <- weights / sum(weights)

    out <- list(
      n = as.integer(n),
      nodes = as.numeric(nodes),
      weights = as.numeric(weights)
    )
    assign(key, out, envir = cache)
    out
  }
})

bayes_logitnormal_unit_mean <- function(mu, sigma, n = 20L) {
  sigma <- bayes_validate_positive_scalar(sigma, "sigma")
  quad <- bayes_logitnormal_mean_quad(n)
  sum(quad$weights * stats::plogis(mu + sqrt(2) * sigma * quad$nodes))
}

bayes_logitnormal_density <- function(x, mu, sigma, score_min, score_max) {
  sigma <- bayes_validate_positive_scalar(sigma, "sigma")
  score_range <- score_max - score_min
  x_unit <- (x - score_min) / score_range
  in_bounds <- x_unit > 0 & x_unit < 1
  x_eval <- pmin(pmax(x_unit, .Machine$double.eps^0.25), 1 - .Machine$double.eps^0.25)
  out <- stats::dnorm(stats::qlogis(x_eval), mean = mu, sd = sigma) /
    (x_eval * (1 - x_eval) * score_range)
  out[!in_bounds] <- 0
  out
}

bayes_logitnormal_cdf <- function(x, mu, sigma, score_min, score_max) {
  sigma <- bayes_validate_positive_scalar(sigma, "sigma")
  x_unit <- (x - score_min) / (score_max - score_min)
  out <- numeric(length(x_unit))
  out[x_unit >= 1] <- 1
  inside <- x_unit > 0 & x_unit < 1
  out[inside] <- stats::pnorm(stats::qlogis(x_unit[inside]), mean = mu, sd = sigma)
  out
}

bayes_logitnormal_random <- function(mu, sigma, score_min, score_max) {
  if(!(is.numeric(mu) && is.numeric(sigma) && length(mu) == length(sigma))) {
    stop("mu and sigma must be numeric vectors with the same length.")
  }
  if(!(all(is.finite(mu)) && all(is.finite(sigma)) && all(sigma > 0))) {
    stop("mu must be finite and sigma must be positive and finite.")
  }

  score_min + (score_max - score_min) *
    stats::plogis(stats::rnorm(length(mu), mean = mu, sd = sigma))
}

bayes_logitnormal_interval_prob <- function(lower, upper, mu, sigma) {
  sigma <- bayes_validate_positive_scalar(sigma, "sigma")

  if(upper <= lower) {
    return(0)
  }

  if(lower <= 0 && upper >= 1) {
    return(1)
  }

  if(lower <= 0) {
    return(stats::pnorm(stats::qlogis(upper), mean = mu, sd = sigma))
  }

  if(upper >= 1) {
    return(stats::pnorm(stats::qlogis(lower), mean = mu, sd = sigma, lower.tail = FALSE))
  }

  out <- exp(bayes_log_diff_exp(
    stats::pnorm(stats::qlogis(upper), mean = mu, sd = sigma, log.p = TRUE),
    stats::pnorm(stats::qlogis(lower), mean = mu, sd = sigma, log.p = TRUE)
  ))

  max(0, min(1, out))
}

bayes_score_pmf <- function(weights, m_comp = NULL, phi_comp = NULL, eta,
                            bin_lower, bin_upper, bin_valid,
                            bounded_kernel = c("beta", "logit_normal"),
                            mu_logit_comp = NULL, sigma_logit_comp = NULL) {
  bounded_kernel <- bayes_bounded_kernel(bounded_kernel)
  h_count <- length(weights)
  k_count <- length(eta)
  j_count <- ncol(bin_valid)
  out <- numeric(j_count)

  for(h in seq_len(h_count)) {
    if(identical(bounded_kernel, "beta")) {
      shape1 <- m_comp[[h]] * phi_comp[[h]]
      shape2 <- (1 - m_comp[[h]]) * phi_comp[[h]]
    } else {
      mu_logit <- mu_logit_comp[[h]]
      sigma_logit <- sigma_logit_comp[[h]]
    }
    grid_pmf <- numeric(j_count)

    for(k in seq_len(k_count)) {
      valid_j <- which(bin_valid[k, ] == 1L)
      if(length(valid_j) == 0L) {
        next
      }

      interval_prob <- vapply(
        valid_j,
        function(j) {
          if(identical(bounded_kernel, "beta")) {
            return(bayes_beta_interval_prob(
              lower = bin_lower[k, j],
              upper = bin_upper[k, j],
              shape1 = shape1,
              shape2 = shape2
            ))
          }

          bayes_logitnormal_interval_prob(
            lower = bin_lower[k, j],
            upper = bin_upper[k, j],
            mu = mu_logit,
            sigma = sigma_logit
          )
        },
        numeric(1)
      )
      grid_pmf[valid_j] <- grid_pmf[valid_j] + eta[[k]] * interval_prob
    }

    out <- out + weights[[h]] * grid_pmf
  }

  out
}
