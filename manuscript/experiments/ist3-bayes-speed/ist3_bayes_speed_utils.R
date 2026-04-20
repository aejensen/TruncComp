ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

prior_value <- function(prior, name, default) {
  if (is.null(prior[[name]])) default else prior[[name]]
}

score_grid_tolerance <- function(score_min, score_max, score_step) {
  max(1e-8, 1e-8 * max(1, abs(score_min), abs(score_max), abs(score_step)))
}

is_grid_multiple <- function(x, step, tol) {
  position <- x / step
  abs(position - round(position)) <= tol
}

score_values <- function(score_min, score_max, score_step) {
  tol <- score_grid_tolerance(score_min, score_max, score_step)
  n_steps <- (score_max - score_min) / score_step
  if (!is.finite(n_steps) || !is_grid_multiple(score_max - score_min, score_step, tol)) {
    stop("score_max - score_min must be an integer multiple of score_step.", call. = FALSE)
  }
  values <- score_min + score_step * seq.int(0L, as.integer(round(n_steps)))
  values[length(values)] <- score_max
  as.numeric(values)
}

score_grid_index <- function(y, score_min, score_max, score_step) {
  values <- score_values(score_min, score_max, score_step)
  tol <- score_grid_tolerance(score_min, score_max, score_step)
  index <- as.integer(round((y - score_min) / score_step)) + 1L
  valid <- is.finite(y) &
    index >= 1L &
    index <= length(values) &
    abs(y - values[pmin(pmax(index, 1L), length(values))]) <= tol
  if (any(!valid)) {
    stop("Observed survivor scores must lie exactly on the bounded-score grid.", call. = FALSE)
  }
  index
}

score_compatibility_matrix <- function(values, heaping_grids, score_min, tol) {
  valid <- outer(
    heaping_grids,
    values,
    Vectorize(function(grid, value) is_grid_multiple(value - score_min, grid, tol))
  )
  matrix(as.integer(valid), nrow = length(heaping_grids), ncol = length(values))
}

score_bin_metadata <- function(score_min, score_max, score_step, heaping_grids) {
  values <- score_values(score_min, score_max, score_step)
  score_range <- score_max - score_min
  tol <- score_grid_tolerance(score_min, score_max, score_step)
  bin_valid <- score_compatibility_matrix(values, heaping_grids, score_min, tol)
  bin_lower <- matrix(0, nrow = length(heaping_grids), ncol = length(values))
  bin_upper <- matrix(0, nrow = length(heaping_grids), ncol = length(values))

  for (k in seq_along(heaping_grids)) {
    grid <- heaping_grids[[k]]
    for (j in seq_along(values)) {
      if (bin_valid[k, j] == 0L) next
      lower <- max(score_min, values[[j]] - grid / 2)
      upper <- min(score_max, values[[j]] + grid / 2)
      bin_lower[k, j] <- (lower - score_min) / score_range
      bin_upper[k, j] <- (upper - score_min) / score_range
    }
  }

  list(
    score_value = as.numeric(values),
    bin_lower = bin_lower,
    bin_upper = bin_upper,
    bin_valid = bin_valid
  )
}

beta_interval_prob <- function(lower, upper, shape1, shape2) {
  if (lower <= 0 && upper >= 1) return(1)
  if (lower <= 0) return(stats::pbeta(upper, shape1, shape2))
  if (upper >= 1) return(stats::pbeta(lower, shape1, shape2, lower.tail = FALSE))
  log_upper <- stats::pbeta(upper, shape1, shape2, log.p = TRUE)
  log_lower <- stats::pbeta(lower, shape1, shape2, log.p = TRUE)
  out <- exp(log_upper + log1p(-exp(log_lower - log_upper)))
  max(0, min(1, out))
}

score_pmf <- function(weights, m_comp, phi_comp, eta, bin_lower, bin_upper, bin_valid) {
  out <- numeric(ncol(bin_valid))
  for (h in seq_along(weights)) {
    shape1 <- m_comp[[h]] * phi_comp[[h]]
    shape2 <- (1 - m_comp[[h]]) * phi_comp[[h]]
    grid_pmf <- numeric(ncol(bin_valid))
    for (k in seq_along(eta)) {
      valid_j <- which(bin_valid[k, ] == 1L)
      if (!length(valid_j)) next
      interval_prob <- vapply(
        valid_j,
        function(j) beta_interval_prob(bin_lower[k, j], bin_upper[k, j], shape1, shape2),
        numeric(1)
      )
      grid_pmf[valid_j] <- grid_pmf[valid_j] + eta[[k]] * interval_prob
    }
    out <- out + weights[[h]] * grid_pmf
  }
  pmax(out, .Machine$double.xmin)
}

logitnormal_interval_prob <- function(lower, upper, mu, sigma) {
  if (lower <= 0 && upper >= 1) return(1)
  if (lower <= 0) return(stats::pnorm(stats::qlogis(upper), mean = mu, sd = sigma))
  if (upper >= 1) return(stats::pnorm(stats::qlogis(lower), mean = mu, sd = sigma, lower.tail = FALSE))
  log_upper <- stats::pnorm(stats::qlogis(upper), mean = mu, sd = sigma, log.p = TRUE)
  log_lower <- stats::pnorm(stats::qlogis(lower), mean = mu, sd = sigma, log.p = TRUE)
  out <- exp(log_upper + log1p(-exp(log_lower - log_upper)))
  max(0, min(1, out))
}

score_pmf_logitnormal <- function(weights, mu_logit_comp, sigma_comp, eta,
                                  bin_lower, bin_upper, bin_valid) {
  out <- numeric(ncol(bin_valid))
  for (h in seq_along(weights)) {
    grid_pmf <- numeric(ncol(bin_valid))
    for (k in seq_along(eta)) {
      valid_j <- which(bin_valid[k, ] == 1L)
      if (!length(valid_j)) next
      interval_prob <- vapply(
        valid_j,
        function(j) logitnormal_interval_prob(
          bin_lower[k, j],
          bin_upper[k, j],
          mu_logit_comp[[h]],
          sigma_comp[[h]]
        ),
        numeric(1)
      )
      grid_pmf[valid_j] <- grid_pmf[valid_j] + eta[[k]] * interval_prob
    }
    out <- out + weights[[h]] * grid_pmf
  }
  pmax(out, .Machine$double.xmin)
}

score_pmf_kernel <- function(weights, eta, bin_lower, bin_upper, bin_valid,
                             kernel = "beta",
                             m_comp = NULL, phi_comp = NULL,
                             mu_logit_comp = NULL, sigma_comp = NULL) {
  if (identical(kernel, "logitnormal")) {
    return(score_pmf_logitnormal(
      weights = weights,
      mu_logit_comp = mu_logit_comp,
      sigma_comp = sigma_comp,
      eta = eta,
      bin_lower = bin_lower,
      bin_upper = bin_upper,
      bin_valid = bin_valid
    ))
  }

  score_pmf(
    weights = weights,
    m_comp = m_comp,
    phi_comp = phi_comp,
    eta = eta,
    bin_lower = bin_lower,
    bin_upper = bin_upper,
    bin_valid = bin_valid
  )
}

score_cell_counts <- function(arm_obs, y_obs_index, j_count) {
  counts <- as.data.frame(table(
    arm = factor(arm_obs, levels = 1:2),
    score = factor(y_obs_index, levels = seq_len(j_count))
  ))
  counts <- counts[counts$Freq > 0L, , drop = FALSE]
  data.frame(
    arm = as.integer(as.character(counts$arm)),
    score = as.integer(as.character(counts$score)),
    count = as.integer(counts$Freq)
  )
}

load_ist3_experiment_data <- function(repo_root) {
  path <- file.path(repo_root, "manuscript", "application-data", "ist3-standardized.csv")
  d <- utils::read.csv(path)
  required <- c("Y", "R", "A", "atom")
  missing <- setdiff(required, names(d))
  if (length(missing)) {
    stop("Missing columns in IST-3 standardized CSV: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  d[, required, drop = FALSE]
}

build_ist3_standata <- function(d, prior, H = 2L,
                                score_min = 0, score_max = 100, score_step = 1,
                                heaping_grids = c(1, 5, 10),
                                heaping = "shared",
                                grainsize = 128L) {
  metadata <- score_bin_metadata(score_min, score_max, score_step, heaping_grids)
  observed <- d[d$A == 1L, , drop = FALSE]
  arm_obs <- as.integer(observed$R) + 1L
  y_obs_index <- score_grid_index(observed$Y, score_min, score_max, score_step)
  cells <- score_cell_counts(arm_obs, y_obs_index, length(metadata$score_value))
  eta_groups <- if (identical(heaping, "shared")) 1L else 2L
  eta_group_by_arm <- if (identical(heaping, "shared")) c(1L, 1L) else c(1L, 2L)

  list(
    H = as.integer(H),
    J = as.integer(length(metadata$score_value)),
    score_value = as.numeric(metadata$score_value),
    K = as.integer(length(heaping_grids)),
    bin_lower = metadata$bin_lower,
    bin_upper = metadata$bin_upper,
    bin_valid = metadata$bin_valid,
    eta_groups = as.integer(eta_groups),
    eta_group_by_arm = as.integer(eta_group_by_arm),
    N_score_cells = as.integer(nrow(cells)),
    score_cell_arm = as.integer(cells$arm),
    score_cell_index = as.integer(cells$score),
    score_cell_count = as.integer(cells$count),
    n_arm = as.integer(vapply(0:1, function(r) sum(d$R == r), integer(1))),
    n_obs_arm = as.integer(vapply(0:1, function(r) sum(d$R == r & d$A == 1L), integer(1))),
    grainsize = as.integer(grainsize),
    atom = as.numeric(unique(d$atom)),
    rho_prior_alpha = prior$rho_alpha,
    rho_prior_beta = prior$rho_beta,
    alpha_prior_shape = prior$alpha_shape,
    alpha_prior_rate = prior$alpha_rate,
    m_prior_alpha = prior$m_alpha,
    m_prior_beta = prior$m_beta,
    phi_prior_meanlog = prior$phi_meanlog,
    phi_prior_sdlog = prior$phi_sdlog,
    gap_prior_alpha = prior_value(prior, "gap_prior_alpha", 1),
    gap_prior_beta = prior_value(prior, "gap_prior_beta", 1),
    use_fixed_v_prior = as.integer(isTRUE(prior_value(prior, "use_fixed_v_prior", FALSE))),
    v_prior_alpha = prior_value(prior, "v_prior_alpha", 2),
    v_prior_beta = prior_value(prior, "v_prior_beta", 2),
    eta_prior = eta_prior_for_k(prior$eta_prior, length(heaping_grids)),
    .kernel = prior_value(prior, "kernel", "beta"),
    .observed = observed,
    .arm_obs = arm_obs,
    .y_obs_index = y_obs_index,
    .cells = cells
  )
}

stan_data_only <- function(standata) {
  standata[!startsWith(names(standata), ".")]
}

gauss_legendre_01 <- function(n) {
  if (n < 2L) {
    stop("Gauss-Legendre quadrature needs at least two nodes.", call. = FALSE)
  }
  i <- seq_len(n - 1L)
  off_diag <- i / sqrt(4 * i^2 - 1)
  jacobi <- matrix(0, nrow = n, ncol = n)
  jacobi[cbind(i, i + 1L)] <- off_diag
  jacobi[cbind(i + 1L, i)] <- off_diag
  eig <- eigen(jacobi, symmetric = TRUE)
  ord <- order(eig$values)
  x <- eig$values[ord]
  w <- 2 * eig$vectors[1L, ord]^2
  list(
    node = as.numeric((x + 1) / 2),
    weight = as.numeric(w / 2)
  )
}

add_stick_quadrature <- function(standata, n = 64L) {
  quad <- gauss_legendre_01(as.integer(n))
  standata$v_quad_n <- as.integer(n)
  standata$v_quad_node <- quad$node
  standata$v_quad_log_weight <- log(quad$weight)
  standata
}

eta_prior_for_k <- function(eta_prior, k) {
  eta_prior <- as.numeric(eta_prior)
  if (length(eta_prior) == 1L) {
    return(rep(eta_prior, k))
  }
  if (length(eta_prior) == k) {
    return(eta_prior)
  }
  if (k == 1L) {
    return(eta_prior[[1]])
  }
  stop("eta_prior must be length 1 or match the number of heaping grids.", call. = FALSE)
}

ist3_prior_candidates <- function() {
  list(
    regularized_balanced = list(
      label = "regularized_balanced",
      rho_alpha = 1, rho_beta = 1,
      alpha_shape = 2, alpha_rate = 1,
      m_alpha = 2, m_beta = 2,
      phi_meanlog = log(8), phi_sdlog = 0.5,
      eta_prior = c(2, 2, 2),
      gap_prior_alpha = 1, gap_prior_beta = 1,
      use_fixed_v_prior = FALSE,
      v_prior_alpha = 2, v_prior_beta = 2
    ),
    regularized_shifted = list(
      label = "regularized_shifted",
      rho_alpha = 1, rho_beta = 1,
      alpha_shape = 2, alpha_rate = 1,
      m_alpha = 2.5, m_beta = 2,
      phi_meanlog = log(6), phi_sdlog = 0.5,
      eta_prior = c(2, 2, 2),
      gap_prior_alpha = 1, gap_prior_beta = 1,
      use_fixed_v_prior = FALSE,
      v_prior_alpha = 2, v_prior_beta = 2
    ),
    balanced_alpha_stable = list(
      label = "balanced_alpha_stable",
      rho_alpha = 1, rho_beta = 1,
      alpha_shape = 4, alpha_rate = 4,
      m_alpha = 2, m_beta = 2,
      phi_meanlog = log(8), phi_sdlog = 0.5,
      eta_prior = c(2, 2, 2),
      gap_prior_alpha = 1, gap_prior_beta = 1,
      use_fixed_v_prior = FALSE,
      v_prior_alpha = 2, v_prior_beta = 2
    ),
    strong_regularized = list(
      label = "strong_regularized",
      rho_alpha = 1, rho_beta = 1,
      alpha_shape = 4, alpha_rate = 4,
      m_alpha = 2.5, m_beta = 2,
      phi_meanlog = log(5), phi_sdlog = 0.4,
      eta_prior = c(3, 3, 3),
      gap_prior_alpha = 1, gap_prior_beta = 1,
      use_fixed_v_prior = FALSE,
      v_prior_alpha = 2, v_prior_beta = 2
    )
  )
}

validate_aggregation <- function(standata) {
  stopifnot(sum(standata$score_cell_count) == nrow(standata$.observed))
  set.seed(9001)
  H <- standata$H
  K <- standata$K
  kernel <- prior_value(standata, ".kernel", "beta")
  v <- matrix(stats::rbeta(2 * (H - 1), 2, 2), nrow = 2)
  w <- matrix(NA_real_, nrow = 2, ncol = H)
  for (r in 1:2) {
    remaining <- 1
    for (h in seq_len(H - 1)) {
      w[r, h] <- v[r, h] * remaining
      remaining <- remaining * (1 - v[r, h])
    }
    w[r, H] <- remaining
  }
  m_comp_raw <- matrix(stats::rbeta(2 * H, 2, 2), nrow = 2)
  m_comp <- matrix(NA_real_, nrow = 2, ncol = H)
  for (r in 1:2) {
    m_comp[r, ] <- sort(m_comp_raw[r, ])
  }
  phi_comp <- matrix(stats::rlnorm(2 * H, log(8), 0.5), nrow = 2)
  mu_logit_comp <- matrix(stats::rnorm(2 * H, 0, 1.2), nrow = 2)
  sigma_comp <- matrix(stats::rlnorm(2 * H, log(1.2), 0.25), nrow = 2)
  eta <- matrix(stats::rgamma(standata$eta_groups * K, 2, 1), nrow = standata$eta_groups)
  eta <- eta / rowSums(eta)
  pi <- c(0.7, 0.7)

  pmf <- lapply(1:2, function(r) {
    score_pmf_kernel(
      weights = w[r, ],
      kernel = kernel,
      m_comp = m_comp[r, ],
      phi_comp = phi_comp[r, ],
      mu_logit_comp = mu_logit_comp[r, ],
      sigma_comp = sigma_comp[r, ],
      eta = eta[standata$eta_group_by_arm[[r]], ],
      bin_lower = standata$bin_lower,
      bin_upper = standata$bin_upper,
      bin_valid = standata$bin_valid
    )
  })

  obs_ll <- sum(log(pmf[[1]][standata$.y_obs_index[standata$.arm_obs == 1L]])) +
    sum(log(pmf[[2]][standata$.y_obs_index[standata$.arm_obs == 2L]]))
  agg_ll <- sum(vapply(seq_len(standata$N_score_cells), function(i) {
    standata$score_cell_count[[i]] *
      log(pmf[[standata$score_cell_arm[[i]]]][[standata$score_cell_index[[i]]]])
  }, numeric(1)))

  atom_obs_kernel <- sum(vapply(1:2, function(r) {
    standata$n_obs_arm[[r]] * log(pi[[r]]) +
      (standata$n_arm[[r]] - standata$n_obs_arm[[r]]) * log1p(-pi[[r]])
  }, numeric(1)))

  list(
    n_observed = nrow(standata$.observed),
    n_score_cells = standata$N_score_cells,
    continuous_loglik_diff = unname(obs_ll - agg_ll),
    atom_kernel_loglik = unname(atom_obs_kernel),
    old_beta_interval_work = nrow(standata$.observed) * standata$H * standata$K,
    new_beta_interval_work = standata$N_score_cells * standata$H * standata$K
  )
}

prior_predictive_summary <- function(standata, prior, ndraws = 4000L) {
  set.seed(20260419)
  H <- standata$H
  K <- standata$K
  scores <- standata$score_value
  out <- matrix(NA_real_, nrow = ndraws, ncol = 5)
  colnames(out) <- c("mean", "sd", "p0", "p100", "max_p")
  for (i in seq_len(ndraws)) {
    alpha <- stats::rgamma(1, prior$alpha_shape, prior$alpha_rate)
    if (isTRUE(prior_value(prior, "use_fixed_v_prior", FALSE))) {
      v <- stats::rbeta(
        H - 1,
        prior_value(prior, "v_prior_alpha", 2),
        prior_value(prior, "v_prior_beta", 2)
      )
    } else {
      v <- stats::rbeta(H - 1, 1, alpha)
    }
    w <- numeric(H)
    remaining <- 1
    for (h in seq_len(H - 1)) {
      w[h] <- v[h] * remaining
      remaining <- remaining * (1 - v[h])
    }
    w[H] <- remaining
    if (!identical(prior_value(prior, "m_parameterization", "ordered_gap"), "ordered_logit") &&
        H == 2L &&
        (!identical(prior_value(prior, "gap_prior_alpha", 1), 1) ||
         !identical(prior_value(prior, "gap_prior_beta", 1), 1))) {
      m_low <- stats::rbeta(1, prior$m_alpha, prior$m_beta)
      gap <- stats::rbeta(
        1,
        prior_value(prior, "gap_prior_alpha", 1),
        prior_value(prior, "gap_prior_beta", 1)
      )
      m <- c(m_low, m_low + (1 - 1e-5 - m_low) * gap)
    } else {
      m <- sort(stats::rbeta(H, prior$m_alpha, prior$m_beta))
    }
    phi <- stats::rlnorm(H, prior$phi_meanlog, prior$phi_sdlog)
    mu_logit <- stats::rnorm(H, 0, 2)
    sigma_logit <- stats::rlnorm(H, log(1.2), 0.5)
    eta <- stats::rgamma(K, eta_prior_for_k(prior$eta_prior, K), 1)
    eta <- eta / sum(eta)
    pmf <- score_pmf_kernel(
      weights = w,
      kernel = prior_value(prior, "kernel", "beta"),
      m_comp = m,
      phi_comp = phi,
      mu_logit_comp = mu_logit,
      sigma_comp = sigma_logit,
      eta = eta,
      bin_lower = standata$bin_lower,
      bin_upper = standata$bin_upper,
      bin_valid = standata$bin_valid
    )
    mu <- sum(scores * pmf)
    out[i, ] <- c(mu, sqrt(sum((scores - mu)^2 * pmf)), pmf[[1]], pmf[[length(pmf)]], max(pmf))
  }
  data.frame(
    label = prior$label,
    mean_q05 = stats::quantile(out[, "mean"], 0.05),
    mean_q50 = stats::quantile(out[, "mean"], 0.50),
    mean_q95 = stats::quantile(out[, "mean"], 0.95),
    sd_q05 = stats::quantile(out[, "sd"], 0.05),
    sd_q50 = stats::quantile(out[, "sd"], 0.50),
    sd_q95 = stats::quantile(out[, "sd"], 0.95),
    p0_gt_10 = mean(out[, "p0"] > 0.1),
    p100_gt_10 = mean(out[, "p100"] > 0.1),
    max_p_gt_25 = mean(out[, "max_p"] > 0.25),
    stringsAsFactors = FALSE
  )
}

draw_col <- function(draws, name) {
  if (!name %in% names(draws)) {
    stop("Missing draw column: ", name, call. = FALSE)
  }
  draws[[name]]
}

compute_estimand_draws <- function(fit, standata) {
  draws <- posterior::as_draws_df(fit$draws(format = "draws_df"))
  H <- standata$H
  K <- standata$K
  scores <- standata$score_value
  n <- nrow(draws)
  mu <- matrix(NA_real_, nrow = n, ncol = 2)
  kernel <- prior_value(standata, ".kernel", "beta")
  has_sampled_sticks <- H <= 1L || sprintf("v[%d,%d]", 1L, 1L) %in% names(draws)
  has_marginal_sticks <- !has_sampled_sticks &&
    H == 2L &&
    all(c("v_quad_node", "v_quad_log_weight") %in% names(standata))
  if (!has_sampled_sticks && !has_marginal_sticks) {
    stop("Cannot compute estimands without sampled sticks or stick quadrature metadata.", call. = FALSE)
  }

  for (i in seq_len(n)) {
    for (r in 1:2) {
      if (identical(kernel, "logitnormal")) {
        m <- NULL
        phi <- NULL
        mu_logit <- vapply(seq_len(H), function(h) draws[[sprintf("mu_logit_comp[%d,%d]", r, h)]][[i]], numeric(1))
        sigma_logit <- exp(vapply(seq_len(H), function(h) draws[[sprintf("log_sigma_comp[%d,%d]", r, h)]][[i]], numeric(1)))
      } else {
        m <- vapply(seq_len(H), function(h) draws[[sprintf("m_comp[%d,%d]", r, h)]][[i]], numeric(1))
        phi <- exp(vapply(seq_len(H), function(h) draws[[sprintf("log_phi_comp[%d,%d]", r, h)]][[i]], numeric(1)))
        mu_logit <- NULL
        sigma_logit <- NULL
      }
      eta <- vapply(seq_len(K), function(k) draws[[sprintf("eta[%d,%d]", standata$eta_group_by_arm[[r]], k)]][[i]], numeric(1))

      if (has_sampled_sticks) {
        v <- if (H > 1) vapply(seq_len(H - 1), function(h) draws[[sprintf("v[%d,%d]", r, h)]][[i]], numeric(1)) else numeric(0)
        w <- numeric(H)
        remaining <- 1
        for (h in seq_len(H - 1)) {
          w[h] <- v[[h]] * remaining
          remaining <- remaining * (1 - v[[h]])
        }
        w[H] <- remaining
        pmf <- score_pmf_kernel(
          weights = w,
          kernel = kernel,
          m_comp = m,
          phi_comp = phi,
          mu_logit_comp = mu_logit,
          sigma_comp = sigma_logit,
          eta = eta,
          bin_lower = standata$bin_lower,
          bin_upper = standata$bin_upper,
          bin_valid = standata$bin_valid
        )
      } else {
        component_pmf <- lapply(seq_len(H), function(h) {
          weight <- rep(0, H)
          weight[[h]] <- 1
          score_pmf_kernel(
            weights = weight,
            kernel = kernel,
            m_comp = m,
            phi_comp = phi,
            mu_logit_comp = mu_logit,
            sigma_comp = sigma_logit,
            eta = eta,
            bin_lower = standata$bin_lower,
            bin_upper = standata$bin_upper,
            bin_valid = standata$bin_valid
          )
        })
        arm_cells <- which(standata$score_cell_arm == r)
        alpha_r <- draws[[sprintf("alpha[%d]", r)]][[i]]
        log_q <- standata$v_quad_log_weight +
          stats::dbeta(standata$v_quad_node, 1, alpha_r, log = TRUE)
        for (cell in arm_cells) {
          j <- standata$score_cell_index[[cell]]
          p1 <- component_pmf[[1L]][[j]]
          p2 <- component_pmf[[2L]][[j]]
          log_q <- log_q + standata$score_cell_count[[cell]] *
            log(standata$v_quad_node * p1 + (1 - standata$v_quad_node) * p2)
        }
        q_weight <- exp(log_q - max(log_q))
        q_weight <- q_weight / sum(q_weight)
        w <- c(sum(standata$v_quad_node * q_weight), 1 - sum(standata$v_quad_node * q_weight))
        pmf <- w[[1L]] * component_pmf[[1L]] + w[[2L]] * component_pmf[[2L]]
      }
      mu[i, r] <- sum(scores * pmf)
    }
  }

  rho_0 <- draw_col(draws, "rho[1]")
  rho_1 <- draw_col(draws, "rho[2]")
  pi_0 <- 1 - rho_0
  pi_1 <- 1 - rho_1
  out <- posterior::as_draws_df(data.frame(
    rho_0 = rho_0,
    rho_1 = rho_1,
    pi_0 = pi_0,
    pi_1 = pi_1,
    mu_0_c = mu[, 1],
    mu_1_c = mu[, 2],
    delta_atom = rho_1 - rho_0,
    mu_delta = mu[, 2] - mu[, 1],
    alpha_delta = (pi_1 / rho_1) / (pi_0 / rho_0),
    delta = (standata$atom * rho_1 + pi_1 * mu[, 2]) -
      (standata$atom * rho_0 + pi_0 * mu[, 1]),
    .chain = draws$.chain,
    .iteration = draws$.iteration,
    .draw = draws$.draw
  ))
  out
}

bfmi_by_chain <- function(sampler_df) {
  chains <- split(sampler_df$energy__, sampler_df$.chain)
  vapply(chains, function(e) {
    if (length(e) < 3 || stats::var(e) <= 0) return(NA_real_)
    mean(diff(e)^2) / stats::var(e)
  }, numeric(1))
}

summarize_fit <- function(fit, standata, label, stage, adapt_delta, elapsed, iter_warmup, iter_sampling) {
  sampler <- fit$sampler_diagnostics(format = "df")
  estimands <- compute_estimand_draws(fit, standata)
  key_summary <- posterior::summarise_draws(
    estimands,
    posterior::default_summary_measures(),
    posterior::default_convergence_measures()
  )
  key_summary <- as.data.frame(key_summary)
  bfmi <- bfmi_by_chain(sampler)
  divergences <- sum(sampler$divergent__)
  max_treedepth <- max(sampler$treedepth__)
  min_bfmi <- suppressWarnings(min(bfmi, na.rm = TRUE))
  if (!is.finite(min_bfmi)) min_bfmi <- NA_real_
  max_rhat <- suppressWarnings(max(key_summary$rhat, na.rm = TRUE))
  min_bulk <- suppressWarnings(min(key_summary$ess_bulk, na.rm = TRUE))
  min_tail <- suppressWarnings(min(key_summary$ess_tail, na.rm = TRUE))
  if (!is.finite(max_rhat)) max_rhat <- NA_real_
  if (!is.finite(min_bulk)) min_bulk <- NA_real_
  if (!is.finite(min_tail)) min_tail <- NA_real_
  step_sizes <- unique(sampler$stepsize__)

  metrics <- data.frame(
    stage = stage,
    label = label,
    adapt_delta = adapt_delta,
    chains = length(unique(sampler$.chain)),
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    elapsed_sec = unname(elapsed[["elapsed"]]),
    seconds_per_draw = unname(elapsed[["elapsed"]]) / (length(unique(sampler$.chain)) * iter_sampling),
    divergences = divergences,
    max_treedepth = max_treedepth,
    mean_treedepth = mean(sampler$treedepth__),
    median_leapfrog = stats::median(sampler$n_leapfrog__),
    mean_leapfrog = mean(sampler$n_leapfrog__),
    max_leapfrog = max(sampler$n_leapfrog__),
    min_stepsize = min(step_sizes),
    median_stepsize = stats::median(step_sizes),
    max_stepsize = max(step_sizes),
    min_bfmi = min_bfmi,
    max_rhat = max_rhat,
    min_bulk_ess = min_bulk,
    min_tail_ess = min_tail,
    seconds_per_bulk_ess = if (is.finite(min_bulk) && min_bulk > 0) unname(elapsed[["elapsed"]]) / min_bulk else NA_real_,
    seconds_per_tail_ess = if (is.finite(min_tail) && min_tail > 0) unname(elapsed[["elapsed"]]) / min_tail else NA_real_,
    stringsAsFactors = FALSE
  )

  list(
    metrics = metrics,
    estimands = estimands,
    estimand_summary = key_summary,
    component_diagnostics = component_diagnostics(fit, standata, estimands),
    sampler = sampler
  )
}

summarise_values_by_chain <- function(variable, values, chains) {
  pieces <- split(values, chains)
  do.call(rbind, lapply(names(pieces), function(chain) {
    x <- pieces[[chain]]
    data.frame(
      variable = variable,
      chain = as.integer(chain),
      mean = mean(x),
      sd = stats::sd(x),
      q05 = unname(stats::quantile(x, 0.05)),
      q50 = unname(stats::quantile(x, 0.50)),
      q95 = unname(stats::quantile(x, 0.95)),
      stringsAsFactors = FALSE
    )
  }))
}

component_diagnostics <- function(fit, standata, estimands = NULL) {
  draws <- posterior::as_draws_df(fit$draws(format = "draws_df"))
  rows <- list()
  add_draw_var <- function(variable, values = NULL, chain = draws$.chain) {
    if (is.null(values)) {
      if (!variable %in% names(draws)) return(invisible(NULL))
      values <- draws[[variable]]
    }
    rows[[length(rows) + 1L]] <<- summarise_values_by_chain(variable, values, chain)
    invisible(NULL)
  }

  for (r in 1:2) {
    add_draw_var(sprintf("rho[%d]", r))
    add_draw_var(sprintf("alpha[%d]", r))
    for (h in seq_len(standata$H)) {
      add_draw_var(sprintf("m_comp[%d,%d]", r, h))
      add_draw_var(sprintf("phi_comp[%d,%d]", r, h))
      log_phi_name <- sprintf("log_phi_comp[%d,%d]", r, h)
      if (!sprintf("phi_comp[%d,%d]", r, h) %in% names(draws) && log_phi_name %in% names(draws)) {
        add_draw_var(sprintf("phi_comp[%d,%d]", r, h), exp(draws[[log_phi_name]]))
      }
      add_draw_var(sprintf("mu_logit_comp[%d,%d]", r, h))
      add_draw_var(sprintf("sigma_comp[%d,%d]", r, h))
      log_sigma_name <- sprintf("log_sigma_comp[%d,%d]", r, h)
      if (!sprintf("sigma_comp[%d,%d]", r, h) %in% names(draws) && log_sigma_name %in% names(draws)) {
        add_draw_var(sprintf("sigma_comp[%d,%d]", r, h), exp(draws[[log_sigma_name]]))
      }
    }
    if (standata$H > 1L) {
      for (h in seq_len(standata$H - 1L)) {
        v_name <- sprintf("v[%d,%d]", r, h)
        add_draw_var(v_name)
        add_draw_var(sprintf("u_stick[%d,%d]", r, h))
      }
      if (standata$H == 2L) {
        v_name <- sprintf("v[%d,1]", r)
        if (v_name %in% names(draws)) {
          add_draw_var(sprintf("w[%d,1]", r), draws[[v_name]])
          add_draw_var(sprintf("w[%d,2]", r), 1 - draws[[v_name]])
        }
        add_draw_var(sprintf("m_gap_prop[%d]", r))
      }
      for (h in seq_len(standata$H)) {
        add_draw_var(sprintf("m_logit[%d,%d]", r, h))
      }
    }
  }

  if (!is.null(estimands)) {
    estimand_draws <- posterior::as_draws_df(estimands)
    for (variable in c("mu_0_c", "mu_1_c", "mu_delta", "delta")) {
      if (variable %in% names(estimand_draws)) {
        add_draw_var(variable, estimand_draws[[variable]], estimand_draws$.chain)
      }
    }
  }

  if (!length(rows)) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

make_init_fun <- function(standata, prior, init_strategy = "empirical") {
  observed <- standata$.observed
  atom_rate <- 1 - standata$n_obs_arm / standata$n_arm
  y <- observed$Y
  kernel <- prior_value(prior, "kernel", "beta")
  component_centers <- stats::quantile(y / max(standata$score_value), probs = seq(0.35, 0.75, length.out = standata$H), names = FALSE)
  component_centers <- pmin(pmax(as.numeric(component_centers), 0.05), 0.95)
  function(chain_id = 1) {
    jitter <- function(x, amount) pmin(pmax(x + stats::runif(length(x), -amount, amount), 1e-4), 1 - 1e-4)
    jitter_positive <- function(x, amount) pmax(x * exp(stats::runif(length(x), -amount, amount)), 0.01)
    eta_init <- eta_prior_for_k(prior$eta_prior, standata$K)
    eta_init <- eta_init / sum(eta_init)
    common <- list(
      eta = matrix(rep(eta_init, standata$eta_groups), nrow = standata$eta_groups, byrow = TRUE)
    )
    if (!identical(prior_value(prior, "rho_parameterization", "sampled"), "collapsed")) {
      common$rho <- jitter(atom_rate, 0.01)
    }
    if (identical(kernel, "beta")) {
      common$log_phi_comp <- matrix(stats::rnorm(2 * standata$H, prior$phi_meanlog, 0.05), nrow = 2)
    }
    if (standata$H == 1L) {
      if (identical(kernel, "logitnormal")) {
        arm_centers <- vapply(1:2, function(r) {
          values <- observed$Y[as.integer(observed$R) + 1L == r] / max(standata$score_value)
          stats::median(values, na.rm = TRUE)
        }, numeric(1))
        arm_centers <- pmin(pmax(arm_centers, 0.05), 0.95)
        return(c(
          list(
            mu_logit_comp = matrix(stats::qlogis(jitter(arm_centers, 0.02)), nrow = 2),
            log_sigma_comp = matrix(log(jitter_positive(rep(1.2, 2), 0.05)), nrow = 2)
          ),
          common
        ))
      }
      common$m_comp <- matrix(jitter(rep(component_centers, 2), 0.03), nrow = 2)
      return(common)
    }
    if (identical(kernel, "logitnormal")) {
      if (identical(init_strategy, "control_low_main")) {
        mu_init_unit <- rbind(
          c(0.50, 0.63),
          c(0.52, 0.67)
        )
        sigma_init <- rbind(
          c(1.25, 1.00),
          c(1.20, 0.95)
        )
        v_init <- matrix(c(0.70, 0.70), nrow = 2)
      } else if (identical(init_strategy, "control_high_main")) {
        mu_init_unit <- rbind(
          c(0.48, 0.61),
          c(0.50, 0.65)
        )
        sigma_init <- rbind(
          c(1.00, 1.25),
          c(0.95, 1.20)
        )
        v_init <- matrix(c(0.30, 0.30), nrow = 2)
      } else {
        mu_init_unit <- t(replicate(2, sort(jitter(component_centers, 0.04))))
        sigma_init <- matrix(jitter_positive(rep(1.2, 2 * standata$H), 0.08), nrow = 2)
        v_init <- matrix(jitter(rep(0.55, 2 * (standata$H - 1)), 0.05), nrow = 2)
      }
      mu_init_unit <- t(apply(mu_init_unit, 1, function(x) sort(jitter(x, 0.01))))
      mu_init <- stats::qlogis(pmin(pmax(mu_init_unit, 1e-4), 1 - 1e-4))
      sigma_init <- matrix(jitter_positive(as.vector(sigma_init), 0.06), nrow = 2)
      alpha_init <- rep(prior$alpha_shape / prior$alpha_rate, 2)
      if (identical(prior_value(prior, "mu_parameterization", "ordered"), "bounded_ordered")) {
        mu_min <- -3
        mu_max <- 3
        mu_low_unit <- (mu_init[, 1] - mu_min) / (mu_max - mu_min)
        mu_gap_prop <- (mu_init[, 2] - mu_init[, 1]) / (mu_max - mu_init[, 1])
        mu_low_unit <- pmin(pmax(mu_low_unit, 1e-4), 1 - 1e-4)
        mu_gap_prop <- pmin(pmax(mu_gap_prop, 1e-3), 1 - 1e-3)
        return(c(
          list(
            v = jitter(v_init, 0.015),
            mu_low_unit = mu_low_unit,
            mu_gap_prop = mu_gap_prop,
            log_sigma_comp = log(sigma_init),
            alpha = alpha_init
          ),
          common
        ))
      }
      if (identical(prior_value(prior, "stick_parameterization", "centered"), "marginal")) {
        return(c(
          list(
            mu_logit_comp = mu_init,
            log_sigma_comp = log(sigma_init),
            alpha = alpha_init
          ),
          common
        ))
      }
      return(c(
        list(
          v = jitter(v_init, 0.015),
          mu_logit_comp = mu_init,
          log_sigma_comp = log(sigma_init),
          alpha = alpha_init
        ),
        common
      ))
    }
    if (identical(init_strategy, "control_low_main")) {
      m_init <- rbind(
        c(0.577, 0.602),
        c(0.560, 0.616)
      )
      phi_init <- rbind(
        c(5.00, 0.75),
        c(0.92, 5.10)
      )
      v_init <- matrix(c(0.85, 0.17), nrow = 2)
      m_init <- t(apply(m_init, 1, function(x) sort(jitter(x, 0.006))))
      phi_init <- matrix(jitter_positive(as.vector(phi_init), 0.06), nrow = 2)
      common$log_phi_comp <- log(phi_init)
    } else if (identical(init_strategy, "control_high_main")) {
      m_init <- rbind(
        c(0.552, 0.583),
        c(0.560, 0.616)
      )
      phi_init <- rbind(
        c(0.75, 5.00),
        c(0.92, 5.10)
      )
      v_init <- matrix(c(0.16, 0.17), nrow = 2)
      m_init <- t(apply(m_init, 1, function(x) sort(jitter(x, 0.006))))
      phi_init <- matrix(jitter_positive(as.vector(phi_init), 0.06), nrow = 2)
      common$log_phi_comp <- log(phi_init)
    } else {
      m_init <- t(replicate(2, sort(jitter(component_centers, 0.03))))
      v_init <- matrix(jitter(rep(0.55, 2 * (standata$H - 1)), 0.05), nrow = 2)
    }
    if (identical(prior_value(prior, "m_parameterization", "ordered_gap"), "ordered_logit")) {
      mean_eps <- 1e-5
      m_unit <- (pmin(pmax(m_init, 1e-4), 1 - 1e-4) - mean_eps) / (1 - 2 * mean_eps)
      m_logit <- t(apply(m_unit, 1, function(x) sort(stats::qlogis(pmin(pmax(x, 1e-4), 1 - 1e-4)))))
      alpha_init <- rep(prior$alpha_shape / prior$alpha_rate, 2)
      v_jittered <- jitter(v_init, 0.015)
      if (identical(prior_value(prior, "stick_parameterization", "centered"), "marginal")) {
        return(c(
          list(
            m_logit = m_logit,
            alpha = alpha_init
          ),
          common
        ))
      }
      if (identical(prior_value(prior, "stick_parameterization", "centered"), "noncentered")) {
        u_stick <- pmin(pmax((1 - v_jittered)^alpha_init, 1e-4), 1 - 1e-4)
        return(c(
          list(
            u_stick = u_stick,
            m_logit = m_logit,
            alpha = alpha_init
          ),
          common
        ))
      }
      return(c(
        list(
          v = v_jittered,
          m_logit = m_logit,
          alpha = alpha_init
        ),
        common
      ))
    }
    m_gap_prop <- (m_init[, 2] - m_init[, 1]) / (1 - 1e-5 - m_init[, 1])
    m_gap_prop <- pmin(pmax(m_gap_prop, 1e-4), 1 - 1e-4)
    c(
      list(
        v = jitter(v_init, 0.015),
        m_low = m_init[, 1],
        m_gap_prop = m_gap_prop,
        alpha = rep(prior$alpha_shape / prior$alpha_rate, 2)
      ),
      common
    )
  }
}

run_cmdstan_fit <- function(model, standata, prior, label, stage, adapt_delta,
                            iter_warmup, iter_sampling, chains, seed,
                            output_dir, max_treedepth = 12L, threads_per_chain = 1L,
                            init_strategy = "empirical") {
  ensure_dir(output_dir)
  started <- Sys.time()
  fit <- NULL
  elapsed <- system.time({
    available_cores <- max(1L, parallel::detectCores(logical = FALSE))
    fit <- model$sample(
      data = stan_data_only(standata),
      chains = chains,
      parallel_chains = min(chains, max(1L, floor(available_cores / threads_per_chain))),
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      seed = seed,
      refresh = 0,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      threads_per_chain = threads_per_chain,
      init = make_init_fun(standata, prior, init_strategy),
      output_dir = output_dir,
      show_messages = FALSE
    )
  })
  out <- summarize_fit(fit, standata, label, stage, adapt_delta, elapsed, iter_warmup, iter_sampling)
  out$fit <- fit
  out$started <- started
  out$finished <- Sys.time()
  out
}

diagnostic_acceptance <- function(metrics, final = FALSE) {
  if (metrics$divergences > 0) return("divergences")
  if (metrics$max_treedepth >= 12) return("treedepth")
  if (is.finite(metrics$min_bfmi) && metrics$min_bfmi < 0.3) return("bfmi")
  if (isTRUE(final)) {
    if (!is.finite(metrics$max_rhat) || metrics$max_rhat > 1.01) return("rhat")
    if (!is.finite(metrics$min_bulk_ess) || metrics$min_bulk_ess < 100) return("bulk_ess")
    if (!is.finite(metrics$min_tail_ess) || metrics$min_tail_ess < 100) return("tail_ess")
    if (metrics$elapsed_sec > 1800) return("runtime")
  }
  "accepted"
}
