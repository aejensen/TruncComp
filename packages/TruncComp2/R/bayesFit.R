bayes_trunc_comp_method_label <- function() {
  "Experimental Bayesian two-part Dirichlet process mixture model"
}

bayes_parameter_names <- function(type = c("all", "contrast", "arm")) {
  type <- match.arg(type)

  switch(
    type,
    all = c(
      "rho_0", "rho_1",
      "pi_0", "pi_1",
      "mu_0_c", "mu_1_c",
      "delta_atom", "mu_delta", "alpha_delta", "delta"
    ),
    contrast = c("delta_atom", "mu_delta", "alpha_delta", "delta"),
    arm = c("rho_0", "rho_1", "pi_0", "pi_1", "mu_0_c", "mu_1_c")
  )
}

bayes_contrast_null_value <- function(parameter) {
  if(identical(parameter, "alpha_delta")) {
    return(1)
  }

  0
}

bayes_parameter_aliases <- function(parameter) {
  parameter <- gsub("^Delta$", "delta", parameter)
  parameter <- gsub("^muDelta$", "mu_delta", parameter)
  parameter <- gsub("^alphaDelta$", "alpha_delta", parameter)
  parameter <- gsub("^deltaAtom$", "delta_atom", parameter)
  parameter
}

validate_bayes_positive_integer <- function(x, name, min_value = 1L, allow_zero = FALSE) {
  lower_bound <- if(allow_zero) 0L else as.integer(min_value)

  if(!(length(x) == 1 &&
       is.numeric(x) &&
       is.finite(x) &&
       x == as.integer(x) &&
       x >= lower_bound &&
       x >= min_value)) {
    comparator <- if(allow_zero && min_value == 0L) "a single non-negative integer" else
      paste0("a single integer >= ", min_value)
    stop(name, " must be ", comparator, ".")
  }

  as.integer(x)
}

validate_bayes_seed <- function(seed) {
  if(is.null(seed)) {
    return(NULL)
  }

  if(!(length(seed) == 1 &&
       is.numeric(seed) &&
       is.finite(seed) &&
       seed == as.integer(seed) &&
       seed >= 0)) {
    stop("seed must be NULL or a single non-negative integer.")
  }

  as.integer(seed)
}

validate_bayes_refresh <- function(refresh) {
  if(!(length(refresh) == 1 &&
       is.numeric(refresh) &&
       is.finite(refresh) &&
       refresh == as.integer(refresh) &&
       refresh >= 0)) {
    stop("refresh must be a single non-negative integer.")
  }

  as.integer(refresh)
}

validate_bayes_control <- function(control) {
  defaults <- list(adapt_delta = 0.95, max_treedepth = 12L)

  if(is.null(control)) {
    control <- list()
  }

  if(!is.list(control)) {
    stop("control must be NULL or a named list.")
  }

  unknown <- setdiff(names(control), names(defaults))
  if(length(unknown) > 0L) {
    stop("Unsupported control fields: ", paste(unknown, collapse = ", "), ".")
  }

  out <- defaults
  for(name in names(control)) {
    out[[name]] <- control[[name]]
  }

  if(!(length(out$adapt_delta) == 1 &&
       is.numeric(out$adapt_delta) &&
       is.finite(out$adapt_delta) &&
       out$adapt_delta > 0 &&
       out$adapt_delta < 1)) {
    stop("control$adapt_delta must be a single number strictly between 0 and 1.")
  }

  out$max_treedepth <- validate_bayes_positive_integer(
    out$max_treedepth,
    "control$max_treedepth",
    min_value = 1L
  )

  out
}

bayes_continuous_support <- function(continuous_support = c(
  "real_line",
  "positive_real",
  "bounded_continuous",
  "bounded_score"
)) {
  if(length(continuous_support) == 0L || is.null(continuous_support)) {
    continuous_support <- "real_line"
  }

  match.arg(
    continuous_support,
    c("real_line", "positive_real", "bounded_continuous", "bounded_score")
  )
}

bayes_bounded_kernel <- function(bounded_kernel = c("beta", "logit_normal")) {
  if(length(bounded_kernel) == 0L || is.null(bounded_kernel)) {
    bounded_kernel <- "beta"
  }

  match.arg(bounded_kernel, c("beta", "logit_normal"))
}

bayes_normalize_bounded_kernel <- function(continuous_support,
                                           bounded_kernel = c("beta", "logit_normal"),
                                           supplied = FALSE) {
  continuous_support <- bayes_continuous_support(continuous_support)
  bounded_kernel <- bayes_bounded_kernel(bounded_kernel)

  if(!bayes_is_bounded_support(continuous_support) && isTRUE(supplied)) {
    stop("bounded_kernel is only supported for bounded Bayesian models.")
  }

  bounded_kernel
}

bayes_model_name <- function(continuous_support = c(
  "real_line",
  "positive_real",
  "bounded_continuous",
  "bounded_score"
),
                             bounded_kernel = c("beta", "logit_normal")) {
  continuous_support <- bayes_continuous_support(continuous_support)
  bounded_kernel <- bayes_normalize_bounded_kernel(
    continuous_support,
    bounded_kernel = bounded_kernel,
    supplied = FALSE
  )

  switch(
    continuous_support,
    real_line = "trunc_comp_bayes",
    positive_real = "trunc_comp_bayes_positive",
    bounded_continuous = if(identical(bounded_kernel, "logit_normal")) {
      "trunc_comp_bayes_bounded_continuous_logit_normal"
    } else {
      "trunc_comp_bayes_bounded_continuous"
    },
    bounded_score = if(identical(bounded_kernel, "logit_normal")) {
      "trunc_comp_bayes_bounded_score_logit_normal"
    } else {
      "trunc_comp_bayes_bounded_score"
    }
  )
}

bayes_default_prior <- function(continuous_support = c(
  "real_line",
  "positive_real",
  "bounded_continuous",
  "bounded_score"
),
                                bounded_kernel = c("beta", "logit_normal")) {
  continuous_support <- bayes_continuous_support(continuous_support)
  bounded_kernel <- bayes_normalize_bounded_kernel(
    continuous_support,
    bounded_kernel = bounded_kernel,
    supplied = FALSE
  )

  common <- list(
    rho_alpha = 1,
    rho_beta = 1,
    alpha_shape = 2,
    alpha_rate = 1
  )

  support_specific <- switch(
    continuous_support,
    real_line = list(
      mu_mean = 0,
      mu_sd = 2.5,
      sigma_meanlog = -0.5,
      sigma_sdlog = 0.5
    ),
    positive_real = list(
      mean_meanlog = 0,
      mean_sdlog = 0.5,
      shape_meanlog = log(2),
      shape_sdlog = 0.5
    ),
    bounded_continuous = if(identical(bounded_kernel, "logit_normal")) {
      list(
        mu_logit_mean = 0,
        mu_logit_sd = 2,
        sigma_logit_meanlog = log(1.2),
        sigma_logit_sdlog = 0.5
      )
    } else {
      list(
        m_alpha = 1,
        m_beta = 1,
        phi_meanlog = log(20),
        phi_sdlog = 1
      )
    },
    bounded_score = c(
      if(identical(bounded_kernel, "logit_normal")) {
        list(
          mu_logit_mean = 0,
          mu_logit_sd = 2,
          sigma_logit_meanlog = log(1.2),
          sigma_logit_sdlog = 0.5
        )
      } else {
        list(
          m_alpha = 1,
          m_beta = 1,
          phi_meanlog = log(20),
          phi_sdlog = 1
        )
      },
      list(eta_prior = 1)
    )
  )

  c(common, support_specific)
}

normalize_bayes_prior <- function(prior,
                                  continuous_support = c(
                                    "real_line",
                                    "positive_real",
                                    "bounded_continuous",
                                    "bounded_score"
                                  ),
                                  bounded_kernel = c("beta", "logit_normal")) {
  continuous_support <- bayes_continuous_support(continuous_support)
  bounded_kernel <- bayes_normalize_bounded_kernel(
    continuous_support,
    bounded_kernel = bounded_kernel,
    supplied = FALSE
  )
  defaults <- bayes_default_prior(continuous_support, bounded_kernel = bounded_kernel)

  if(is.null(prior)) {
    return(defaults)
  }

  if(!is.list(prior)) {
    stop("prior must be NULL or a named list.")
  }

  unknown <- setdiff(names(prior), names(defaults))
  if(length(unknown) > 0L) {
    stop("Unsupported prior fields: ", paste(unknown, collapse = ", "), ".")
  }

  out <- defaults
  for(name in names(prior)) {
    out[[name]] <- prior[[name]]
  }

  positive_names <- c("rho_alpha", "rho_beta", "alpha_shape", "alpha_rate")
  finite_names <- character(0)

  if(identical(continuous_support, "real_line")) {
    positive_names <- c(positive_names, "mu_sd", "sigma_sdlog")
    finite_names <- c("mu_mean", "sigma_meanlog")
  } else if(identical(continuous_support, "positive_real")) {
    positive_names <- c(positive_names, "mean_sdlog", "shape_sdlog")
    finite_names <- c("mean_meanlog", "shape_meanlog")
  } else if(identical(bounded_kernel, "logit_normal")) {
    positive_names <- c(positive_names, "mu_logit_sd", "sigma_logit_sdlog")
    finite_names <- c("mu_logit_mean", "sigma_logit_meanlog")
  } else {
    positive_names <- c(positive_names, "m_alpha", "m_beta", "phi_sdlog")
    finite_names <- c("phi_meanlog")
  }

  for(name in positive_names) {
    value <- out[[name]]
    if(!(length(value) == 1 && is.numeric(value) && is.finite(value) && value > 0)) {
      stop("prior$", name, " must be a single positive finite numeric value.")
    }
  }

  for(name in finite_names) {
    value <- out[[name]]
    if(!(length(value) == 1 && is.numeric(value) && is.finite(value))) {
      stop("prior$", name, " must be a single finite numeric value.")
    }
  }

  if(identical(continuous_support, "bounded_score")) {
    value <- out$eta_prior
    if(!(is.numeric(value) &&
         length(value) >= 1L &&
         all(is.finite(value)) &&
         all(value > 0))) {
      stop("prior$eta_prior must contain positive finite numeric values.")
    }
    out$eta_prior <- as.numeric(value)
  }

  out
}

normalize_bayes_sampling_args <- function(extra_args) {
  if(length(extra_args) == 0L) {
    return(extra_args)
  }

  argument_names <- names(extra_args)
  if(is.null(argument_names)) {
    argument_names <- rep("", length(extra_args))
  }

  if("adjust" %in% argument_names) {
    stop("Covariate adjustment is not implemented for trunc_comp_bayes().")
  }

  reserved <- c("object", "data", "chains", "iter", "warmup", "seed", "refresh", "control")
  conflicts <- intersect(argument_names[nzchar(argument_names)], reserved)
  if(length(conflicts) > 0L) {
    stop("Do not pass reserved Stan arguments through .... Conflicts: ",
         paste(conflicts, collapse = ", "), ".")
  }

  extra_args
}

bayes_outcome_center_scale <- function(y,
                                       continuous_support = c(
                                         "real_line",
                                         "positive_real",
                                         "bounded_continuous",
                                         "bounded_score"
                                       )) {
  continuous_support <- bayes_continuous_support(continuous_support)

  if(identical(continuous_support, "positive_real")) {
    center <- 0
    scale <- mean(y)

    if(!(is.finite(scale) && scale > 0)) {
      scale <- stats::median(y)
    }

    if(!(is.finite(scale) && scale > 0)) {
      scale <- stats::sd(y)
    }

    if(!(is.finite(scale) && scale > 0)) {
      scale <- 1
    }

    return(list(center = 0, scale = as.numeric(scale)))
  }

  center <- mean(y)
  scale <- stats::sd(y)

  if(!(is.finite(scale) && scale > 0)) {
    scale <- stats::IQR(y) / 1.349
  }

  if(!(is.finite(scale) && scale > 0)) {
    scale <- max(abs(y - center))
  }

  if(!(is.finite(scale) && scale > 0)) {
    scale <- 1
  }

  list(center = as.numeric(center), scale = as.numeric(scale))
}

validate_bayes_support_data <- function(data,
                                        atom = NULL,
                                        continuous_support = c(
                                          "real_line",
                                          "positive_real",
                                          "bounded_continuous",
                                          "bounded_score"
                                        ),
                                        support_options = NULL) {
  continuous_support <- bayes_continuous_support(continuous_support)

  if(identical(continuous_support, "real_line")) {
    return(invisible(TRUE))
  }

  if(bayes_is_bounded_support(continuous_support)) {
    if(is.null(support_options)) {
      stop("Bounded Bayesian support settings are missing.")
    }
    bayes_validate_bounded_data(data = data, atom = atom, support_options = support_options)
    return(invisible(TRUE))
  }

  observed <- data$Y[data$A == 1]
  if(any(observed <= 0)) {
    stop(
      "For continuous_support = \"positive_real\", all observed non-atom outcomes must be strictly positive."
    )
  }

  invisible(TRUE)
}

validate_bayes_flag <- function(x, name) {
  if(!(length(x) == 1L && is.logical(x) && !is.na(x))) {
    stop(name, " must be TRUE or FALSE.")
  }

  isTRUE(x)
}

normalize_bayes_mixture_components_max <- function(mixture_components_max,
                                                   mixture_components,
                                                   auto_select_mixture_components = TRUE) {
  auto_select_mixture_components <- validate_bayes_flag(
    auto_select_mixture_components,
    "auto_select_mixture_components"
  )

  if(is.null(mixture_components_max)) {
    if(auto_select_mixture_components) {
      return(as.integer(max(40L, mixture_components)))
    }

    return(as.integer(mixture_components))
  }

  mixture_components_max <- validate_bayes_positive_integer(
    mixture_components_max,
    "mixture_components_max",
    min_value = mixture_components
  )

  if(mixture_components_max < mixture_components) {
    stop("mixture_components_max must be >= mixture_components.")
  }

  mixture_components_max
}

bayes_mixture_component_ladder <- function(mixture_components,
                                           mixture_components_max = mixture_components) {
  mixture_components <- validate_bayes_positive_integer(
    mixture_components,
    "mixture_components",
    min_value = 2L
  )
  mixture_components_max <- validate_bayes_positive_integer(
    mixture_components_max,
    "mixture_components_max",
    min_value = mixture_components
  )

  ladder <- integer(0)
  current <- mixture_components

  while(current <= mixture_components_max) {
    ladder <- c(ladder, current)
    current <- current * 2L
  }

  as.integer(unique(ladder))
}

build_bayes_standata <- function(data, atom, mixture_components, prior,
                                 continuous_support = c(
                                   "real_line",
                                   "positive_real",
                                   "bounded_continuous",
                                   "bounded_score"
                                 ),
                                 bounded_kernel = c("beta", "logit_normal"),
                                 support_options = NULL) {
  continuous_support <- bayes_continuous_support(continuous_support)
  bounded_kernel <- bayes_normalize_bounded_kernel(
    continuous_support,
    bounded_kernel = bounded_kernel,
    supplied = FALSE
  )
  observed <- droplevels(data[data$A == 1, c("Y", "R"), drop = FALSE])
  scaling <- if(bayes_is_bounded_support(continuous_support)) {
    if(is.null(support_options)) {
      stop("Bounded Bayesian support settings are missing.")
    }

    list(
      center = as.numeric(support_options$score_min),
      scale = as.numeric(support_options$score_range)
    )
  } else {
    bayes_outcome_center_scale(observed$Y, continuous_support = continuous_support)
  }

  common_data <- list(
    N = nrow(data),
    A = as.integer(data$A),
    arm = as.integer(data$R) + 1L,
    H = as.integer(mixture_components),
    N_obs = nrow(observed),
    arm_obs = as.integer(observed$R) + 1L,
    atom = as.numeric(atom),
    rho_prior_alpha = as.numeric(prior$rho_alpha),
    rho_prior_beta = as.numeric(prior$rho_beta),
    alpha_prior_shape = as.numeric(prior$alpha_shape),
    alpha_prior_rate = as.numeric(prior$alpha_rate)
  )

  if(identical(continuous_support, "real_line")) {
    y_obs_std <- (observed$Y - scaling$center) / scaling$scale

    return(list(
      stan_data = c(
        common_data,
        list(
          y_obs_std = as.numeric(y_obs_std),
          y_center = scaling$center,
          y_scale = scaling$scale,
          mu_prior_mean = as.numeric(prior$mu_mean),
          mu_prior_sd = as.numeric(prior$mu_sd),
          sigma_prior_meanlog = as.numeric(prior$sigma_meanlog),
          sigma_prior_sdlog = as.numeric(prior$sigma_sdlog)
        )
      ),
      scaling = scaling,
      model_name = bayes_model_name(continuous_support)
    ))
  }

  if(identical(continuous_support, "bounded_continuous")) {
    x_obs <- (observed$Y - support_options$score_min) / support_options$score_range
    kernel_prior <- if(identical(bounded_kernel, "logit_normal")) {
      quad <- bayes_logitnormal_mean_quad()
      list(
        mu_logit_prior_mean = as.numeric(prior$mu_logit_mean),
        mu_logit_prior_sd = as.numeric(prior$mu_logit_sd),
        sigma_logit_prior_meanlog = as.numeric(prior$sigma_logit_meanlog),
        sigma_logit_prior_sdlog = as.numeric(prior$sigma_logit_sdlog),
        mean_quad_n = as.integer(length(quad$nodes)),
        mean_quad_node = as.numeric(quad$nodes),
        mean_quad_weight = as.numeric(quad$weights)
      )
    } else {
      list(
        m_prior_alpha = as.numeric(prior$m_alpha),
        m_prior_beta = as.numeric(prior$m_beta),
        phi_prior_meanlog = as.numeric(prior$phi_meanlog),
        phi_prior_sdlog = as.numeric(prior$phi_sdlog)
      )
    }

    return(list(
      stan_data = c(
        common_data,
        list(
          x_obs = as.numeric(x_obs),
          score_min = as.numeric(support_options$score_min),
          score_max = as.numeric(support_options$score_max)
        ),
        kernel_prior
      ),
      scaling = scaling,
      model_name = bayes_model_name(continuous_support, bounded_kernel = bounded_kernel)
    ))
  }

  if(identical(continuous_support, "bounded_score")) {
    y_obs_index <- bayes_score_grid_index(
      observed$Y,
      score_min = support_options$score_min,
      score_max = support_options$score_max,
      score_step = support_options$score_step
    )
    score_cells <- bayes_score_cell_counts(
      arm_obs = as.integer(observed$R) + 1L,
      y_obs_index = y_obs_index,
      j_count = length(support_options$score_values)
    )
    eta_prior <- bayes_validate_eta_prior(
      prior$eta_prior,
      k = length(support_options$heaping_grids)
    )
    n_arm <- tabulate(as.integer(data$R) + 1L, nbins = 2L)
    n_obs_arm <- tabulate(as.integer(observed$R) + 1L, nbins = 2L)
    bounded_score_common <- common_data[c(
      "H",
      "atom",
      "rho_prior_alpha",
      "rho_prior_beta",
      "alpha_prior_shape",
      "alpha_prior_rate"
    )]
    kernel_prior <- if(identical(bounded_kernel, "logit_normal")) {
      list(
        mu_logit_prior_mean = as.numeric(prior$mu_logit_mean),
        mu_logit_prior_sd = as.numeric(prior$mu_logit_sd),
        sigma_logit_prior_meanlog = as.numeric(prior$sigma_logit_meanlog),
        sigma_logit_prior_sdlog = as.numeric(prior$sigma_logit_sdlog)
      )
    } else {
      list(
        m_prior_alpha = as.numeric(prior$m_alpha),
        m_prior_beta = as.numeric(prior$m_beta),
        phi_prior_meanlog = as.numeric(prior$phi_meanlog),
        phi_prior_sdlog = as.numeric(prior$phi_sdlog)
      )
    }

    return(list(
      stan_data = c(
        bounded_score_common,
        score_cells,
        list(
          n_arm = as.integer(n_arm),
          n_obs_arm = as.integer(n_obs_arm),
          J = as.integer(length(support_options$score_values)),
          score_value = as.numeric(support_options$score_values),
          K = as.integer(length(support_options$heaping_grids)),
          bin_lower = support_options$bin_lower,
          bin_upper = support_options$bin_upper,
          bin_valid = support_options$bin_valid,
          eta_groups = as.integer(support_options$eta_groups),
          eta_group_by_arm = as.integer(support_options$eta_group_by_arm),
          eta_prior = eta_prior
        ),
        kernel_prior
      ),
      scaling = scaling,
      model_name = bayes_model_name(continuous_support, bounded_kernel = bounded_kernel)
    ))
  }

  y_obs_scaled <- observed$Y / scaling$scale

  list(
    stan_data = c(
      common_data,
      list(
        y_obs_scaled = as.numeric(y_obs_scaled),
        y_scale = scaling$scale,
        mean_prior_meanlog = as.numeric(prior$mean_meanlog),
        mean_prior_sdlog = as.numeric(prior$mean_sdlog),
        shape_prior_meanlog = as.numeric(prior$shape_meanlog),
        shape_prior_sdlog = as.numeric(prior$shape_sdlog)
      )
    ),
    scaling = scaling,
    model_name = bayes_model_name(continuous_support)
  )
}

bayes_extract_draws <- function(fit) {
  draws_array <- rstan::extract(
    fit,
    pars = bayes_parameter_names("all"),
    permuted = FALSE,
    inc_warmup = FALSE
  )

  draws <- posterior::as_draws_df(draws_array)
  required <- bayes_parameter_names("all")
  missing <- setdiff(required, names(draws))

  if(length(missing) > 0L) {
    stop("Posterior draws are missing required variables: ",
         paste(missing, collapse = ", "), ".")
  }

  draws
}

bayes_equal_tail_interval <- function(x, conf.level) {
  stats::quantile(
    x,
    probs = c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2),
    names = FALSE
  )
}

bayes_summary_table <- function(draws, conf.level) {
  vars <- bayes_parameter_names("contrast")

  rows <- lapply(vars, function(var) {
    values <- draws[[var]]
    interval <- bayes_equal_tail_interval(values, conf.level)
    data.frame(
      parameter = var,
      estimate = stats::median(values),
      conf.low = interval[[1]],
      conf.high = interval[[2]],
      posterior_prob = mean(values > bayes_contrast_null_value(var)),
      row.names = NULL
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- out$parameter
  out
}

bayes_arm_table <- function(draws, conf.level) {
  vars <- bayes_parameter_names("arm")

  rows <- lapply(vars, function(var) {
    values <- draws[[var]]
    interval <- bayes_equal_tail_interval(values, conf.level)
    data.frame(
      parameter = var,
      estimate = stats::median(values),
      conf.low = interval[[1]],
      conf.high = interval[[2]],
      row.names = NULL
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- out$parameter
  out
}

bayes_select_draws <- function(draws, variables) {
  metadata <- intersect(c(".chain", ".iteration", ".draw"), names(draws))
  draws[, c(metadata, variables), drop = FALSE]
}

bayes_truncation_parameter_names <- function(mixture_components) {
  mixture_components <- validate_bayes_positive_integer(
    mixture_components,
    "mixture_components",
    min_value = 2L
  )

  c(
    "alpha[1]",
    "alpha[2]",
    sprintf("w[1,%d]", mixture_components),
    sprintf("w[2,%d]", mixture_components)
  )
}

bayes_extract_truncation_draws <- function(fit, mixture_components) {
  mixture_components <- validate_bayes_positive_integer(
    mixture_components,
    "mixture_components",
    min_value = 2L
  )

  draws_array <- rstan::extract(
    fit,
    pars = c("alpha", "w"),
    permuted = FALSE,
    inc_warmup = FALSE
  )

  draws <- posterior::as_draws_df(draws_array)
  required <- bayes_truncation_parameter_names(mixture_components)
  missing <- setdiff(required, names(draws))

  if(length(missing) > 0L) {
    stop("Posterior draws are missing required truncation variables: ",
         paste(missing, collapse = ", "), ".")
  }

  bayes_select_draws(draws, required)
}

bayes_convergence_ok <- function(max_rhat, min_bulk_ess, min_tail_ess) {
  is.finite(max_rhat) &&
    max_rhat <= 1.01 &&
    is.finite(min_bulk_ess) &&
    min_bulk_ess >= 400 &&
    is.finite(min_tail_ess) &&
    min_tail_ess >= 400
}

bayes_convergence_summary <- function(draws, parameters) {
  diagnostic_draws <- bayes_select_draws(draws, parameters)
  convergence <- posterior::summarise_draws(
    diagnostic_draws,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  convergence <- as.data.frame(convergence)
  rownames(convergence) <- convergence$variable

  max_rhat <- max(convergence$rhat, na.rm = TRUE)
  min_bulk_ess <- min(convergence$ess_bulk, na.rm = TRUE)
  min_tail_ess <- min(convergence$ess_tail, na.rm = TRUE)

  list(
    max_rhat = as.numeric(max_rhat),
    min_bulk_ess = as.numeric(min_bulk_ess),
    min_tail_ess = as.numeric(min_tail_ess),
    convergence_ok = bayes_convergence_ok(max_rhat, min_bulk_ess, min_tail_ess),
    parameter_table = convergence
  )
}

bayes_truncation_ok <- function(convergence_ok, divergences = 0L) {
  isTRUE(divergences == 0L) &&
    isTRUE(convergence_ok)
}

bayes_truncation_diagnostics <- function(truncation_draws, divergences = 0L) {
  parameters <- setdiff(names(truncation_draws), c(".chain", ".iteration", ".draw"))
  truncation_summary <- bayes_convergence_summary(
    truncation_draws,
    parameters
  )

  list(
    parameter_table = truncation_summary$parameter_table,
    max_rhat = truncation_summary$max_rhat,
    min_bulk_ess = truncation_summary$min_bulk_ess,
    min_tail_ess = truncation_summary$min_tail_ess,
    convergence_ok = truncation_summary$convergence_ok,
    truncation_ok = bayes_truncation_ok(
      convergence_ok = truncation_summary$convergence_ok,
      divergences = divergences
    )
  )
}

bayes_sampler_divergences <- function(fit) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  as.integer(sum(vapply(
    sampler_params,
    function(x) sum(x[, "divergent__"]),
    numeric(1)
  )))
}

bayes_diagnostics <- function(fit, draws, truncation_draws,
                              parameters = bayes_parameter_names("contrast")) {
  divergences <- bayes_sampler_divergences(fit)
  core_summary <- bayes_convergence_summary(draws, parameters)
  truncation_summary <- bayes_truncation_diagnostics(
    truncation_draws = truncation_draws,
    divergences = divergences
  )
  core_ok <- isTRUE(divergences == 0L) && isTRUE(core_summary$convergence_ok)
  truncation_ok <- isTRUE(truncation_summary$truncation_ok)

  list(
    divergences = divergences,
    max_rhat = core_summary$max_rhat,
    min_bulk_ess = core_summary$min_bulk_ess,
    min_tail_ess = core_summary$min_tail_ess,
    parameter_table = core_summary$parameter_table,
    core_parameter_table = core_summary$parameter_table,
    core_ok = core_ok,
    truncation = truncation_summary,
    truncation_ok = truncation_ok,
    diagnostic_ok = isTRUE(core_ok) && isTRUE(truncation_ok)
  )
}

bayes_mixture_selection_history_row <- function(mixture_components,
                                                diagnostics = NULL,
                                                accepted = FALSE,
                                                error = NA_character_) {
  if(is.null(diagnostics)) {
    return(data.frame(
      mixture_components = as.integer(mixture_components),
      core_ok = NA,
      truncation_ok = FALSE,
      diagnostic_ok = FALSE,
      accepted = isTRUE(accepted),
      divergences = NA_integer_,
      max_rhat_core = NA_real_,
      min_bulk_ess_core = NA_real_,
      min_tail_ess_core = NA_real_,
      max_rhat_truncation = NA_real_,
      min_bulk_ess_truncation = NA_real_,
      min_tail_ess_truncation = NA_real_,
      fit_failed = TRUE,
      error = as.character(error),
      row.names = NULL
    ))
  }

  data.frame(
    mixture_components = as.integer(mixture_components),
    core_ok = isTRUE(diagnostics$core_ok),
    truncation_ok = isTRUE(diagnostics$truncation_ok),
    diagnostic_ok = isTRUE(diagnostics$diagnostic_ok),
    accepted = isTRUE(accepted),
    divergences = as.integer(diagnostics$divergences),
    max_rhat_core = as.numeric(diagnostics$max_rhat),
    min_bulk_ess_core = as.numeric(diagnostics$min_bulk_ess),
    min_tail_ess_core = as.numeric(diagnostics$min_tail_ess),
    max_rhat_truncation = as.numeric(diagnostics$truncation$max_rhat),
    min_bulk_ess_truncation = as.numeric(diagnostics$truncation$min_bulk_ess),
    min_tail_ess_truncation = as.numeric(diagnostics$truncation$min_tail_ess),
    fit_failed = FALSE,
    error = as.character(error),
    row.names = NULL
  )
}

new_trunc_comp_bayes_fit <- function(fit = NULL, draws = NULL,
                                     summary_table = NULL, arm_table = NULL,
                                     diagnostics = NULL, settings = NULL,
                                     ppc_table = NULL, ppc_settings = NULL,
                                     conf.level, success, error = "",
                                     data = NULL, atom = NULL, call = NULL) {
  out <- list(
    fit = fit,
    draws = draws,
    summary_table = summary_table,
    arm_table = arm_table,
    diagnostics = diagnostics,
    settings = settings,
    ppc_table = ppc_table,
    ppc_settings = ppc_settings,
    conf.level = conf.level,
    success = success,
    error = error,
    data = data,
    atom = atom,
    call = call
  )

  class(out) <- c("trunc_comp_bayes_fit", "list")
  out
}

new_failed_trunc_comp_bayes_fit <- function(error, conf.level,
                                            data = NULL, atom = NULL,
                                            call = NULL, settings = NULL) {
  new_trunc_comp_bayes_fit(
    conf.level = conf.level,
    success = FALSE,
    error = error,
    data = data,
    atom = atom,
    call = call,
    settings = settings
  )
}

bayes_package_stanmodel <- function(model_name) {
  if(!exists("stanmodels", inherits = TRUE)) {
    return(NULL)
  }

  package_stanmodels <- get("stanmodels", inherits = TRUE)
  package_stanmodels[[model_name]]
}

bayes_finalize_trunc_comp_bayes_fit <- function(fit_object) {
  default_ppc <- tryCatch(
    bayes_ppc_summary(fit_object),
    error = identity
  )

  if(!inherits(default_ppc, "error")) {
    fit_object$ppc_table <- default_ppc$table
    fit_object$ppc_settings <- default_ppc$settings
  }

  fit_object
}

fit_trunc_comp_bayes_once <- function(data, atom, conf.level = 0.95,
                                      continuous_support = c(
                                        "real_line",
                                        "positive_real",
                                        "bounded_continuous",
                                        "bounded_score"
                                      ),
                                      bounded_kernel = c("beta", "logit_normal"),
                                      mixture_components = 10,
                                      chains = 4, iter_warmup = 1000,
                                      iter_sampling = 1000, seed = NULL,
                                      refresh = 0, control = list(adapt_delta = 0.95, max_treedepth = 12),
                                      prior = NULL, call = NULL,
                                      extra_args = list(),
                                      support_options = NULL,
                                      model_object = NULL,
                                      settings = NULL) {
  standata <- build_bayes_standata(
    data,
    atom,
    mixture_components,
    prior,
    continuous_support = continuous_support,
    bounded_kernel = bounded_kernel,
    support_options = support_options
  )

  settings$mixture_components <- as.integer(mixture_components)
  settings$y_center <- standata$scaling$center
  settings$y_scale <- standata$scaling$scale

  sampling_args <- c(
    list(
      object = model_object,
      data = standata$stan_data,
      chains = chains,
      iter = iter_warmup + iter_sampling,
      warmup = iter_warmup,
      seed = seed,
      refresh = refresh,
      control = control
    ),
    extra_args
  )

  fit <- tryCatch(
    do.call(rstan::sampling, sampling_args),
    error = identity
  )

  if(inherits(fit, "error")) {
    return(new_failed_trunc_comp_bayes_fit(
      fit$message,
      conf.level = conf.level,
      data = data,
      atom = atom,
      call = call,
      settings = settings
    ))
  }

  draws <- tryCatch(bayes_extract_draws(fit), error = identity)
  if(inherits(draws, "error") || nrow(draws) == 0L) {
    return(new_failed_trunc_comp_bayes_fit(
      "Stan sampling completed but no posterior draws could be extracted.",
      conf.level = conf.level,
      data = data,
      atom = atom,
      call = call,
      settings = settings
    ))
  }

  truncation_draws <- tryCatch(
    bayes_extract_truncation_draws(fit, mixture_components = mixture_components),
    error = identity
  )
  if(inherits(truncation_draws, "error") || nrow(truncation_draws) == 0L) {
    return(new_failed_trunc_comp_bayes_fit(
      "Stan sampling completed but the truncation diagnostics could not be extracted.",
      conf.level = conf.level,
      data = data,
      atom = atom,
      call = call,
      settings = settings
    ))
  }

  summary_table <- bayes_summary_table(draws, conf.level)
  arm_table <- bayes_arm_table(draws, conf.level)
  diagnostics <- bayes_diagnostics(
    fit = fit,
    draws = draws,
    truncation_draws = truncation_draws
  )

  new_trunc_comp_bayes_fit(
    fit = fit,
    draws = draws,
    summary_table = summary_table,
    arm_table = arm_table,
    diagnostics = diagnostics,
    settings = settings,
    conf.level = conf.level,
    success = TRUE,
    data = data,
    atom = atom,
    call = call
  )
}

fit_trunc_comp_bayes <- function(data, atom, conf.level = 0.95,
                                 continuous_support = c(
                                   "real_line",
                                   "positive_real",
                                   "bounded_continuous",
                                   "bounded_score"
                                 ),
                                 bounded_kernel = c("beta", "logit_normal"),
                                 bounded_kernel_supplied = FALSE,
                                 mixture_components = 10,
                                 auto_select_mixture_components = TRUE,
                                 mixture_components_max = NULL,
                                 chains = 4, iter_warmup = 1000,
                                 iter_sampling = 1000, seed = NULL,
                                 refresh = 0, control = list(adapt_delta = 0.95, max_treedepth = 12),
                                 prior = NULL, call = NULL,
                                 score_min = NULL,
                                 score_max = NULL,
                                 score_step = 1,
                                 heaping_grids = 1,
                                 heaping = "shared",
                                 support_supplied = bayes_support_supplied_defaults(),
                                 extra_args = list()) {
  conf.level <- validateConfidenceLevel(conf.level)
  continuous_support <- bayes_continuous_support(continuous_support)
  bounded_kernel <- bayes_normalize_bounded_kernel(
    continuous_support,
    bounded_kernel = bounded_kernel,
    supplied = bounded_kernel_supplied
  )
  support_options <- bayes_normalize_support_options(
    continuous_support = continuous_support,
    score_min = score_min,
    score_max = score_max,
    score_step = score_step,
    heaping_grids = heaping_grids,
    heaping = heaping,
    supplied = support_supplied
  )
  mixture_components <- validate_bayes_positive_integer(
    mixture_components,
    "mixture_components",
    min_value = 2L
  )
  auto_select_mixture_components <- validate_bayes_flag(
    auto_select_mixture_components,
    "auto_select_mixture_components"
  )
  mixture_components_max <- normalize_bayes_mixture_components_max(
    mixture_components_max = mixture_components_max,
    mixture_components = mixture_components,
    auto_select_mixture_components = auto_select_mixture_components
  )
  chains <- validate_bayes_positive_integer(chains, "chains", min_value = 1L)
  iter_warmup <- validate_bayes_positive_integer(iter_warmup, "iter_warmup", min_value = 1L)
  iter_sampling <- validate_bayes_positive_integer(iter_sampling, "iter_sampling", min_value = 1L)
  seed <- validate_bayes_seed(seed)
  refresh <- validate_bayes_refresh(refresh)
  control <- validate_bayes_control(control)
  prior <- normalize_bayes_prior(
    prior,
    continuous_support = continuous_support,
    bounded_kernel = bounded_kernel
  )
  if(identical(continuous_support, "bounded_score")) {
    prior$eta_prior <- bayes_validate_eta_prior(
      prior$eta_prior,
      k = length(support_options$heaping_grids)
    )
  }
  extra_args <- normalize_bayes_sampling_args(extra_args)
  model_name <- bayes_model_name(continuous_support, bounded_kernel = bounded_kernel)

  settings <- list(
    continuous_support = continuous_support,
    bounded_kernel = bounded_kernel,
    model_name = model_name,
    mixture_components = mixture_components,
    auto_select_mixture_components = auto_select_mixture_components,
    mixture_components_initial = mixture_components,
    mixture_components_final = mixture_components,
    mixture_components_max = mixture_components_max,
    mixture_component_path = bayes_mixture_component_ladder(
      mixture_components = mixture_components,
      mixture_components_max = mixture_components_max
    ),
    mixture_selection_history = NULL,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = refresh,
    control = control,
    prior = prior,
    experimental = TRUE
  )
  settings <- c(
    settings,
    support_options[setdiff(names(support_options), "continuous_support")]
  )

  if(!isDataOkay(data)) {
    return(new_failed_trunc_comp_bayes_fit(
      "Estimation failed due to data error.",
      conf.level = conf.level,
      data = data,
      atom = atom,
      call = call,
      settings = settings
    ))
  }

  validate_bayes_support_data(
    data,
    atom = atom,
    continuous_support = continuous_support,
    support_options = support_options
  )

  model_object <- bayes_package_stanmodel(model_name)
  if(is.null(model_object)) {
    return(new_failed_trunc_comp_bayes_fit(
      paste0(
        "The packaged Stan model '",
        model_name,
        "' is not available. Reinstall TruncComp2 after generating Stan exports."
      ),
      conf.level = conf.level,
      data = data,
      atom = atom,
      call = call,
      settings = settings
    ))
  }

  candidate_path <- if(auto_select_mixture_components) {
    settings$mixture_component_path
  } else {
    as.integer(mixture_components)
  }
  history_rows <- vector("list", length(candidate_path))
  successful_fits <- list()
  selected_fit <- NULL
  selection_note <- NULL

  for(index in seq_along(candidate_path)) {
    candidate_H <- candidate_path[[index]]
    candidate_settings <- settings
    candidate_settings$mixture_components <- as.integer(candidate_H)

    candidate_fit <- fit_trunc_comp_bayes_once(
      data = data,
      atom = atom,
      conf.level = conf.level,
      continuous_support = continuous_support,
      bounded_kernel = bounded_kernel,
      mixture_components = candidate_H,
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      seed = seed,
      refresh = refresh,
      control = control,
      prior = prior,
      call = call,
      extra_args = extra_args,
      support_options = support_options,
      model_object = model_object,
      settings = candidate_settings
    )

    if(!isTRUE(candidate_fit$success)) {
      history_rows[[index]] <- bayes_mixture_selection_history_row(
        mixture_components = candidate_H,
        diagnostics = NULL,
        accepted = FALSE,
        error = candidate_fit$error
      )

      if(index == 1L) {
        candidate_fit$settings$mixture_selection_history <- do.call(
          rbind,
          history_rows[seq_len(index)]
        )
        candidate_fit$settings$mixture_component_path <- as.integer(
          candidate_fit$settings$mixture_selection_history$mixture_components
        )
        candidate_fit$settings$mixture_selection_note <- if(auto_select_mixture_components) {
          paste0(
            "The initial Bayesian fit at H = ",
            candidate_H,
            " failed before automatic truncation selection could proceed."
          )
        } else {
          paste0(
            "The Bayesian fit at H = ",
            candidate_H,
            " failed."
          )
        }
        return(candidate_fit)
      }

      selected_fit <- successful_fits[[length(successful_fits)]]
      selection_note <- paste0(
        "Automatic mixture-component selection stopped after a failed refit at H = ",
        candidate_H,
        ". Returning the last successful fit at H = ",
        selected_fit$settings$mixture_components,
        "."
      )
      selected_fit$diagnostics$truncation_ok <- FALSE
      selected_fit$diagnostics$diagnostic_ok <- FALSE
      selected_fit$diagnostics$truncation$truncation_ok <- FALSE
      break
    }

    successful_fits[[length(successful_fits) + 1L]] <- candidate_fit
    accepted <- isTRUE(candidate_fit$diagnostics$diagnostic_ok)
    history_rows[[index]] <- bayes_mixture_selection_history_row(
      mixture_components = candidate_H,
      diagnostics = candidate_fit$diagnostics,
      accepted = accepted
    )

    if(!auto_select_mixture_components || accepted) {
      selected_fit <- candidate_fit
      break
    }
  }

  if(is.null(selected_fit)) {
    selected_fit <- successful_fits[[length(successful_fits)]]
    selection_note <- paste0(
      "No candidate H up to ",
      mixture_components_max,
      " satisfied the truncation diagnostics. Returning the largest successful fit at H = ",
      selected_fit$settings$mixture_components,
      "."
    )
    selected_fit$diagnostics$truncation_ok <- FALSE
    selected_fit$diagnostics$diagnostic_ok <- FALSE
    selected_fit$diagnostics$truncation$truncation_ok <- FALSE
  }

  selection_history <- do.call(
    rbind,
    history_rows[!vapply(history_rows, is.null, logical(1))]
  )

  selected_fit$settings$auto_select_mixture_components <- auto_select_mixture_components
  selected_fit$settings$mixture_components <- as.integer(selected_fit$settings$mixture_components)
  selected_fit$settings$mixture_components_initial <- as.integer(mixture_components)
  selected_fit$settings$mixture_components_final <- as.integer(selected_fit$settings$mixture_components)
  selected_fit$settings$mixture_components_max <- as.integer(mixture_components_max)
  selected_fit$settings$mixture_component_path <- as.integer(selection_history$mixture_components)
  selected_fit$settings$mixture_selection_history <- selection_history
  selected_fit$settings$mixture_selection_note <- selection_note

  bayes_finalize_trunc_comp_bayes_fit(selected_fit)
}
