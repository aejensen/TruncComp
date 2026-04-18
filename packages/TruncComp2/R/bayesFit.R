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

bayes_continuous_support <- function(continuous_support = c("real_line", "positive_real")) {
  if(length(continuous_support) == 0L || is.null(continuous_support)) {
    continuous_support <- "real_line"
  }

  match.arg(continuous_support, c("real_line", "positive_real"))
}

bayes_model_name <- function(continuous_support = c("real_line", "positive_real")) {
  continuous_support <- bayes_continuous_support(continuous_support)

  switch(
    continuous_support,
    real_line = "trunc_comp_bayes",
    positive_real = "trunc_comp_bayes_positive"
  )
}

bayes_default_prior <- function(continuous_support = c("real_line", "positive_real")) {
  continuous_support <- bayes_continuous_support(continuous_support)

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
    )
  )

  c(common, support_specific)
}

normalize_bayes_prior <- function(prior,
                                  continuous_support = c("real_line", "positive_real")) {
  continuous_support <- bayes_continuous_support(continuous_support)
  defaults <- bayes_default_prior(continuous_support)

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
  } else {
    positive_names <- c(positive_names, "mean_sdlog", "shape_sdlog")
    finite_names <- c("mean_meanlog", "shape_meanlog")
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
                                       continuous_support = c("real_line", "positive_real")) {
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
                                        continuous_support = c("real_line", "positive_real")) {
  continuous_support <- bayes_continuous_support(continuous_support)

  if(!identical(continuous_support, "positive_real")) {
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

build_bayes_standata <- function(data, atom, mixture_components, prior,
                                 continuous_support = c("real_line", "positive_real")) {
  continuous_support <- bayes_continuous_support(continuous_support)
  observed <- droplevels(data[data$A == 1, c("Y", "R"), drop = FALSE])
  scaling <- bayes_outcome_center_scale(observed$Y, continuous_support = continuous_support)

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

bayes_diagnostics <- function(fit, draws, parameters = bayes_parameter_names("contrast")) {
  diagnostic_draws <- bayes_select_draws(draws, parameters)
  convergence <- posterior::summarise_draws(
    diagnostic_draws,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  convergence <- as.data.frame(convergence)
  rownames(convergence) <- convergence$variable

  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(vapply(
    sampler_params,
    function(x) sum(x[, "divergent__"]),
    numeric(1)
  ))

  max_rhat <- max(convergence$rhat, na.rm = TRUE)
  min_bulk_ess <- min(convergence$ess_bulk, na.rm = TRUE)
  min_tail_ess <- min(convergence$ess_tail, na.rm = TRUE)

  list(
    divergences = as.integer(divergences),
    max_rhat = as.numeric(max_rhat),
    min_bulk_ess = as.numeric(min_bulk_ess),
    min_tail_ess = as.numeric(min_tail_ess),
    diagnostic_ok = isTRUE(divergences == 0) &&
      is.finite(max_rhat) &&
      max_rhat <= 1.01 &&
      is.finite(min_bulk_ess) &&
      min_bulk_ess >= 400 &&
      is.finite(min_tail_ess) &&
      min_tail_ess >= 400,
    parameter_table = convergence
  )
}

new_trunc_comp_bayes_fit <- function(fit = NULL, draws = NULL,
                                     summary_table = NULL, arm_table = NULL,
                                     diagnostics = NULL, settings = NULL,
                                     conf.level, success, error = "",
                                     data = NULL, atom = NULL, call = NULL) {
  out <- list(
    fit = fit,
    draws = draws,
    summary_table = summary_table,
    arm_table = arm_table,
    diagnostics = diagnostics,
    settings = settings,
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

fit_trunc_comp_bayes <- function(data, atom, conf.level = 0.95,
                                 continuous_support = c("real_line", "positive_real"),
                                 mixture_components = 10,
                                 chains = 4, iter_warmup = 1000,
                                 iter_sampling = 1000, seed = NULL,
                                 refresh = 0, control = list(adapt_delta = 0.95, max_treedepth = 12),
                                 prior = NULL, call = NULL,
                                 extra_args = list()) {
  conf.level <- validateConfidenceLevel(conf.level)
  continuous_support <- bayes_continuous_support(continuous_support)
  mixture_components <- validate_bayes_positive_integer(
    mixture_components,
    "mixture_components",
    min_value = 2L
  )
  chains <- validate_bayes_positive_integer(chains, "chains", min_value = 1L)
  iter_warmup <- validate_bayes_positive_integer(iter_warmup, "iter_warmup", min_value = 1L)
  iter_sampling <- validate_bayes_positive_integer(iter_sampling, "iter_sampling", min_value = 1L)
  seed <- validate_bayes_seed(seed)
  refresh <- validate_bayes_refresh(refresh)
  control <- validate_bayes_control(control)
  prior <- normalize_bayes_prior(prior, continuous_support = continuous_support)
  extra_args <- normalize_bayes_sampling_args(extra_args)
  model_name <- bayes_model_name(continuous_support)

  settings <- list(
    continuous_support = continuous_support,
    model_name = model_name,
    mixture_components = mixture_components,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = refresh,
    control = control,
    prior = prior,
    experimental = TRUE
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

  validate_bayes_support_data(data, continuous_support = continuous_support)

  if(!exists("stanmodels", inherits = TRUE)) {
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

  package_stanmodels <- get("stanmodels", inherits = TRUE)
  if(is.null(package_stanmodels[[model_name]])) {
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

  standata <- build_bayes_standata(
    data,
    atom,
    mixture_components,
    prior,
    continuous_support = continuous_support
  )
  settings$y_center <- standata$scaling$center
  settings$y_scale <- standata$scaling$scale

  sampling_args <- c(
    list(
      object = package_stanmodels[[standata$model_name]],
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

  summary_table <- bayes_summary_table(draws, conf.level)
  arm_table <- bayes_arm_table(draws, conf.level)
  diagnostics <- bayes_diagnostics(fit, draws)

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
