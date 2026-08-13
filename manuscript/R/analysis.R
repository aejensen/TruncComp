application_summary_stats <- function(d, metadata) {
  groups <- 0:1
  out <- lapply(groups, function(r) {
    arm <- d[d$R == r, , drop = FALSE]
    non_atom <- arm[arm$A == 1L, , drop = FALSE]
    q <- stats::quantile(non_atom$Y, probs = c(0.25, 0.75), names = FALSE, type = 2)

    data.frame(
      R = r,
      group = metadata$group_labels[[r + 1L]],
      n = nrow(arm),
      atom_n = sum(arm$A == 0L),
      atom_prop = mean(arm$A == 0L),
      non_atom_n = nrow(non_atom),
      non_atom_prop = mean(arm$A == 1L),
      survivor_mean = mean(non_atom$Y),
      survivor_sd = stats::sd(non_atom$Y),
      survivor_median = stats::median(non_atom$Y),
      survivor_q1 = q[[1]],
      survivor_q3 = q[[2]],
      survivor_min = min(non_atom$Y),
      survivor_max = max(non_atom$Y),
      combined_mean = mean(arm$Y),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

application_contrasts <- function(summary_stats) {
  control <- summary_stats[summary_stats$R == 0L, , drop = FALSE]
  treatment <- summary_stats[summary_stats$R == 1L, , drop = FALSE]

  data.frame(
    contrast = c(
      "Death risk difference",
      "Survival/non-atom probability difference",
      "Alive-with-EQ-VAS mean difference",
      "Coded mean contrast Delta"
    ),
    estimate = c(
      treatment$atom_prop - control$atom_prop,
      treatment$non_atom_prop - control$non_atom_prop,
      treatment$survivor_mean - control$survivor_mean,
      treatment$combined_mean - control$combined_mean
    ),
    interpretation = c(
      "rt-PA minus control probability of death by 6 months",
      "rt-PA minus control probability of being alive with observed EQ-VAS",
      "rt-PA minus control mean EQ-VAS among participants alive with observed EQ-VAS",
      "rt-PA minus control mean of death-coded combined endpoint"
    ),
    stringsAsFactors = FALSE
  )
}

run_standard_application_tests <- function(d) {
  non_atom <- d[d$A == 1L, , drop = FALSE]
  atom_tab <- table(factor(d$R, levels = 0:1), factor(d$A == 0L, levels = c(FALSE, TRUE)))

  data.frame(
    analysis = c(
      "Coded endpoint",
      "Coded endpoint",
      "Alive with EQ-VAS",
      "Alive with EQ-VAS",
      "Death atom"
    ),
    method = c(
      "Welch t-test",
      "Wilcoxon rank-sum",
      "Welch t-test",
      "Wilcoxon rank-sum",
      "Fisher exact test"
    ),
    statistic = c(
      -unname(stats::t.test(Y ~ R, data = d)$statistic),
      unname(suppressWarnings(stats::wilcox.test(Y ~ R, data = d)$statistic)),
      -unname(stats::t.test(Y ~ R, data = non_atom)$statistic),
      unname(suppressWarnings(stats::wilcox.test(Y ~ R, data = non_atom)$statistic)),
      NA_real_
    ),
    p_value = c(
      stats::t.test(Y ~ R, data = d)$p.value,
      suppressWarnings(stats::wilcox.test(Y ~ R, data = d)$p.value),
      stats::t.test(Y ~ R, data = non_atom)$p.value,
      suppressWarnings(stats::wilcox.test(Y ~ R, data = non_atom)$p.value),
      stats::fisher.test(atom_tab)$p.value
    ),
    stringsAsFactors = FALSE
  )
}

.safe_confint <- function(expr) {
  tryCatch(
    {
      out <- NULL
      utils::capture.output(out <- suppressMessages(suppressWarnings(expr)))
      out
    },
    error = function(e) structure(list(error = conditionMessage(e)), class = "application_ci_error")
  )
}

application_delta_intervals <- function(model) {
  methods <- c("welch", "profile", "projected")
  rows <- lapply(methods, function(method) {
    ci <- .safe_confint(stats::confint(
      model,
      parameter = "delta",
      method = method
    ))

    if (inherits(ci, "application_ci_error")) {
      return(data.frame(
        method = method,
        lower = NA_real_,
        upper = NA_real_,
        status = ci$error,
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      method = method,
      lower = unname(ci[1, 1]),
      upper = unname(ci[1, 2]),
      status = "estimated",
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

application_bayes_cache_path <- function(manuscript_dir) {
  file.path(.application_data_dir(manuscript_dir), "ist3-bayes-cache.rds")
}

application_bayes_candidate_cache_path <- function(manuscript_dir) {
  file.path(.application_data_dir(manuscript_dir), "ist3-bayes-cache-candidate.rds")
}

application_bayes_model_path <- function(manuscript_dir) {
  file.path(
    manuscript_dir,
    "stan",
    "ist3_ordered_bounded_score_logit_normal_h2.stan"
  )
}

application_bayes_settings <- function() {
  score_min <- 0
  score_max <- 100
  score_step <- 1
  heaping_grids <- c(1, 5, 10)
  list(
    continuous_support = "bounded_score",
    bounded_kernel = "logit_normal",
    score_min = score_min,
    score_max = score_max,
    score_step = score_step,
    heaping_grids = heaping_grids,
    heaping = "shared",
    mixture_components = 2L,
    auto_select_mixture_components = FALSE,
    init_strategy = "aligned_logit_normal_h2",
    chains = 4L,
    iter_warmup = 1000L,
    iter_sampling = 2000L,
    seed = 20260425L,
    control = list(adapt_delta = 0.995, max_treedepth = 12L),
    prior = list(
      rho_alpha = 1,
      rho_beta = 1,
      alpha_shape = 4,
      alpha_rate = 4,
      mu_logit_mean = 0,
      mu_logit_sd = 2,
      sigma_logit_meanlog = log(1.2),
      sigma_logit_sdlog = 0.5,
      eta_prior = rep(2, length(heaping_grids))
    )
  )
}

application_bayes_sampled_parameter_names <- function() {
  c(
    sprintf("v[%d,1]", 1:2),
    as.vector(outer(
      1:2,
      1:2,
      function(r, h) sprintf("mu_logit_comp[%d,%d]", r, h)
    )),
    as.vector(outer(
      1:2,
      1:2,
      function(r, h) sprintf("log_sigma_comp[%d,%d]", r, h)
    )),
    sprintf("alpha[%d]", 1:2),
    sprintf("rho[%d]", 1:2),
    sprintf("eta[1,%d]", 1:3)
  )
}

application_bayes_sampled_parameter_diagnostics <- function(fit) {
  if (!requireNamespace("rstan", quietly = TRUE) ||
      !requireNamespace("posterior", quietly = TRUE)) {
    stop(
      "The rstan and posterior packages are required to compute application diagnostics.",
      call. = FALSE
    )
  }

  stan_fit <- if (methods::is(fit, "stanfit")) fit else fit$fit
  if (!methods::is(stan_fit, "stanfit")) {
    stop("The application diagnostics require a fitted Stan model.", call. = FALSE)
  }

  parameter_bases <- c(
    "v",
    "mu_logit_comp",
    "log_sigma_comp",
    "alpha",
    "rho",
    "eta"
  )
  draws_array <- rstan::extract(
    stan_fit,
    pars = parameter_bases,
    permuted = FALSE,
    inc_warmup = FALSE
  )
  draws <- posterior::as_draws_array(draws_array)
  convergence <- as.data.frame(posterior::summarise_draws(
    draws,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  ))

  expected <- application_bayes_sampled_parameter_names()
  observed <- as.character(convergence$variable)
  if (!setequal(observed, expected) || length(observed) != length(expected)) {
    stop(
      paste(
        "The fitted application model parameters do not match the locked",
        "IST-3 parameter contract."
      ),
      call. = FALSE
    )
  }

  metric_names <- c("rhat", "ess_bulk", "ess_tail")
  if (any(!is.finite(as.matrix(convergence[metric_names])))) {
    stop("The application convergence diagnostics are not finite.", call. = FALSE)
  }

  list(
    max_rhat = max(convergence$rhat),
    min_bulk_ess = min(convergence$ess_bulk),
    min_tail_ess = min(convergence$ess_tail),
    parameter_table = convergence[match(expected, observed), , drop = FALSE]
  )
}

.application_sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required to validate Bayesian application artifacts.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("Cannot fingerprint missing file: %s", path), call. = FALSE)
  }

  digest::digest(file = path, algo = "sha256")
}

.application_sha256_object <- function(object) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required to validate Bayesian application artifacts.", call. = FALSE)
  }

  digest::digest(
    object,
    algo = "sha256",
    serialize = TRUE,
    serializeVersion = 2
  )
}

.application_sha256_text <- function(text) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required to validate Bayesian application artifacts.", call. = FALSE)
  }

  digest::digest(text, algo = "sha256", serialize = FALSE)
}

application_bayes_canonical_model_code <- function(code) {
  if (!(length(code) == 1L && is.character(code) && !is.na(code))) {
    stop("Stan model code must be a single character string.", call. = FALSE)
  }

  lines <- strsplit(gsub("\\r\\n?", "\n", code), "\n", fixed = TRUE)[[1L]]
  lines <- sub("[[:space:]]+$", "", lines)
  lines <- lines[nzchar(trimws(lines))]
  paste(lines, collapse = "\n")
}

application_bayes_model_code_fingerprint <- function(code) {
  .application_sha256_text(application_bayes_canonical_model_code(code))
}

application_bayes_normalize_data <- function(d) {
  required <- c("Y", "A", "R")
  if (!is.data.frame(d) || !all(required %in% names(d))) {
    stop("IST-3 Bayesian data must be a data frame containing Y, A, and R.", call. = FALSE)
  }

  data.frame(
    Y = as.numeric(d$Y),
    A = as.integer(d$A),
    R = as.integer(d$R)
  )
}

application_bayes_contract <- function(manuscript_dir) {
  list(
    model_id = "ist3_ordered_bounded_score_logit_normal_h2",
    model_path = application_bayes_model_path(manuscript_dir),
    locked_cache_sha256 = "6700cc423b91c2543b307b035a6f41748249d7befe3905004b956c703b356def",
    settings_sha256 = "14afe4fd3ed055f3a86e88a93e42dd7c36fcf933fef3e5ca60aefc6061f24d70",
    data_sha256 = "bd7ee54a9c59abcc72a27bfef7e76697c7122cfe6c89e0f4077ef1434f94432f",
    model_code_sha256 = "13f779073bc06cd585d91a61bd702aadf74065757f247dfec7955da17277dc90"
  )
}

application_bayes_init_function <- function(d, atom, settings) {
  if (!identical(settings$continuous_support, "bounded_score") ||
      !identical(settings$bounded_kernel, "logit_normal")) {
    return(NULL)
  }
  if (!identical(settings$mixture_components, 2L)) {
    stop("The dedicated IST-3 initialization requires exactly H = 2.", call. = FALSE)
  }

  k <- length(settings$heaping_grids)
  eta_groups <- if (identical(settings$heaping, "arm_specific")) 2L else 1L
  atom_rate <- as.numeric(tapply(d$Y == atom, factor(d$R, levels = 0:1), mean))
  eta_init <- settings$prior$eta_prior
  if (length(eta_init) == 1L) {
    eta_init <- rep(eta_init, k)
  }
  eta_init <- eta_init / sum(eta_init)

  function(chain_id = 1) {
    jitter_unit <- function(x, amount) {
      pmin(pmax(x + stats::runif(length(x), -amount, amount), 1e-4), 1 - 1e-4)
    }
    jitter_positive <- function(x, amount) {
      pmax(x * exp(stats::runif(length(x), -amount, amount)), 0.01)
    }

    common <- list(
      rho = jitter_unit(atom_rate, 0.005),
      eta = matrix(rep(eta_init, eta_groups), nrow = eta_groups, byrow = TRUE),
      alpha = rep(settings$prior$alpha_shape / settings$prior$alpha_rate, 2)
    )

    mu_unit <- rbind(c(0.50, 0.63), c(0.52, 0.67))
    sigma <- rbind(c(1.25, 1.00), c(1.20, 0.95))
    v <- matrix(c(0.70, 0.70), nrow = 2)
    mu_unit <- t(apply(mu_unit, 1, function(x) sort(jitter_unit(x, 0.01))))
    list(
      v = jitter_unit(v, 0.015),
      mu_logit_comp = matrix(stats::qlogis(as.vector(mu_unit)), nrow = 2),
      log_sigma_comp = log(matrix(jitter_positive(as.vector(sigma), 0.06), nrow = 2)),
      alpha = common$alpha,
      rho = common$rho,
      eta = common$eta
    )
  }
}

application_bayes_cached_model_code <- function(cached) {
  tryCatch(
    {
      stan_fit <- cached$fit$fit
      stan_model <- methods::slot(stan_fit, "stanmodel")
      code <- methods::slot(stan_model, "model_code")
      if (length(code) == 1L && is.character(code) && !is.na(code)) code else NULL
    },
    error = function(e) NULL
  )
}

application_bayes_validate_cache <- function(cached, d, cache_path, manuscript_dir) {
  contract <- application_bayes_contract(manuscript_dir)
  settings <- application_bayes_settings()
  errors <- character(0)
  add_error <- function(message) {
    errors <<- c(errors, message)
  }

  cache_sha256 <- if (file.exists(cache_path)) {
    .application_sha256_file(cache_path)
  } else {
    add_error("the locked cache file is missing")
    NA_character_
  }
  if (!is.na(cache_sha256) && !identical(cache_sha256, contract$locked_cache_sha256)) {
    add_error("the locked cache SHA-256 fingerprint does not match")
  }

  settings_sha256 <- .application_sha256_object(settings)
  if (!identical(settings_sha256, contract$settings_sha256)) {
    add_error("the configured application settings fingerprint does not match the locked contract")
  }
  if (!is.list(cached) || !identical(cached$settings, settings)) {
    add_error("the cached application settings do not exactly match the locked settings")
  }

  current_data <- tryCatch(application_bayes_normalize_data(d), error = identity)
  current_data_sha256 <- if (inherits(current_data, "error")) {
    add_error(conditionMessage(current_data))
    NA_character_
  } else {
    .application_sha256_object(current_data)
  }
  if (!is.na(current_data_sha256) && !identical(current_data_sha256, contract$data_sha256)) {
    add_error("the current IST-3 analysis data fingerprint does not match the locked contract")
  }

  cached_data <- tryCatch(
    application_bayes_normalize_data(cached$fit$data),
    error = identity
  )
  cached_data_sha256 <- if (inherits(cached_data, "error")) {
    add_error("the locked cache does not contain valid fitted data")
    NA_character_
  } else {
    .application_sha256_object(cached_data)
  }
  if (!is.na(cached_data_sha256) && !identical(cached_data_sha256, contract$data_sha256)) {
    add_error("the fitted data stored in the locked cache do not match the data fingerprint")
  }
  if (!inherits(current_data, "error") &&
      !inherits(cached_data, "error") &&
      !identical(current_data, cached_data)) {
    add_error("the current IST-3 data and the fitted data stored in the cache differ")
  }

  model_source <- if (file.exists(contract$model_path)) {
    paste(readLines(contract$model_path, warn = FALSE), collapse = "\n")
  } else {
    add_error("the dedicated IST-3 Stan model source is missing")
    NULL
  }
  model_source_sha256 <- if (is.null(model_source)) {
    NA_character_
  } else {
    application_bayes_model_code_fingerprint(model_source)
  }
  if (!is.na(model_source_sha256) &&
      !identical(model_source_sha256, contract$model_code_sha256)) {
    add_error("the dedicated IST-3 Stan source fingerprint does not match the locked contract")
  }

  cached_model_code <- application_bayes_cached_model_code(cached)
  cached_model_code_sha256 <- if (is.null(cached_model_code)) {
    add_error("the locked cache does not contain embedded Stan model code")
    NA_character_
  } else {
    application_bayes_model_code_fingerprint(cached_model_code)
  }
  if (!is.na(cached_model_code_sha256) &&
      !identical(cached_model_code_sha256, contract$model_code_sha256)) {
    add_error("the Stan model code embedded in the locked cache does not match the model fingerprint")
  }
  if (!(is.list(cached) &&
        isTRUE(cached$success) &&
        is.list(cached$fit) &&
        isTRUE(cached$fit$success))) {
    add_error("the locked cache does not contain a successful fit")
  }

  list(
    valid = length(errors) == 0L,
    errors = errors,
    fingerprints = list(
      cache_sha256 = cache_sha256,
      settings_sha256 = settings_sha256,
      current_data_sha256 = current_data_sha256,
      cached_data_sha256 = cached_data_sha256,
      model_source_sha256 = model_source_sha256,
      cached_model_code_sha256 = cached_model_code_sha256
    )
  )
}

application_bayes_validate_candidate <- function(candidate, d, manuscript_dir) {
  contract <- application_bayes_contract(manuscript_dir)
  settings <- application_bayes_settings()
  errors <- character(0)
  add_error <- function(message) {
    errors <<- c(errors, message)
  }

  if (!(is.list(candidate) &&
        isTRUE(candidate$success) &&
        is.list(candidate$fit) &&
        isTRUE(candidate$fit$success))) {
    add_error("the candidate does not contain a successful fit")
  }
  if (!identical(candidate$model_id, contract$model_id)) {
    add_error("the candidate model identifier does not match the application model")
  }
  if (!identical(candidate$fit$settings$model_name, contract$model_id)) {
    add_error("the fitted model identifier does not match the application model")
  }

  settings_sha256 <- .application_sha256_object(settings)
  if (!identical(settings_sha256, contract$settings_sha256)) {
    add_error("the configured application settings fingerprint does not match the contract")
  }
  if (!identical(candidate$settings, settings)) {
    add_error("the candidate settings do not exactly match the application settings")
  }

  current_data <- tryCatch(application_bayes_normalize_data(d), error = identity)
  current_data_sha256 <- if (inherits(current_data, "error")) {
    add_error(conditionMessage(current_data))
    NA_character_
  } else {
    .application_sha256_object(current_data)
  }
  if (!is.na(current_data_sha256) && !identical(current_data_sha256, contract$data_sha256)) {
    add_error("the current IST-3 analysis data fingerprint does not match the contract")
  }

  fitted_data <- tryCatch(
    application_bayes_normalize_data(candidate$fit$data),
    error = identity
  )
  fitted_data_sha256 <- if (inherits(fitted_data, "error")) {
    add_error("the candidate does not contain valid fitted data")
    NA_character_
  } else {
    .application_sha256_object(fitted_data)
  }
  if (!is.na(fitted_data_sha256) && !identical(fitted_data_sha256, contract$data_sha256)) {
    add_error("the candidate fitted-data fingerprint does not match the contract")
  }
  if (!inherits(current_data, "error") &&
      !inherits(fitted_data, "error") &&
      !identical(current_data, fitted_data)) {
    add_error("the current IST-3 data and the candidate fitted data differ")
  }

  model_source <- if (file.exists(contract$model_path)) {
    paste(readLines(contract$model_path, warn = FALSE), collapse = "\n")
  } else {
    add_error("the dedicated IST-3 Stan model source is missing")
    NULL
  }
  model_source_sha256 <- if (is.null(model_source)) {
    NA_character_
  } else {
    application_bayes_model_code_fingerprint(model_source)
  }
  if (!is.na(model_source_sha256) &&
      !identical(model_source_sha256, contract$model_code_sha256)) {
    add_error("the dedicated IST-3 Stan source fingerprint does not match the contract")
  }

  fitted_model_code <- application_bayes_cached_model_code(candidate)
  fitted_model_code_sha256 <- if (is.null(fitted_model_code)) {
    add_error("the candidate does not contain embedded Stan model code")
    NA_character_
  } else {
    application_bayes_model_code_fingerprint(fitted_model_code)
  }
  if (!is.na(fitted_model_code_sha256) &&
      !identical(fitted_model_code_sha256, contract$model_code_sha256)) {
    add_error("the candidate embedded Stan model code does not match the model fingerprint")
  }

  contract_fields <- c(
    "model_id",
    "settings_sha256",
    "data_sha256",
    "model_code_sha256"
  )
  if (!is.list(candidate$contract) ||
      !identical(candidate$contract[contract_fields], contract[contract_fields])) {
    add_error("the candidate contract metadata does not match the application contract")
  }

  list(
    valid = length(errors) == 0L,
    errors = errors,
    fingerprints = list(
      settings_sha256 = settings_sha256,
      current_data_sha256 = current_data_sha256,
      fitted_data_sha256 = fitted_data_sha256,
      model_source_sha256 = model_source_sha256,
      fitted_model_code_sha256 = fitted_model_code_sha256
    )
  )
}

application_bayes_fit_ordered_h2 <- function(d, atom, settings, manuscript_dir) {
  contract <- application_bayes_contract(manuscript_dir)
  if (!identical(settings$mixture_components, 2L)) {
    stop("The dedicated IST-3 application model requires exactly H = 2.", call. = FALSE)
  }
  if (!file.exists(contract$model_path)) {
    stop("The dedicated IST-3 Stan model source is missing.", call. = FALSE)
  }

  source_code <- paste(readLines(contract$model_path, warn = FALSE), collapse = "\n")
  source_sha256 <- application_bayes_model_code_fingerprint(source_code)
  if (!identical(source_sha256, contract$model_code_sha256)) {
    stop("The dedicated IST-3 Stan model source does not match the locked model fingerprint.", call. = FALSE)
  }
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("The rstan package is required for an explicitly requested IST-3 refresh.", call. = FALSE)
  }

  fit_function <- getFromNamespace("fit_trunc_comp_bayes", "TruncComp")
  required_formals <- c("model_object", "model_name_override")
  if (!all(required_formals %in% names(formals(fit_function)))) {
    stop(
      "The local TruncComp source is required for the dedicated IST-3 application model.",
      call. = FALSE
    )
  }

  previous_mc_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = previous_mc_cores), add = TRUE)
  options(mc.cores = min(settings$chains, max(1L, parallel::detectCores(logical = FALSE))))
  model_object <- rstan::stan_model(
    file = contract$model_path,
    model_name = contract$model_id
  )

  fit_function(
    data = application_bayes_normalize_data(d),
    atom = atom,
    conf.level = 0.95,
    continuous_support = settings$continuous_support,
    bounded_kernel = settings$bounded_kernel,
    bounded_kernel_supplied = TRUE,
    score_min = settings$score_min,
    score_max = settings$score_max,
    score_step = settings$score_step,
    heaping_grids = settings$heaping_grids,
    heaping = settings$heaping,
    support_supplied = list(
      score_min = TRUE,
      score_max = TRUE,
      score_step = TRUE,
      heaping_grids = TRUE,
      heaping = TRUE
    ),
    mixture_components = settings$mixture_components,
    auto_select_mixture_components = FALSE,
    mixture_components_max = settings$mixture_components,
    chains = settings$chains,
    iter_warmup = settings$iter_warmup,
    iter_sampling = settings$iter_sampling,
    seed = settings$seed,
    refresh = 0,
    control = settings$control,
    prior = settings$prior,
    call = match.call(),
    extra_args = list(init = application_bayes_init_function(d, atom, settings)),
    model_object = model_object,
    model_name_override = contract$model_id
  )
}

application_bayes_result_from_fit <- function(fit, settings, manuscript_dir) {
  if (inherits(fit, "error")) {
    stop(conditionMessage(fit), call. = FALSE)
  }
  if (!isTRUE(fit$success)) {
    stop(sprintf("The explicitly requested IST-3 fit failed: %s", fit$error), call. = FALSE)
  }

  ppc <- tryCatch(
    TruncComp::posterior_predictive_pvalues(
      fit,
      ndraws = 200,
      seed = settings$seed + 1L
    ),
    error = function(e) e
  )
  draws <- fit$draws
  probabilities <- c(
    delta_gt_0 = mean(draws$delta > 0),
    mu_delta_gt_0 = mean(draws$mu_delta > 0),
    alpha_delta_gt_1 = mean(draws$alpha_delta > 1),
    death_risk_reduction = mean(draws$delta_atom < 0),
    survival_probability_increase = mean((draws$pi_1 - draws$pi_0) > 0)
  )
  contract <- application_bayes_contract(manuscript_dir)

  list(
    success = TRUE,
    fit = fit,
    summary_table = fit$summary_table,
    arm_table = fit$arm_table,
    diagnostics = fit$diagnostics,
    ppc_table = if (inherits(ppc, "error")) NULL else ppc,
    ppc_error = if (inherits(ppc, "error")) conditionMessage(ppc) else NULL,
    probabilities = probabilities,
    settings = settings,
    cache_path = application_bayes_candidate_cache_path(manuscript_dir),
    from_cache = FALSE,
    model_id = contract$model_id,
    contract = contract
  )
}

.application_env_flag <- function(name) {
  identical(tolower(trimws(Sys.getenv(name, "false"))), "true")
}

fit_or_load_application_bayes <- function(d, metadata, manuscript_dir) {
  cache_path <- application_bayes_cache_path(manuscript_dir)
  candidate_path <- application_bayes_candidate_cache_path(manuscript_dir)
  refresh <- .application_env_flag("TRUNCCOMP_REFRESH_IST3_BAYES")
  skip <- .application_env_flag("TRUNCCOMP_SKIP_IST3_BAYES")
  settings <- application_bayes_settings()
  contract <- application_bayes_contract(manuscript_dir)

  if (skip) {
    stop(
      "The manuscript requires its locked IST-3 Bayesian cache and cannot be built with TRUNCCOMP_SKIP_IST3_BAYES=true.",
      call. = FALSE
    )
  }

  if (!refresh) {
    if (!file.exists(cache_path)) {
      stop(
        paste(
          "The locked IST-3 Bayesian cache is missing.",
          "No Bayesian fit was started because refresh was not explicitly requested."
        ),
        call. = FALSE
      )
    }

    cached <- tryCatch(readRDS(cache_path), error = identity)
    if (inherits(cached, "error")) {
      stop(
        paste(
          "The locked IST-3 Bayesian cache could not be read.",
          "No Bayesian fit was started because refresh was not explicitly requested."
        ),
        call. = FALSE
      )
    }
    validation <- application_bayes_validate_cache(
      cached = cached,
      d = d,
      cache_path = cache_path,
      manuscript_dir = manuscript_dir
    )
    if (!validation$valid) {
      stop(
        paste0(
          "The locked IST-3 Bayesian cache failed validation: ",
          paste(validation$errors, collapse = "; "),
          ". No Bayesian fit was started because refresh was not explicitly requested."
        ),
        call. = FALSE
      )
    }

    cached$cache_path <- cache_path
    cached$from_cache <- TRUE
    cached$model_id <- contract$model_id
    cached$contract <- c(contract, validation$fingerprints)
    return(cached)
  }

  fit <- tryCatch(
    application_bayes_fit_ordered_h2(
      d = d,
      atom = metadata$atom,
      settings = settings,
      manuscript_dir = manuscript_dir
    ),
    error = identity
  )
  out <- application_bayes_result_from_fit(fit, settings, manuscript_dir)
  in_memory_validation <- application_bayes_validate_candidate(
    candidate = out,
    d = d,
    manuscript_dir = manuscript_dir
  )
  if (!in_memory_validation$valid) {
    stop(
      paste0(
        "The explicitly refreshed IST-3 candidate failed in-memory validation: ",
        paste(in_memory_validation$errors, collapse = "; "),
        "."
      ),
      call. = FALSE
    )
  }

  saveRDS(out, candidate_path)
  saved_candidate <- tryCatch(readRDS(candidate_path), error = identity)
  if (inherits(saved_candidate, "error")) {
    unlink(candidate_path)
    stop("The saved IST-3 candidate cache could not be read and was removed.", call. = FALSE)
  }
  saved_validation <- application_bayes_validate_candidate(
    candidate = saved_candidate,
    d = d,
    manuscript_dir = manuscript_dir
  )
  if (!saved_validation$valid) {
    unlink(candidate_path)
    stop(
      paste0(
        "The saved IST-3 candidate cache failed validation and was removed: ",
        paste(saved_validation$errors, collapse = "; "),
        "."
      ),
      call. = FALSE
    )
  }
  message(
    "Saved the explicitly refreshed IST-3 fit as a candidate cache at ",
    candidate_path,
    ". The locked publication cache was not modified."
  )
  saved_candidate
}

analyze_application_data <- function(application_data, manuscript_dir) {
  d <- application_data$data
  metadata <- application_data$metadata
  atom <- metadata$atom

  model_lrt <- TruncComp::trunc_comp(Y ~ R, atom = atom, data = d, method = "lrt")
  model_splrt <- TruncComp::trunc_comp(Y ~ R, atom = atom, data = d, method = "splrt")
  surface <- suppressMessages(stats::confint(
    model_splrt,
    parameter = "joint",
    resolution = metadata$surface_resolution,
    plot = FALSE
  ))
  delta_intervals <- application_delta_intervals(model_splrt)
  summary_stats <- application_summary_stats(d, metadata)
  contrasts <- application_contrasts(summary_stats)
  standard_tests <- run_standard_application_tests(d)
  bayes <- fit_or_load_application_bayes(d, metadata, manuscript_dir)

  list(
    data = d,
    raw = application_data$raw,
    metadata = metadata,
    summary_stats = summary_stats,
    contrasts = contrasts,
    standard_tests = standard_tests,
    model_lrt = model_lrt,
    model_splrt = model_splrt,
    surface = surface,
    delta_intervals = delta_intervals,
    bayes = bayes
  )
}

compute_application_results <- function(manuscript_dir) {
  analyze_application_data(load_application_data(manuscript_dir), manuscript_dir)
}
