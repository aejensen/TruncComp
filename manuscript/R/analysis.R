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
      "Combined mean contrast Delta"
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
      "Combined endpoint",
      "Combined endpoint",
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
      unname(stats::t.test(Y ~ R, data = d)$statistic),
      unname(suppressWarnings(stats::wilcox.test(Y ~ R, data = d)$statistic)),
      unname(stats::t.test(Y ~ R, data = non_atom)$statistic),
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

application_delta_intervals <- function(model, resolution = 45) {
  methods <- c("welch", "profile", "projected")
  rows <- lapply(methods, function(method) {
    ci <- .safe_confint(stats::confint(
      model,
      parameter = "delta",
      method = method,
      algorithm = "optimize",
      resolution = resolution
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

application_bayes_settings <- function() {
  as_int <- function(name, default) {
    value <- Sys.getenv(name, as.character(default))
    out <- suppressWarnings(as.integer(value))
    if (is.na(out)) default else out
  }

  score_min <- 0
  score_max <- 100
  score_step <- 1
  heaping_grids <- c(1, 5, 10)
  mixture_components <- as_int("TRUNCCOMP_IST3_BAYES_COMPONENTS", 2L)
  list(
    continuous_support = "bounded_score",
    bounded_kernel = "logit_normal",
    score_min = score_min,
    score_max = score_max,
    score_step = score_step,
    heaping_grids = heaping_grids,
    heaping = "shared",
    mixture_components = mixture_components,
    auto_select_mixture_components = FALSE,
    init_strategy = "aligned_logit_normal_h2",
    chains = as_int("TRUNCCOMP_IST3_BAYES_CHAINS", 4L),
    iter_warmup = as_int("TRUNCCOMP_IST3_BAYES_WARMUP", 1000L),
    iter_sampling = as_int("TRUNCCOMP_IST3_BAYES_SAMPLING", 2000L),
    seed = as_int("TRUNCCOMP_IST3_BAYES_SEED", 20260425L),
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

application_bayes_init_function <- function(d, atom, settings) {
  if (!identical(settings$continuous_support, "bounded_score") ||
      !identical(settings$bounded_kernel, "logit_normal")) {
    return(NULL)
  }

  h <- settings$mixture_components
  k <- length(settings$heaping_grids)
  eta_groups <- if (identical(settings$heaping, "arm_specific")) 2L else 1L
  score_range <- settings$score_max - settings$score_min
  x_unit <- (d$Y[d$Y != atom] - settings$score_min) / score_range
  x_unit <- pmin(pmax(x_unit, 0.05), 0.95)
  centers <- as.numeric(stats::quantile(x_unit, probs = seq(0.35, 0.75, length.out = h), names = FALSE))
  centers <- pmin(pmax(centers, 0.05), 0.95)
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

    if (h == 1L) {
      arm_centers <- vapply(0:1, function(r) {
        values <- (d$Y[d$R == r & d$Y != atom] - settings$score_min) / score_range
        stats::median(values, na.rm = TRUE)
      }, numeric(1))
      return(c(
        list(
          mu_logit_comp = matrix(stats::qlogis(jitter_unit(arm_centers, 0.02)), nrow = 2),
          log_sigma_comp = matrix(log(jitter_positive(rep(1.2, 2), 0.05)), nrow = 2)
        ),
        common
      ))
    }

    if (h == 2L) {
      mu_unit <- rbind(c(0.50, 0.63), c(0.52, 0.67))
      sigma <- rbind(c(1.25, 1.00), c(1.20, 0.95))
      v <- matrix(c(0.70, 0.70), nrow = 2)
    } else {
      mu_unit <- t(replicate(2, sort(jitter_unit(centers, 0.04))))
      sigma <- matrix(jitter_positive(rep(1.2, 2 * h), 0.08), nrow = 2)
      v <- matrix(jitter_unit(rep(0.55, 2 * (h - 1L)), 0.05), nrow = 2)
    }

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

application_bayes_cache_matches <- function(cached, settings) {
  if (!is.list(cached) || is.null(cached$settings)) {
    return(FALSE)
  }

  identical(cached$settings, settings)
}

fit_or_load_application_bayes <- function(d, metadata, manuscript_dir) {
  cache_path <- application_bayes_cache_path(manuscript_dir)
  refresh <- identical(tolower(Sys.getenv("TRUNCCOMP_REFRESH_IST3_BAYES", "false")), "true")
  skip <- identical(tolower(Sys.getenv("TRUNCCOMP_SKIP_IST3_BAYES", "false")), "true")
  settings <- application_bayes_settings()

  if (file.exists(cache_path) && !refresh) {
    cached <- readRDS(cache_path)
    if (application_bayes_cache_matches(cached, settings)) {
      cached$cache_path <- cache_path
      cached$from_cache <- TRUE
      return(cached)
    }
    message("Ignoring stale IST-3 Bayesian cache because the model settings have changed.")
  }
  if (skip) {
    return(list(
      success = FALSE,
      error = "Bayesian fitting skipped because TRUNCCOMP_SKIP_IST3_BAYES=true.",
      cache_path = cache_path,
      from_cache = FALSE
    ))
  }

  if (!requireNamespace("rstan", quietly = TRUE)) {
    return(list(
      success = FALSE,
      error = "The rstan package is not available.",
      cache_path = cache_path,
      from_cache = FALSE
    ))
  }

  options(mc.cores = min(settings$chains, max(1L, parallel::detectCores(logical = FALSE))))
  rstan::rstan_options(auto_write = TRUE)

  fit <- tryCatch(
    TruncComp2::trunc_comp_bayes(
      Y ~ R,
      atom = metadata$atom,
      data = d,
      continuous_support = settings$continuous_support,
      bounded_kernel = settings$bounded_kernel,
      score_min = settings$score_min,
      score_max = settings$score_max,
      score_step = settings$score_step,
      heaping_grids = settings$heaping_grids,
      heaping = settings$heaping,
      mixture_components = settings$mixture_components,
      auto_select_mixture_components = settings$auto_select_mixture_components,
      chains = settings$chains,
      iter_warmup = settings$iter_warmup,
      iter_sampling = settings$iter_sampling,
      seed = settings$seed,
      refresh = 0,
      control = settings$control,
      prior = settings$prior,
      init = application_bayes_init_function(d, metadata$atom, settings)
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    out <- list(
      success = FALSE,
      error = conditionMessage(fit),
      settings = settings,
      cache_path = cache_path,
      from_cache = FALSE
    )
    saveRDS(out, cache_path)
    return(out)
  }

  ppc <- tryCatch(
    TruncComp2::posterior_predictive_pvalues(fit, ndraws = 200, seed = settings$seed + 1L),
    error = function(e) e
  )
  draws <- fit$draws
  probabilities <- if (isTRUE(fit$success)) {
    c(
      delta_gt_0 = mean(draws$delta > 0),
      mu_delta_gt_0 = mean(draws$mu_delta > 0),
      alpha_delta_gt_1 = mean(draws$alpha_delta > 1),
      death_risk_reduction = mean(draws$delta_atom < 0),
      survival_probability_increase = mean((draws$pi_1 - draws$pi_0) > 0)
    )
  } else {
    rep(NA_real_, 5)
  }

  out <- list(
    success = isTRUE(fit$success),
    fit = fit,
    summary_table = fit$summary_table,
    arm_table = fit$arm_table,
    diagnostics = fit$diagnostics,
    ppc_table = if (inherits(ppc, "error")) NULL else ppc,
    ppc_error = if (inherits(ppc, "error")) conditionMessage(ppc) else NULL,
    probabilities = probabilities,
    settings = settings,
    cache_path = cache_path,
    from_cache = FALSE
  )
  saveRDS(out, cache_path)
  out
}

analyze_application_data <- function(application_data, manuscript_dir) {
  d <- application_data$data
  metadata <- application_data$metadata
  atom <- metadata$atom

  model_lrt <- TruncComp2::trunc_comp(Y ~ R, atom = atom, data = d, method = "lrt")
  model_splrt <- TruncComp2::trunc_comp(Y ~ R, atom = atom, data = d, method = "splrt")
  surface <- suppressMessages(stats::confint(
    model_splrt,
    parameter = "joint",
    resolution = metadata$surface_resolution,
    plot = FALSE
  ))
  delta_intervals <- application_delta_intervals(model_splrt, resolution = metadata$surface_resolution)
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

liver_appendix_summary_stats <- function(d, metadata) {
  groups <- 0:1
  out <- lapply(groups, function(r) {
    arm <- d[d$R == r, , drop = FALSE]
    non_atom <- arm[arm$A == 1L, , drop = FALSE]
    q <- stats::quantile(non_atom$Y, probs = c(0.25, 0.75), names = FALSE, type = 2)
    lag_q <- stats::quantile(non_atom$measurement_lag, probs = c(0.25, 0.75), names = FALSE, type = 2)

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
      lag_mean = mean(non_atom$measurement_lag),
      lag_sd = stats::sd(non_atom$measurement_lag),
      lag_median = stats::median(non_atom$measurement_lag),
      lag_q1 = lag_q[[1]],
      lag_q3 = lag_q[[2]],
      lag_min = min(non_atom$measurement_lag),
      lag_max = max(non_atom$measurement_lag),
      combined_mean = mean(arm$Y),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

liver_appendix_contrasts <- function(summary_stats) {
  control <- summary_stats[summary_stats$R == 0L, , drop = FALSE]
  treatment <- summary_stats[summary_stats$R == 1L, , drop = FALSE]

  data.frame(
    contrast = c(
      "Death risk difference",
      "Non-atom probability difference",
      "Prothrombin mean difference",
      "Combined mean contrast Delta"
    ),
    estimate = c(
      treatment$atom_prop - control$atom_prop,
      treatment$non_atom_prop - control$non_atom_prop,
      treatment$survivor_mean - control$survivor_mean,
      treatment$combined_mean - control$combined_mean
    ),
    interpretation = c(
      "prednisone minus placebo probability of death by 2 years",
      "prednisone minus placebo probability of being known alive at 2 years with prothrombin observed",
      "prednisone minus placebo mean prothrombin index among known 2-year survivors",
      "prednisone minus placebo mean of the death-coded two-year endpoint"
    ),
    stringsAsFactors = FALSE
  )
}

run_standard_liver_appendix_tests <- function(d) {
  non_atom <- d[d$A == 1L, , drop = FALSE]
  atom_tab <- table(factor(d$R, levels = 0:1), factor(d$A == 0L, levels = c(FALSE, TRUE)))

  data.frame(
    analysis = c(
      "Combined endpoint",
      "Combined endpoint",
      "Non-atom component",
      "Non-atom component",
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
      unname(stats::t.test(Y ~ R, data = d)$statistic),
      unname(suppressWarnings(stats::wilcox.test(Y ~ R, data = d)$statistic)),
      unname(stats::t.test(Y ~ R, data = non_atom)$statistic),
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

liver_appendix_bayes_cache_path <- function(manuscript_dir) {
  file.path(.application_data_dir(manuscript_dir), "liver-appendix-bayes-cache.rds")
}

liver_appendix_bayes_settings <- function() {
  as_int <- function(name, default) {
    value <- Sys.getenv(name, as.character(default))
    out <- suppressWarnings(as.integer(value))
    if (is.na(out)) default else out
  }

  list(
    continuous_support = "positive_real",
    mixture_components = as_int("TRUNCCOMP_LIVER_BAYES_COMPONENTS", 4L),
    auto_select_mixture_components = FALSE,
    chains = as_int("TRUNCCOMP_LIVER_BAYES_CHAINS", 4L),
    iter_warmup = as_int("TRUNCCOMP_LIVER_BAYES_WARMUP", 300L),
    iter_sampling = as_int("TRUNCCOMP_LIVER_BAYES_SAMPLING", 300L),
    seed = as_int("TRUNCCOMP_LIVER_BAYES_SEED", 20260420L),
    control = list(adapt_delta = 0.99, max_treedepth = 12L),
    prior = list(
      rho_alpha = 1,
      rho_beta = 1,
      alpha_shape = 2,
      alpha_rate = 1,
      mean_meanlog = 0,
      mean_sdlog = 0.5,
      shape_meanlog = log(2),
      shape_sdlog = 0.5
    )
  )
}

fit_or_load_liver_appendix_bayes <- function(d, metadata, manuscript_dir) {
  cache_path <- liver_appendix_bayes_cache_path(manuscript_dir)
  refresh <- identical(tolower(Sys.getenv("TRUNCCOMP_REFRESH_LIVER_BAYES", "false")), "true")
  skip <- identical(tolower(Sys.getenv("TRUNCCOMP_SKIP_LIVER_BAYES", "false")), "true")

  if (file.exists(cache_path) && !refresh) {
    cached <- readRDS(cache_path)
    cached$cache_path <- cache_path
    cached$from_cache <- TRUE
    return(cached)
  }
  if (skip) {
    return(list(
      success = FALSE,
      error = "Bayesian fitting skipped because TRUNCCOMP_SKIP_LIVER_BAYES=true.",
      cache_path = cache_path,
      from_cache = FALSE
    ))
  }

  if (!requireNamespace("rstan", quietly = TRUE)) {
    return(list(
      success = FALSE,
      error = "The rstan package is not available.",
      cache_path = cache_path,
      from_cache = FALSE
    ))
  }

  settings <- liver_appendix_bayes_settings()
  options(mc.cores = min(settings$chains, max(1L, parallel::detectCores(logical = FALSE))))
  rstan::rstan_options(auto_write = TRUE)

  fit <- tryCatch(
    TruncComp2::trunc_comp_bayes(
      Y ~ R,
      atom = metadata$atom,
      data = d,
      continuous_support = settings$continuous_support,
      mixture_components = settings$mixture_components,
      auto_select_mixture_components = settings$auto_select_mixture_components,
      chains = settings$chains,
      iter_warmup = settings$iter_warmup,
      iter_sampling = settings$iter_sampling,
      seed = settings$seed,
      refresh = 0,
      control = settings$control,
      prior = settings$prior
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    out <- list(
      success = FALSE,
      error = conditionMessage(fit),
      settings = settings,
      cache_path = cache_path,
      from_cache = FALSE
    )
    saveRDS(out, cache_path)
    return(out)
  }

  ppc <- tryCatch(
    TruncComp2::posterior_predictive_pvalues(fit, ndraws = 200, seed = settings$seed + 1L),
    error = function(e) e
  )
  draws <- fit$draws
  probabilities <- if (isTRUE(fit$success)) {
    c(
      delta_gt_0 = mean(draws$delta > 0),
      mu_delta_gt_0 = mean(draws$mu_delta > 0),
      alpha_delta_gt_1 = mean(draws$alpha_delta > 1),
      death_risk_reduction = mean(draws$delta_atom < 0),
      non_atom_probability_increase = mean((draws$pi_1 - draws$pi_0) > 0)
    )
  } else {
    rep(NA_real_, 5)
  }

  out <- list(
    success = isTRUE(fit$success),
    fit = fit,
    summary_table = fit$summary_table,
    arm_table = fit$arm_table,
    diagnostics = fit$diagnostics,
    ppc_table = if (inherits(ppc, "error")) NULL else ppc,
    ppc_error = if (inherits(ppc, "error")) conditionMessage(ppc) else NULL,
    probabilities = probabilities,
    settings = settings,
    cache_path = cache_path,
    from_cache = FALSE
  )
  saveRDS(out, cache_path)
  out
}

analyze_liver_appendix_data <- function(liver_data, manuscript_dir) {
  d <- liver_data$data
  metadata <- liver_data$metadata
  atom <- metadata$atom

  model_lrt <- TruncComp2::trunc_comp(Y ~ R, atom = atom, data = d, method = "lrt")
  model_splrt <- TruncComp2::trunc_comp(Y ~ R, atom = atom, data = d, method = "splrt")
  surface <- suppressMessages(stats::confint(
    model_splrt,
    parameter = "joint",
    resolution = metadata$surface_resolution,
    plot = FALSE
  ))
  delta_intervals <- application_delta_intervals(model_splrt, resolution = metadata$surface_resolution)
  summary_stats <- liver_appendix_summary_stats(d, metadata)
  contrasts <- liver_appendix_contrasts(summary_stats)
  standard_tests <- run_standard_liver_appendix_tests(d)
  bayes <- fit_or_load_liver_appendix_bayes(d, metadata, manuscript_dir)

  list(
    data = d,
    raw = liver_data$raw,
    subjects = liver_data$subjects,
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

compute_liver_appendix_results <- function(manuscript_dir) {
  analyze_liver_appendix_data(load_liver_appendix_data(manuscript_dir), manuscript_dir)
}

licorice_appendix_summary_stats <- function(d, metadata) {
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
      survivor_distinct = length(unique(non_atom$Y)),
      combined_mean = mean(arm$Y),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

licorice_appendix_contrasts <- function(summary_stats) {
  control <- summary_stats[summary_stats$R == 0L, , drop = FALSE]
  treatment <- summary_stats[summary_stats$R == 1L, , drop = FALSE]

  data.frame(
    contrast = c(
      "No-pain atom probability difference",
      "Any-pain probability difference",
      "Positive pain mean difference",
      "Combined mean contrast Delta"
    ),
    estimate = c(
      treatment$atom_prop - control$atom_prop,
      treatment$non_atom_prop - control$non_atom_prop,
      treatment$survivor_mean - control$survivor_mean,
      treatment$combined_mean - control$combined_mean
    ),
    interpretation = c(
      "licorice minus sugar-water probability of no swallowing pain at 30 minutes",
      "licorice minus sugar-water probability of any swallowing pain at 30 minutes",
      "licorice minus sugar-water mean swallowing-pain score among participants with positive pain",
      "licorice minus sugar-water mean of the no-pain-at-zero combined endpoint"
    ),
    stringsAsFactors = FALSE
  )
}

run_standard_licorice_appendix_tests <- function(d) {
  non_atom <- d[d$A == 1L, , drop = FALSE]
  atom_tab <- table(factor(d$R, levels = 0:1), factor(d$A == 0L, levels = c(FALSE, TRUE)))

  data.frame(
    analysis = c(
      "Combined endpoint",
      "Combined endpoint",
      "Positive pain component",
      "Positive pain component",
      "No-pain atom"
    ),
    method = c(
      "Welch t-test",
      "Wilcoxon rank-sum",
      "Welch t-test",
      "Wilcoxon rank-sum",
      "Fisher exact test"
    ),
    statistic = c(
      unname(stats::t.test(Y ~ R, data = d)$statistic),
      unname(suppressWarnings(stats::wilcox.test(Y ~ R, data = d)$statistic)),
      unname(stats::t.test(Y ~ R, data = non_atom)$statistic),
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

licorice_appendix_bayes_cache_path <- function(manuscript_dir) {
  file.path(.application_data_dir(manuscript_dir), "licorice-appendix-bayes-cache.rds")
}

licorice_appendix_bayes_settings <- function() {
  as_int <- function(name, default) {
    value <- Sys.getenv(name, as.character(default))
    out <- suppressWarnings(as.integer(value))
    if (is.na(out)) default else out
  }

  list(
    continuous_support = "positive_real",
    mixture_components = as_int("TRUNCCOMP_LICORICE_BAYES_COMPONENTS", 2L),
    auto_select_mixture_components = FALSE,
    chains = as_int("TRUNCCOMP_LICORICE_BAYES_CHAINS", 4L),
    iter_warmup = as_int("TRUNCCOMP_LICORICE_BAYES_WARMUP", 300L),
    iter_sampling = as_int("TRUNCCOMP_LICORICE_BAYES_SAMPLING", 300L),
    seed = as_int("TRUNCCOMP_LICORICE_BAYES_SEED", 20260421L),
    control = list(adapt_delta = 0.99, max_treedepth = 12L),
    prior = list(
      rho_alpha = 1,
      rho_beta = 1,
      alpha_shape = 2,
      alpha_rate = 1,
      mean_meanlog = log(2),
      mean_sdlog = 0.75,
      shape_meanlog = log(2),
      shape_sdlog = 0.75
    )
  )
}

fit_or_load_licorice_appendix_bayes <- function(d, metadata, manuscript_dir) {
  cache_path <- licorice_appendix_bayes_cache_path(manuscript_dir)
  refresh <- identical(tolower(Sys.getenv("TRUNCCOMP_REFRESH_LICORICE_BAYES", "false")), "true")
  skip <- identical(tolower(Sys.getenv("TRUNCCOMP_SKIP_LICORICE_BAYES", "false")), "true")

  if (file.exists(cache_path) && !refresh) {
    cached <- readRDS(cache_path)
    cached$cache_path <- cache_path
    cached$from_cache <- TRUE
    return(cached)
  }
  if (skip) {
    return(list(
      success = FALSE,
      error = "Bayesian fitting skipped because TRUNCCOMP_SKIP_LICORICE_BAYES=true.",
      cache_path = cache_path,
      from_cache = FALSE
    ))
  }

  if (!requireNamespace("rstan", quietly = TRUE)) {
    return(list(
      success = FALSE,
      error = "The rstan package is not available.",
      cache_path = cache_path,
      from_cache = FALSE
    ))
  }

  settings <- licorice_appendix_bayes_settings()
  options(mc.cores = min(settings$chains, max(1L, parallel::detectCores(logical = FALSE))))
  rstan::rstan_options(auto_write = TRUE)

  fit <- tryCatch(
    TruncComp2::trunc_comp_bayes(
      Y ~ R,
      atom = metadata$atom,
      data = d,
      continuous_support = settings$continuous_support,
      mixture_components = settings$mixture_components,
      auto_select_mixture_components = settings$auto_select_mixture_components,
      chains = settings$chains,
      iter_warmup = settings$iter_warmup,
      iter_sampling = settings$iter_sampling,
      seed = settings$seed,
      refresh = 0,
      control = settings$control,
      prior = settings$prior
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    out <- list(
      success = FALSE,
      error = conditionMessage(fit),
      settings = settings,
      cache_path = cache_path,
      from_cache = FALSE
    )
    saveRDS(out, cache_path)
    return(out)
  }

  ppc <- tryCatch(
    TruncComp2::posterior_predictive_pvalues(fit, ndraws = 200, seed = settings$seed + 1L),
    error = function(e) e
  )
  draws <- fit$draws
  probabilities <- if (isTRUE(fit$success)) {
    c(
      delta_lt_0 = mean(draws$delta < 0),
      mu_delta_lt_0 = mean(draws$mu_delta < 0),
      alpha_delta_lt_1 = mean(draws$alpha_delta < 1),
      no_pain_probability_increase = mean(draws$delta_atom > 0),
      any_pain_probability_reduction = mean((draws$pi_1 - draws$pi_0) < 0)
    )
  } else {
    c(
      delta_lt_0 = NA_real_,
      mu_delta_lt_0 = NA_real_,
      alpha_delta_lt_1 = NA_real_,
      no_pain_probability_increase = NA_real_,
      any_pain_probability_reduction = NA_real_
    )
  }

  out <- list(
    success = isTRUE(fit$success),
    fit = fit,
    summary_table = fit$summary_table,
    arm_table = fit$arm_table,
    diagnostics = fit$diagnostics,
    ppc_table = if (inherits(ppc, "error")) NULL else ppc,
    ppc_error = if (inherits(ppc, "error")) conditionMessage(ppc) else NULL,
    probabilities = probabilities,
    settings = settings,
    cache_path = cache_path,
    from_cache = FALSE
  )
  saveRDS(out, cache_path)
  out
}

analyze_licorice_appendix_data <- function(licorice_data, manuscript_dir) {
  d <- licorice_data$data
  metadata <- licorice_data$metadata
  atom <- metadata$atom

  model_lrt <- TruncComp2::trunc_comp(Y ~ R, atom = atom, data = d, method = "lrt")
  model_splrt <- TruncComp2::trunc_comp(Y ~ R, atom = atom, data = d, method = "splrt")
  surface <- suppressMessages(stats::confint(
    model_splrt,
    parameter = "joint",
    resolution = metadata$surface_resolution,
    plot = FALSE
  ))
  delta_intervals <- application_delta_intervals(model_splrt, resolution = metadata$surface_resolution)
  summary_stats <- licorice_appendix_summary_stats(d, metadata)
  contrasts <- licorice_appendix_contrasts(summary_stats)
  standard_tests <- run_standard_licorice_appendix_tests(d)
  bayes <- fit_or_load_licorice_appendix_bayes(d, metadata, manuscript_dir)

  list(
    data = d,
    raw = licorice_data$raw,
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

compute_licorice_appendix_results <- function(manuscript_dir) {
  analyze_licorice_appendix_data(load_licorice_appendix_data(manuscript_dir), manuscript_dir)
}
