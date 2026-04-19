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
      "Survivor EQ-VAS mean difference",
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
      "rt-PA minus control mean EQ-VAS among 6-month survivors",
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
      "Survivors only",
      "Survivors only",
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

  list(
    continuous_support = "real_line",
    mixture_components = as_int("TRUNCCOMP_IST3_BAYES_COMPONENTS", 2L),
    auto_select_mixture_components = FALSE,
    chains = as_int("TRUNCCOMP_IST3_BAYES_CHAINS", 4L),
    iter_warmup = as_int("TRUNCCOMP_IST3_BAYES_WARMUP", 300L),
    iter_sampling = as_int("TRUNCCOMP_IST3_BAYES_SAMPLING", 300L),
    seed = as_int("TRUNCCOMP_IST3_BAYES_SEED", 20260419L),
    control = list(adapt_delta = 0.99, max_treedepth = 12L),
    prior = list(
      rho_alpha = 1,
      rho_beta = 1,
      mu_mean = 0,
      mu_sd = 2.5,
      sigma_meanlog = -0.5,
      sigma_sdlog = 0.5,
      alpha_shape = 2,
      alpha_rate = 1
    )
  )
}

fit_or_load_application_bayes <- function(d, metadata, manuscript_dir) {
  cache_path <- application_bayes_cache_path(manuscript_dir)
  refresh <- identical(tolower(Sys.getenv("TRUNCCOMP_REFRESH_IST3_BAYES", "false")), "true")
  skip <- identical(tolower(Sys.getenv("TRUNCCOMP_SKIP_IST3_BAYES", "false")), "true")

  if (file.exists(cache_path) && !refresh) {
    cached <- readRDS(cache_path)
    cached$cache_path <- cache_path
    cached$from_cache <- TRUE
    return(cached)
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

  settings <- application_bayes_settings()
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
