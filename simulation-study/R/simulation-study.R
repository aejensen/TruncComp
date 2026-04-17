simulation_study_root <- function(repo_root) {
  ensure_dir(file.path(repo_root, "simulation-study"))
}

simulation_study_output_dir <- function(repo_root) {
  ensure_dir(file.path(simulation_study_root(repo_root), "results", "trunccomp2-study"))
}

simulation_study_results_path <- function(repo_root) {
  file.path(simulation_study_output_dir(repo_root), "simulation-study.rds")
}

simulation_study_config_path <- function(output_dir) {
  file.path(output_dir, "config.rds")
}

simulation_study_manifest_path <- function(output_dir) {
  file.path(output_dir, "pending-cells.rds")
}

simulation_study_cells_dir <- function(output_dir) {
  ensure_dir(file.path(output_dir, "cells"))
}

simulation_study_slurm_dir <- function(output_dir) {
  ensure_dir(file.path(output_dir, "slurm"))
}

simulation_study_methods <- function() {
  c("wilcoxon", "t_test", "lrt", "splrt")
}

simulation_study_method_labels <- function() {
  c(
    wilcoxon = "Wilcoxon",
    t_test = "T-test",
    lrt = "Parametric LRT",
    splrt = "Semi-parametric LRT"
  )
}

load_local_trunccomp2 <- function(repo_root) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The pkgload package is required to load the local TruncComp2 package.", call. = FALSE)
  }

  pkgload::load_all(
    file.path(repo_root, "packages", "TruncComp2"),
    quiet = TRUE,
    export_all = FALSE,
    helpers = FALSE,
    attach_testthat = FALSE
  )

  invisible(TRUE)
}

simulation_study_default_config <- function(reps = 10000L,
                                            n_seq = seq.int(50L, 500L, by = 50L),
                                            effect_levels = 0:3,
                                            alpha = 0.05,
                                            atom = 0,
                                            base_seed = 20260416L,
                                            power_effect_n = 250L) {
  list(
    version = "trunccomp2-simulation-study-v1",
    reps = as.integer(reps),
    n_seq = as.integer(n_seq),
    effect_levels = as.integer(effect_levels),
    alpha = alpha,
    atom = atom,
    base_seed = as.integer(base_seed),
    power_effect_n = as.integer(power_effect_n),
    main_text_scenarios = c("S2", "S3", "S5", "S6"),
    supplementary_scenarios = c("S1", "S4")
  )
}

.simulation_study_contaminated_shift <- function(n, shift) {
  major_component <- stats::rbinom(n, 1, 0.95) == 1
  means <- ifelse(major_component, 3 + shift, 10 + shift)
  sds <- ifelse(major_component, 0.35, 0.5)
  stats::rnorm(n, mean = means, sd = sds)
}

simulation_study_scenarios <- function() {
  list(
    S1 = list(
      scenario_id = "S1",
      scenario_order = 1L,
      short_label = "S1 Normal survivor shift",
      title = "Normal survivor shift, low atom",
      purpose = paste(
        "Baseline parametric reference where the survivor distribution is normal,",
        "the atom probability is unchanged, and both TruncComp2 methods should track each other closely."
      ),
      atom_latex = "$\\pi_0 = \\pi_1 = 0.85$",
      survivor_latex = paste(
        "$Y \\mid A = 1, R = 0 \\sim N(3, 1)$,",
        "$Y \\mid A = 1, R = 1 \\sim N(3 + \\delta_h, 1)$,",
        "$\\delta_h \\in \\{0, 0.20, 0.35, 0.50\\}$"
      ),
      effect_parameters = function(h) {
        delta_values <- c(0, 0.20, 0.35, 0.50)
        delta <- delta_values[[h + 1L]]
        list(
          pi0 = 0.85,
          pi1 = 0.85,
          f0 = function(n) stats::rnorm(n, 3, 1),
          f1 = function(n) stats::rnorm(n, 3 + delta, 1),
          effect_text = sprintf("$\\\\delta = %.2f$", delta),
          survivor_mean0 = 3,
          survivor_mean1 = 3 + delta,
          expected_combined_delta = 0.85 * delta
        )
      }
    ),
    S2 = list(
      scenario_id = "S2",
      scenario_order = 2L,
      short_label = "S2 Survival-only effect",
      title = "Survival-only effect",
      purpose = paste(
        "Two-degree-of-freedom penalty case where treatment acts only through the binary observation component."
      ),
      atom_latex = "$\\pi_0 = 0.45$, $\\pi_1 \\in \\{0.45, 0.50, 0.55, 0.60\\}$",
      survivor_latex = "$Y \\mid A = 1, R = 0 \\sim Y \\mid A = 1, R = 1 \\sim N(3, 1)$",
      effect_parameters = function(h) {
        pi1_values <- c(0.45, 0.50, 0.55, 0.60)
        pi1 <- pi1_values[[h + 1L]]
        list(
          pi0 = 0.45,
          pi1 = pi1,
          f0 = function(n) stats::rnorm(n, 3, 1),
          f1 = function(n) stats::rnorm(n, 3, 1),
          effect_text = sprintf("$\\\\pi_1 = %.2f$", pi1),
          survivor_mean0 = 3,
          survivor_mean1 = 3,
          expected_combined_delta = (pi1 - 0.45) * 3
        )
      }
    ),
    S3 = list(
      scenario_id = "S3",
      scenario_order = 3L,
      short_label = "S3 Antagonistic cancellation",
      title = "Antagonistic cancellation",
      purpose = paste(
        "Signature case where survival and survivor means move in opposite directions while the combined-outcome mean stays constant."
      ),
      atom_latex = "$\\pi_0 = 0.50$, $\\pi_1 \\in \\{0.50, 0.55, 0.60, 0.65\\}$",
      survivor_latex = paste(
        "$Y \\mid A = 1, R = 0 \\sim N(3, 1)$,",
        "$Y \\mid A = 1, R = 1 \\sim N(1.5 / \\pi_1, 1)$"
      ),
      effect_parameters = function(h) {
        pi1_values <- c(0.50, 0.55, 0.60, 0.65)
        pi1 <- pi1_values[[h + 1L]]
        treatment_mean <- 1.5 / pi1
        list(
          pi0 = 0.50,
          pi1 = pi1,
          f0 = function(n) stats::rnorm(n, 3, 1),
          f1 = function(n) stats::rnorm(n, treatment_mean, 1),
          effect_text = sprintf("$\\\\pi_1 = %.2f$, $E[Y \\\\mid A = 1, R = 1] = %.3f$", pi1, treatment_mean),
          survivor_mean0 = 3,
          survivor_mean1 = treatment_mean,
          expected_combined_delta = 0
        )
      }
    ),
    S4 = list(
      scenario_id = "S4",
      scenario_order = 4L,
      short_label = "S4 Concordant joint improvement",
      title = "Concordant joint improvement",
      purpose = paste(
        "Configuration where both components move in the same helpful direction, making one-dimensional tests more competitive."
      ),
      atom_latex = "$\\pi_0 = 0.50$, $\\pi_1 \\in \\{0.50, 0.55, 0.60, 0.65\\}$",
      survivor_latex = paste(
        "$Y \\mid A = 1, R = 0 \\sim N(3, 1)$,",
        "$Y \\mid A = 1, R = 1 \\sim N(3 + 0.15h, 1)$"
      ),
      effect_parameters = function(h) {
        pi1_values <- c(0.50, 0.55, 0.60, 0.65)
        pi1 <- pi1_values[[h + 1L]]
        treatment_mean <- 3 + 0.15 * h
        list(
          pi0 = 0.50,
          pi1 = pi1,
          f0 = function(n) stats::rnorm(n, 3, 1),
          f1 = function(n) stats::rnorm(n, treatment_mean, 1),
          effect_text = sprintf("$\\\\pi_1 = %.2f$, $E[Y \\\\mid A = 1, R = 1] = %.2f$", pi1, treatment_mean),
          survivor_mean0 = 3,
          survivor_mean1 = treatment_mean,
          expected_combined_delta = pi1 * treatment_mean - 1.5
        )
      }
    ),
    S5 = list(
      scenario_id = "S5",
      scenario_order = 5L,
      short_label = "S5 Gamma shape change",
      title = "Gamma shape-change with mean shift",
      purpose = paste(
        "Non-normal survivor case with concurrent mean, variance, and skewness changes designed to separate semi-parametric and parametric performance."
      ),
      atom_latex = "$\\pi_0 = \\pi_1 = 0.65$",
      survivor_latex = paste(
        "$Y \\mid A = 1, R = 0 \\sim \\Gamma(9, 3/9)$,",
        "$Y \\mid A = 1, R = 1 \\sim \\Gamma(k_h, m_h / k_h)$,",
        "$(k_h, m_h) \\in \\{(9, 3.0), (6, 3.2), (4, 3.4), (2.5, 3.6)\\}$"
      ),
      effect_parameters = function(h) {
        shapes <- c(9, 6, 4, 2.5)
        means <- c(3.0, 3.2, 3.4, 3.6)
        shape1 <- shapes[[h + 1L]]
        mean1 <- means[[h + 1L]]
        list(
          pi0 = 0.65,
          pi1 = 0.65,
          f0 = function(n) stats::rgamma(n, shape = 9, scale = 3 / 9),
          f1 = function(n) stats::rgamma(n, shape = shape1, scale = mean1 / shape1),
          effect_text = sprintf("$k_1 = %.1f$, $E[Y \\\\mid A = 1, R = 1] = %.1f$", shape1, mean1),
          survivor_mean0 = 3,
          survivor_mean1 = mean1,
          expected_combined_delta = 0.65 * (mean1 - 3)
        )
      }
    ),
    S6 = list(
      scenario_id = "S6",
      scenario_order = 6L,
      short_label = "S6 Rare-outlier contamination",
      title = "Rare-outlier contamination",
      purpose = paste(
        "Mixture case with rare extreme survivors that is deliberately friendly to rank-based robustness."
      ),
      atom_latex = "$\\pi_0 = \\pi_1 = 0.80$",
      survivor_latex = paste(
        "$Y \\mid A = 1, R = r$ is",
        "$0.95 N(3 + \\delta_h r, 0.35^2) + 0.05 N(10 + \\delta_h r, 0.5^2)$,",
        "$\\delta_h \\in \\{0, 0.10, 0.20, 0.35\\}$"
      ),
      effect_parameters = function(h) {
        delta_values <- c(0, 0.10, 0.20, 0.35)
        delta <- delta_values[[h + 1L]]
        list(
          pi0 = 0.80,
          pi1 = 0.80,
          f0 = function(n) .simulation_study_contaminated_shift(n, shift = 0),
          f1 = function(n) .simulation_study_contaminated_shift(n, shift = delta),
          effect_text = sprintf("$\\\\delta = %.2f$", delta),
          survivor_mean0 = 0.95 * 3 + 0.05 * 10,
          survivor_mean1 = 0.95 * (3 + delta) + 0.05 * (10 + delta),
          expected_combined_delta = 0.80 * delta
        )
      }
    )
  )
}

simulation_study_scenario_table <- function(config = simulation_study_default_config()) {
  scenario_ids <- names(simulation_study_scenarios())
  rows <- lapply(scenario_ids, function(scenario_id) {
    scenario <- simulation_study_scenarios()[[scenario_id]]
    data.frame(
      scenario_id = scenario$scenario_id,
      scenario_order = scenario$scenario_order,
      short_label = scenario$short_label,
      title = scenario$title,
      purpose = scenario$purpose,
      atom_latex = scenario$atom_latex,
      survivor_latex = scenario$survivor_latex,
      in_main_text = scenario$scenario_id %in% config$main_text_scenarios,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

simulation_study_design <- function(config = simulation_study_default_config()) {
  scenarios <- simulation_study_scenarios()
  design_rows <- list()
  cell_id <- 1L

  for (scenario_id in names(scenarios)) {
    scenario <- scenarios[[scenario_id]]
    for (h in config$effect_levels) {
      params <- scenario$effect_parameters(h)
      for (n_value in config$n_seq) {
        design_rows[[length(design_rows) + 1L]] <- data.frame(
          cell_id = cell_id,
          scenario_id = scenario$scenario_id,
          scenario_order = scenario$scenario_order,
          scenario_label = paste0(scenario$scenario_id, ": ", scenario$title),
          short_label = scenario$short_label,
          purpose = scenario$purpose,
          in_main_text = scenario$scenario_id %in% config$main_text_scenarios,
          n = as.integer(n_value),
          h = as.integer(h),
          is_null = h == 0L,
          reps = config$reps,
          alpha = config$alpha,
          atom = config$atom,
          seed = as.integer(config$base_seed + cell_id - 1L),
          pi0 = params$pi0,
          pi1 = params$pi1,
          effect_text = params$effect_text,
          survivor_mean0 = params$survivor_mean0,
          survivor_mean1 = params$survivor_mean1,
          expected_combined_delta = params$expected_combined_delta,
          stringsAsFactors = FALSE
        )
        cell_id <- cell_id + 1L
      }
    }
  }

  do.call(rbind, design_rows)
}

simulation_study_cell_path <- function(output_dir, cell) {
  file.path(
    simulation_study_cells_dir(output_dir),
    sprintf("%s-h%d-n%03d.rds", cell$scenario_id, cell$h, cell$n)
  )
}

simulation_study_prepare_run <- function(output_dir,
                                         reps = 10000L,
                                         scenario_ids = NULL,
                                         n_seq = NULL,
                                         effect_levels = NULL,
                                         overwrite = FALSE) {
  config <- simulation_study_default_config(reps = reps)
  if (!is.null(n_seq)) {
    config$n_seq <- as.integer(n_seq)
  }
  if (!is.null(effect_levels)) {
    config$effect_levels <- as.integer(effect_levels)
  }

  output_dir <- ensure_dir(output_dir)
  saveRDS(config, simulation_study_config_path(output_dir))

  design <- simulation_study_design(config)
  if (!is.null(scenario_ids)) {
    design <- design[design$scenario_id %in% scenario_ids, , drop = FALSE]
  }

  if (!nrow(design)) {
    stop("No simulation-study cells matched the requested filters.", call. = FALSE)
  }

  target_paths <- vapply(
    split(design, seq_len(nrow(design))),
    function(cell) simulation_study_cell_path(output_dir, cell),
    character(1)
  )

  pending_design <- design
  if (!overwrite) {
    pending_design <- pending_design[!file.exists(target_paths), , drop = FALSE]
  }

  list(
    config = config,
    output_dir = output_dir,
    design = design,
    pending_design = pending_design,
    target_paths = target_paths
  )
}

.simulation_study_measure_pvalue <- function(expr) {
  start <- proc.time()[["elapsed"]]
  p_value <- tryCatch(
    suppressWarnings(force(expr)),
    error = function(e) NA_real_
  )
  elapsed <- proc.time()[["elapsed"]] - start

  if (!(length(p_value) == 1L && is.numeric(p_value) && is.finite(p_value))) {
    p_value <- NA_real_
  }

  list(p_value = as.numeric(p_value), elapsed = as.numeric(elapsed))
}

.simulation_study_evaluate_methods <- function(data, atom) {
  results <- vector("list", length(simulation_study_methods()))
  names(results) <- simulation_study_methods()

  results$wilcoxon <- .simulation_study_measure_pvalue(
    stats::wilcox.test(Y ~ R, data = data, exact = FALSE)$p.value
  )
  results$t_test <- .simulation_study_measure_pvalue(
    stats::t.test(Y ~ R, data = data)$p.value
  )
  results$lrt <- .simulation_study_measure_pvalue(
    TruncComp2::truncComp(Y ~ R, atom = atom, data = data, method = "LRT")$p
  )
  results$splrt <- .simulation_study_measure_pvalue(
    TruncComp2::truncComp(Y ~ R, atom = atom, data = data, method = "SPLRT")$p
  )

  results
}

.simulation_study_write_rds_atomic <- function(object, path) {
  tmp_path <- tempfile(pattern = "trunccomp2-sim-", tmpdir = dirname(path), fileext = ".rds")
  saveRDS(object, file = tmp_path)

  if (!file.rename(tmp_path, path)) {
    file.copy(tmp_path, path, overwrite = TRUE)
    unlink(tmp_path)
  }

  invisible(path)
}

simulation_study_write_manifest <- function(design, output_dir) {
  .simulation_study_write_rds_atomic(design, simulation_study_manifest_path(output_dir))
}

simulation_study_run_cell <- function(cell, output_dir) {
  cell <- cell[1, , drop = FALSE]
  scenario <- simulation_study_scenarios()[[cell$scenario_id]]
  params <- scenario$effect_parameters(cell$h)
  method_names <- simulation_study_methods()

  set.seed(cell$seed)

  p_values <- matrix(NA_real_, nrow = cell$reps, ncol = length(method_names))
  colnames(p_values) <- method_names
  runtimes <- matrix(NA_real_, nrow = cell$reps, ncol = length(method_names))
  colnames(runtimes) <- method_names
  combined_deltas <- numeric(cell$reps)

  for (rep_index in seq_len(cell$reps)) {
    data <- TruncComp2::simulateTruncatedData(
      n = cell$n,
      f0 = params$f0,
      f1 = params$f1,
      pi0 = params$pi0,
      pi1 = params$pi1,
      atom = cell$atom
    )

    combined_deltas[[rep_index]] <- mean(data$Y[data$R == 1]) - mean(data$Y[data$R == 0])
    method_results <- .simulation_study_evaluate_methods(data, atom = cell$atom)

    for (method_name in method_names) {
      p_values[rep_index, method_name] <- method_results[[method_name]]$p_value
      runtimes[rep_index, method_name] <- method_results[[method_name]]$elapsed
    }
  }

  method_metrics <- lapply(method_names, function(method_name) {
    valid <- is.finite(p_values[, method_name])
    successes <- sum(valid)
    reject_rate <- if (successes) {
      mean(p_values[valid, method_name] < cell$alpha)
    } else {
      NA_real_
    }

    data.frame(
      method = method_name,
      reject_rate = reject_rate,
      mcse = if (successes) sqrt(reject_rate * (1 - reject_rate) / successes) else NA_real_,
      fail_rate = mean(!valid),
      median_runtime_sec = if (any(is.finite(runtimes[, method_name]))) {
        stats::median(runtimes[, method_name], na.rm = TRUE)
      } else {
        NA_real_
      },
      successful_reps = successes,
      stringsAsFactors = FALSE
    )
  })
  method_metrics <- do.call(rbind, method_metrics)

  cell_summary <- data.frame(
    cell_id = cell$cell_id,
    scenario_id = cell$scenario_id,
    scenario_order = cell$scenario_order,
    scenario_label = cell$scenario_label,
    short_label = cell$short_label,
    purpose = cell$purpose,
    in_main_text = cell$in_main_text,
    n = cell$n,
    h = cell$h,
    is_null = cell$is_null,
    reps = cell$reps,
    alpha = cell$alpha,
    atom = cell$atom,
    seed = cell$seed,
    effect_text = cell$effect_text,
    pi0 = cell$pi0,
    pi1 = cell$pi1,
    expected_combined_delta = cell$expected_combined_delta,
    empirical_combined_delta = mean(combined_deltas),
    empirical_combined_delta_sd = stats::sd(combined_deltas),
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE
  )

  result <- list(
    version = "trunccomp2-simulation-study-cell-v1",
    cell_summary = cell_summary,
    method_metrics = method_metrics
  )

  .simulation_study_write_rds_atomic(result, simulation_study_cell_path(output_dir, cell))
  result
}

aggregate_simulation_study_results <- function(output_dir, config = NULL) {
  if (is.null(config)) {
    config_path <- simulation_study_config_path(output_dir)
    if (!file.exists(config_path)) {
      stop("No simulation-study config found. Run the simulation driver first.", call. = FALSE)
    }
    config <- readRDS(config_path)
  }

  design <- simulation_study_design(config)
  cell_paths <- vapply(
    split(design, seq_len(nrow(design))),
    function(cell) simulation_study_cell_path(output_dir, cell),
    character(1)
  )
  existing <- file.exists(cell_paths)

  if (!any(existing)) {
    stop("No simulation-study cell files were found. Run the simulation driver first.", call. = FALSE)
  }

  cell_results <- lapply(cell_paths[existing], readRDS)
  cell_metrics <- do.call(rbind, lapply(cell_results, function(x) x$cell_summary))
  method_metrics <- do.call(rbind, lapply(cell_results, function(x) {
    cbind(
      x$cell_summary[rep(1, nrow(x$method_metrics)), c(
        "cell_id", "scenario_id", "scenario_order", "scenario_label", "short_label",
        "purpose", "in_main_text", "n", "h", "is_null", "reps", "alpha", "atom",
        "seed", "effect_text", "pi0", "pi1", "expected_combined_delta",
        "empirical_combined_delta", "empirical_combined_delta_sd"
      )],
      x$method_metrics,
      stringsAsFactors = FALSE
    )
  }))

  rownames(cell_metrics) <- NULL
  rownames(method_metrics) <- NULL
  method_metrics$method_label <- unname(simulation_study_method_labels()[method_metrics$method])
  method_metrics$method_label <- factor(
    method_metrics$method_label,
    levels = unname(simulation_study_method_labels())
  )

  results <- list(
    version = config$version,
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    config = config,
    complete = all(existing),
    expected_cells = nrow(design),
    completed_cells = sum(existing),
    scenarios = simulation_study_scenario_table(config),
    design = design,
    cell_metrics = cell_metrics[order(cell_metrics$scenario_order, cell_metrics$h, cell_metrics$n), , drop = FALSE],
    method_metrics = method_metrics[order(method_metrics$scenario_order, method_metrics$h, method_metrics$n, method_metrics$method_label), , drop = FALSE]
  )

  saveRDS(results, file = file.path(output_dir, "simulation-study.rds"))
  utils::write.csv(
    results$cell_metrics,
    file = file.path(output_dir, "simulation-study-cell-metrics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    results$method_metrics,
    file = file.path(output_dir, "simulation-study-method-metrics.csv"),
    row.names = FALSE
  )

  results
}

run_simulation_study <- function(repo_root,
                                 output_dir,
                                 reps = 10000L,
                                 workers = 1L,
                                 scenario_ids = NULL,
                                 n_seq = NULL,
                                 effect_levels = NULL,
                                 overwrite = FALSE) {
  prepared <- simulation_study_prepare_run(
    output_dir = output_dir,
    reps = reps,
    scenario_ids = scenario_ids,
    n_seq = n_seq,
    effect_levels = effect_levels,
    overwrite = overwrite
  )

  load_local_trunccomp2(repo_root)

  if (nrow(prepared$pending_design)) {
    cell_tasks <- split(prepared$pending_design, seq_len(nrow(prepared$pending_design)))
    worker_fun <- function(cell) {
      message(sprintf("Running %s, h = %d, n = %d", cell$scenario_id, cell$h, cell$n))
      simulation_study_run_cell(cell, output_dir = prepared$output_dir)
    }

    if (.Platform$OS.type != "windows" && workers > 1L) {
      parallel::mclapply(cell_tasks, worker_fun, mc.cores = workers, mc.preschedule = FALSE)
    } else {
      lapply(cell_tasks, worker_fun)
    }
  }

  aggregate_simulation_study_results(prepared$output_dir, config = prepared$config)
}

simulation_study_metric_wide <- function(simulation_results, value_col) {
  method_metrics <- simulation_results$method_metrics
  wide <- stats::reshape(
    method_metrics[, c("cell_id", "method", value_col)],
    idvar = "cell_id",
    timevar = "method",
    direction = "wide"
  )
  names(wide) <- sub(paste0("^", value_col, "\\."), "", names(wide))

  meta <- unique(method_metrics[, c(
    "cell_id", "scenario_id", "scenario_order", "scenario_label", "short_label",
    "n", "h", "is_null", "effect_text", "expected_combined_delta", "empirical_combined_delta"
  )])

  out <- merge(meta, wide, by = "cell_id", sort = FALSE)
  out[order(out$scenario_order, out$h, out$n), , drop = FALSE]
}

simulation_study_null_flags <- function(simulation_results, lower = 0.04, upper = 0.06, fail_cutoff = 0.01) {
  reject_wide <- simulation_study_metric_wide(simulation_results, "reject_rate")
  fail_wide <- simulation_study_metric_wide(simulation_results, "fail_rate")

  reject_wide <- reject_wide[reject_wide$is_null, , drop = FALSE]
  fail_wide <- fail_wide[fail_wide$is_null, , drop = FALSE]
  merged <- merge(
    reject_wide,
    fail_wide[, c("cell_id", simulation_study_methods())],
    by = "cell_id",
    suffixes = c("_reject", "_fail"),
    sort = FALSE
  )

  merged$size_flag <- apply(
    merged[, paste0(simulation_study_methods(), "_reject"), drop = FALSE],
    1,
    function(x) any(is.finite(x) & (x < lower | x > upper))
  )
  merged$fail_flag <- apply(
    merged[, paste0(simulation_study_methods(), "_fail"), drop = FALSE],
    1,
    function(x) any(is.finite(x) & x > fail_cutoff)
  )

  merged
}

simulation_study_finalize_results <- function(repo_root,
                                              output_dir,
                                              simulation_results,
                                              refresh_manuscript_assets = TRUE) {
  if (refresh_manuscript_assets && isTRUE(simulation_results$complete)) {
    simulation_study_build_manuscript_assets(repo_root, simulation_results)
    message("Simulation-study manuscript assets written to manuscript/build.")
  } else if (refresh_manuscript_assets) {
    message(
      sprintf(
        "Simulation study is incomplete (%d/%d cells), so manuscript assets were not refreshed.",
        simulation_results$completed_cells,
        simulation_results$expected_cells
      )
    )
  }

  message(
    sprintf(
      "Simulation study complete: %d/%d cells available. Aggregated results written to %s",
      simulation_results$completed_cells,
      simulation_results$expected_cells,
      file.path(output_dir, "simulation-study.rds")
    )
  )

  invisible(simulation_results)
}
