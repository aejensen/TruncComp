ci_study_root <- function(repo_root) {
  simulation_study_root(repo_root)
}

ci_study_output_dir <- function(repo_root) {
  ensure_dir(file.path(ci_study_root(repo_root), "results", "trunccomp2-ci-study"))
}

ci_study_config_path <- function(output_dir) {
  file.path(output_dir, "config.rds")
}

ci_study_manifest_path <- function(output_dir) {
  file.path(output_dir, "pending-chunks.rds")
}

ci_study_chunks_dir <- function(output_dir) {
  ensure_dir(file.path(output_dir, "cells"))
}

ci_study_slurm_dir <- function(output_dir) {
  ensure_dir(file.path(output_dir, "slurm"))
}

ci_study_slurm_submission_path <- function(output_dir) {
  file.path(ci_study_slurm_dir(output_dir), "latest-submission.txt")
}

ci_study_existing_chunk_paths <- function(output_dir) {
  chunks_dir <- file.path(output_dir, "cells")
  if (!dir.exists(chunks_dir)) {
    return(character())
  }

  sort(list.files(chunks_dir, pattern = "\\.rds$", full.names = TRUE))
}

ci_study_methods <- function() {
  c("lrt", "splrt")
}

ci_study_method_labels <- function() {
  c(
    lrt = "Parametric LRT",
    splrt = "Semi-parametric LRT"
  )
}

ci_study_interval_targets <- function() {
  data.frame(
    interval_target = c(
      "mu_delta",
      "alpha_delta",
      "log_or_delta",
      "joint_region",
      "delta_welch",
      "delta_projected_optimize",
      "delta_profile_optimize"
    ),
    interval_label = c(
      "mu_delta",
      "alpha_delta",
      "log_or_delta",
      "Joint region",
      "delta Welch",
      "delta projected (optimize)",
      "delta profile (optimize)"
    ),
    interval_order = seq_len(7L),
    stringsAsFactors = FALSE
  )
}

ci_study_default_config <- function(reps = 10000L,
                                    chunk_reps = 500L,
                                    scenario_ids = c("S1", "S2", "S3", "S5", "S6"),
                                    n_seq = c(50L, 100L, 150L, 500L),
                                    effect_levels = c(0L, 2L, 3L),
                                    conf_level = 0.95,
                                    atom = 0,
                                    base_seed = 20260701L,
                                    cell_seed_stride = 1000000L,
                                    joint_resolution = 35L) {
  reps <- as.integer(reps)
  chunk_reps <- as.integer(chunk_reps)
  n_seq <- as.integer(n_seq)
  effect_levels <- as.integer(effect_levels)
  joint_resolution <- as.integer(joint_resolution)
  cell_seed_stride <- as.integer(cell_seed_stride)

  if (!(length(reps) == 1L && is.finite(reps) && reps >= 1L)) {
    stop("reps must be a single positive integer.", call. = FALSE)
  }
  if (!(length(chunk_reps) == 1L && is.finite(chunk_reps) && chunk_reps >= 1L)) {
    stop("chunk_reps must be a single positive integer.", call. = FALSE)
  }
  if (!(length(joint_resolution) == 1L && is.finite(joint_resolution) && joint_resolution >= 2L)) {
    stop("joint_resolution must be a single integer of at least 2.", call. = FALSE)
  }
  if (!(length(conf_level) == 1L && is.numeric(conf_level) && is.finite(conf_level) &&
        conf_level > 0 && conf_level < 1)) {
    stop("conf_level must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  if (!length(scenario_ids)) {
    stop("scenario_ids must contain at least one scenario.", call. = FALSE)
  }
  if (!length(n_seq)) {
    stop("n_seq must contain at least one sample size.", call. = FALSE)
  }
  if (!length(effect_levels)) {
    stop("effect_levels must contain at least one effect level.", call. = FALSE)
  }

  list(
    version = "trunccomp2-ci-study-v1",
    cell_version = "trunccomp2-ci-study-cell-v1",
    reps = reps,
    chunk_reps = chunk_reps,
    scenario_ids = as.character(scenario_ids),
    n_seq = n_seq,
    effect_levels = effect_levels,
    conf_level = as.numeric(conf_level),
    atom = atom,
    base_seed = as.integer(base_seed),
    cell_seed_stride = cell_seed_stride,
    joint_resolution = joint_resolution,
    methods = ci_study_methods(),
    interval_targets = ci_study_interval_targets()$interval_target
  )
}

ci_study_scenario_table <- function(config = ci_study_default_config()) {
  scenarios <- simulation_study_scenarios()
  scenario_ids <- config$scenario_ids
  missing <- setdiff(scenario_ids, names(scenarios))
  if (length(missing)) {
    stop("Unknown CI-study scenario id(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  rows <- lapply(scenario_ids, function(scenario_id) {
    scenario <- scenarios[[scenario_id]]
    data.frame(
      scenario_id = scenario$scenario_id,
      scenario_order = scenario$scenario_order,
      short_label = scenario$short_label,
      title = scenario$title,
      purpose = scenario$purpose,
      atom_latex = scenario$atom_latex,
      survivor_latex = scenario$survivor_latex,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

ci_study_design <- function(config = ci_study_default_config()) {
  scenarios <- simulation_study_scenarios()
  missing <- setdiff(config$scenario_ids, names(scenarios))
  if (length(missing)) {
    stop("Unknown CI-study scenario id(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  design_rows <- list()
  cell_id <- 1L

  for (scenario_id in config$scenario_ids) {
    scenario <- scenarios[[scenario_id]]
    for (h in config$effect_levels) {
      params <- scenario$effect_parameters(h)
      true_mu_delta <- params$survivor_mean1 - params$survivor_mean0
      true_log_or_delta <- stats::qlogis(params$pi1) - stats::qlogis(params$pi0)
      true_alpha_delta <- exp(true_log_or_delta)
      true_delta <- params$pi1 * params$survivor_mean1 +
        (1 - params$pi1) * config$atom -
        (params$pi0 * params$survivor_mean0 + (1 - params$pi0) * config$atom)

      for (n_value in config$n_seq) {
        design_rows[[length(design_rows) + 1L]] <- data.frame(
          cell_id = cell_id,
          scenario_id = scenario$scenario_id,
          scenario_order = scenario$scenario_order,
          scenario_label = paste0(scenario$scenario_id, ": ", scenario$title),
          short_label = scenario$short_label,
          purpose = scenario$purpose,
          n = as.integer(n_value),
          h = as.integer(h),
          is_null = h == 0L,
          reps = config$reps,
          chunk_reps = config$chunk_reps,
          expected_chunks = as.integer(ceiling(config$reps / config$chunk_reps)),
          conf_level = config$conf_level,
          atom = config$atom,
          base_seed = config$base_seed,
          seed = as.integer(config$base_seed + (cell_id - 1L) * config$cell_seed_stride),
          joint_resolution = config$joint_resolution,
          pi0 = params$pi0,
          pi1 = params$pi1,
          effect_text = params$effect_text,
          survivor_mean0 = params$survivor_mean0,
          survivor_mean1 = params$survivor_mean1,
          true_mu_delta = true_mu_delta,
          true_log_or_delta = true_log_or_delta,
          true_alpha_delta = true_alpha_delta,
          true_delta = true_delta,
          stringsAsFactors = FALSE
        )
        cell_id <- cell_id + 1L
      }
    }
  }

  do.call(rbind, design_rows)
}

ci_study_chunk_manifest <- function(design, config = ci_study_default_config()) {
  manifest_rows <- list()
  task_id <- 1L

  for (row_index in seq_len(nrow(design))) {
    cell <- design[row_index, , drop = FALSE]
    chunk_count <- as.integer(ceiling(cell$reps / config$chunk_reps))

    for (chunk_id in seq_len(chunk_count)) {
      rep_start <- as.integer((chunk_id - 1L) * config$chunk_reps + 1L)
      rep_end <- as.integer(min(cell$reps, chunk_id * config$chunk_reps))
      manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
        task_id = task_id,
        cell,
        chunk_id = as.integer(chunk_id),
        rep_start = rep_start,
        rep_end = rep_end,
        chunk_reps_actual = as.integer(rep_end - rep_start + 1L),
        stringsAsFactors = FALSE
      )
      task_id <- task_id + 1L
    }
  }

  do.call(rbind, manifest_rows)
}

ci_study_chunk_path <- function(output_dir, chunk) {
  file.path(
    ci_study_chunks_dir(output_dir),
    sprintf(
      "%s-h%d-n%03d-chunk%03d.rds",
      chunk$scenario_id,
      chunk$h,
      chunk$n,
      chunk$chunk_id
    )
  )
}

ci_study_scalar_equal <- function(x, y) {
  if (is.numeric(x) || is.integer(x) || is.numeric(y) || is.integer(y)) {
    return(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 0)))
  }

  identical(as.character(x), as.character(y))
}

ci_study_chunk_file_matches <- function(path, chunk) {
  if (!file.exists(path)) {
    return(FALSE)
  }

  existing <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(existing$chunk_summary) || nrow(existing$chunk_summary) != 1L) {
    return(FALSE)
  }

  summary <- existing$chunk_summary[1, , drop = FALSE]
  fields <- c(
    "cell_id", "scenario_id", "n", "h", "reps", "chunk_reps",
    "chunk_id", "rep_start", "rep_end", "chunk_reps_actual",
    "conf_level", "atom", "seed", "joint_resolution"
  )
  if (!all(fields %in% names(summary)) || !all(fields %in% names(chunk))) {
    return(FALSE)
  }

  all(vapply(fields, function(field) {
    ci_study_scalar_equal(summary[[field]][[1]], chunk[[field]][[1]])
  }, logical(1)))
}

ci_study_write_rds_atomic <- function(object, path) {
  tmp_path <- tempfile(pattern = "trunccomp2-ci-", tmpdir = dirname(path), fileext = ".rds")
  saveRDS(object, file = tmp_path)

  if (!file.rename(tmp_path, path)) {
    file.copy(tmp_path, path, overwrite = TRUE)
    unlink(tmp_path)
  }

  invisible(path)
}

ci_study_write_manifest <- function(manifest, output_dir) {
  ci_study_write_rds_atomic(manifest, ci_study_manifest_path(output_dir))
}

ci_study_prepare_run <- function(output_dir,
                                 reps = 10000L,
                                 chunk_reps = 500L,
                                 scenario_ids = NULL,
                                 n_seq = NULL,
                                 effect_levels = NULL,
                                 conf_level = 0.95,
                                 joint_resolution = 35L,
                                 overwrite = FALSE) {
  config <- ci_study_default_config(
    reps = reps,
    chunk_reps = chunk_reps,
    conf_level = conf_level,
    joint_resolution = joint_resolution
  )
  if (!is.null(scenario_ids)) {
    config$scenario_ids <- as.character(scenario_ids)
  }
  if (!is.null(n_seq)) {
    config$n_seq <- as.integer(n_seq)
  }
  if (!is.null(effect_levels)) {
    config$effect_levels <- as.integer(effect_levels)
  }

  output_dir <- ensure_dir(output_dir)
  saveRDS(config, ci_study_config_path(output_dir))

  design <- ci_study_design(config)
  manifest <- ci_study_chunk_manifest(design, config)
  if (!nrow(manifest)) {
    stop("No CI-study chunks matched the requested filters.", call. = FALSE)
  }

  target_paths <- vapply(
    split(manifest, seq_len(nrow(manifest))),
    function(chunk) ci_study_chunk_path(output_dir, chunk),
    character(1)
  )

  pending_manifest <- manifest
  if (!overwrite) {
    reusable <- vapply(seq_len(nrow(manifest)), function(index) {
      ci_study_chunk_file_matches(target_paths[[index]], manifest[index, , drop = FALSE])
    }, logical(1))
    pending_manifest <- pending_manifest[!reusable, , drop = FALSE]
  }

  list(
    config = config,
    output_dir = output_dir,
    design = design,
    manifest = manifest,
    pending_manifest = pending_manifest,
    target_paths = target_paths
  )
}

ci_study_condition_message <- function(error) {
  if (is.null(error)) {
    return(NA_character_)
  }

  conditionMessage(error)
}

ci_study_measure_expr <- function(expr) {
  start <- proc.time()[["elapsed"]]
  error <- NULL
  value <- tryCatch(
    force(expr),
    error = function(e) {
      error <<- e
      NULL
    }
  )
  elapsed <- proc.time()[["elapsed"]] - start

  list(
    value = value,
    elapsed = as.numeric(elapsed),
    error = error
  )
}

ci_study_silent_confint <- function(fit, ...) {
  out <- NULL
  invisible(utils::capture.output(
    out <- suppressMessages(suppressWarnings(stats::confint(fit, ...)))
  ))
  out
}

ci_study_fit_method <- function(data, method, atom, conf_level) {
  ci_study_measure_expr(
    suppressWarnings(TruncComp2::trunc_comp(
      Y ~ R,
      atom = atom,
      data = data,
      method = method,
      conf.level = conf_level
    ))
  )
}

ci_study_make_record <- function(rep_index,
                                 method,
                                 interval_target,
                                 fit_success,
                                 fit_runtime_sec,
                                 interval_success = FALSE,
                                 covered = NA,
                                 nonfinite_endpoint = NA,
                                 outside_grid = NA,
                                 width = NA_real_,
                                 interval_runtime_sec = NA_real_,
                                 joint_area = NA_real_,
                                 joint_mu_range = NA_real_,
                                 joint_log_or_range = NA_real_,
                                 error_class = "none",
                                 error_message = NA_character_) {
  data.frame(
    rep_index = as.integer(rep_index),
    method = method,
    interval_target = interval_target,
    fit_success = as.logical(fit_success),
    interval_success = as.logical(interval_success),
    covered = as.logical(covered),
    nonfinite_endpoint = as.logical(nonfinite_endpoint),
    outside_grid = as.logical(outside_grid),
    width = as.numeric(width),
    interval_runtime_sec = as.numeric(interval_runtime_sec),
    fit_runtime_sec = as.numeric(fit_runtime_sec),
    joint_area = as.numeric(joint_area),
    joint_mu_range = as.numeric(joint_mu_range),
    joint_log_or_range = as.numeric(joint_log_or_range),
    error_class = error_class,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

ci_study_confint_bounds <- function(value) {
  if (is.null(value)) {
    return(numeric())
  }

  if (is.matrix(value) || is.data.frame(value)) {
    if (nrow(value) < 1L || ncol(value) < 2L) {
      return(numeric())
    }
    return(as.numeric(value[1, 1:2]))
  }

  if (is.numeric(value) && length(value) >= 2L) {
    return(as.numeric(value[1:2]))
  }

  numeric()
}

ci_study_record_from_bounds <- function(rep_index,
                                        method,
                                        interval_target,
                                        truth,
                                        bounds,
                                        fit_runtime_sec,
                                        interval_runtime_sec,
                                        measure_error = NULL) {
  if (!is.null(measure_error)) {
    return(ci_study_make_record(
      rep_index = rep_index,
      method = method,
      interval_target = interval_target,
      fit_success = TRUE,
      fit_runtime_sec = fit_runtime_sec,
      interval_success = FALSE,
      interval_runtime_sec = interval_runtime_sec,
      error_class = "interval_error",
      error_message = ci_study_condition_message(measure_error)
    ))
  }

  if (!(length(bounds) >= 2L && is.numeric(bounds))) {
    return(ci_study_make_record(
      rep_index = rep_index,
      method = method,
      interval_target = interval_target,
      fit_success = TRUE,
      fit_runtime_sec = fit_runtime_sec,
      interval_success = FALSE,
      interval_runtime_sec = interval_runtime_sec,
      error_class = "invalid_interval",
      error_message = "Interval did not contain two numeric endpoints."
    ))
  }

  lower <- as.numeric(bounds[[1]])
  upper <- as.numeric(bounds[[2]])
  if (!is.finite(truth)) {
    return(ci_study_make_record(
      rep_index = rep_index,
      method = method,
      interval_target = interval_target,
      fit_success = TRUE,
      fit_runtime_sec = fit_runtime_sec,
      interval_success = FALSE,
      interval_runtime_sec = interval_runtime_sec,
      error_class = "invalid_truth",
      error_message = "The target truth was not finite."
    ))
  }

  if (!all(is.finite(c(lower, upper)))) {
    return(ci_study_make_record(
      rep_index = rep_index,
      method = method,
      interval_target = interval_target,
      fit_success = TRUE,
      fit_runtime_sec = fit_runtime_sec,
      interval_success = TRUE,
      covered = NA,
      nonfinite_endpoint = TRUE,
      interval_runtime_sec = interval_runtime_sec,
      error_class = "nonfinite_endpoint",
      error_message = "At least one interval endpoint was non-finite."
    ))
  }

  if (lower > upper) {
    return(ci_study_make_record(
      rep_index = rep_index,
      method = method,
      interval_target = interval_target,
      fit_success = TRUE,
      fit_runtime_sec = fit_runtime_sec,
      interval_success = FALSE,
      interval_runtime_sec = interval_runtime_sec,
      error_class = "invalid_endpoint_order",
      error_message = "The lower endpoint exceeded the upper endpoint."
    ))
  }

  ci_study_make_record(
    rep_index = rep_index,
    method = method,
    interval_target = interval_target,
    fit_success = TRUE,
    fit_runtime_sec = fit_runtime_sec,
    interval_success = TRUE,
    covered = lower <= truth && truth <= upper,
    nonfinite_endpoint = FALSE,
    width = upper - lower,
    interval_runtime_sec = interval_runtime_sec
  )
}

ci_study_component_records <- function(fit, method, truths, fit_runtime_sec, rep_index) {
  mu_record <- ci_study_record_from_bounds(
    rep_index = rep_index,
    method = method,
    interval_target = "mu_delta",
    truth = truths$true_mu_delta,
    bounds = fit$mu_delta_ci,
    fit_runtime_sec = fit_runtime_sec,
    interval_runtime_sec = 0
  )

  alpha_record <- ci_study_record_from_bounds(
    rep_index = rep_index,
    method = method,
    interval_target = "alpha_delta",
    truth = truths$true_alpha_delta,
    bounds = fit$alpha_delta_ci,
    fit_runtime_sec = fit_runtime_sec,
    interval_runtime_sec = 0
  )

  alpha_bounds <- fit$alpha_delta_ci
  log_bounds <- if (length(alpha_bounds) >= 2L && all(is.finite(alpha_bounds[1:2])) &&
                    all(alpha_bounds[1:2] > 0)) {
    log(as.numeric(alpha_bounds[1:2]))
  } else {
    c(NA_real_, NA_real_)
  }
  log_record <- ci_study_record_from_bounds(
    rep_index = rep_index,
    method = method,
    interval_target = "log_or_delta",
    truth = truths$true_log_or_delta,
    bounds = log_bounds,
    fit_runtime_sec = fit_runtime_sec,
    interval_runtime_sec = 0
  )

  list(mu_record, alpha_record, log_record)
}

ci_study_delta_record <- function(fit,
                                  method,
                                  interval_target,
                                  delta_method,
                                  truth,
                                  fit_runtime_sec,
                                  rep_index,
                                  algorithm = NULL) {
  measure <- if (is.null(algorithm)) {
    ci_study_measure_expr(
      ci_study_silent_confint(
        fit,
        parameter = "delta",
        method = delta_method,
        conf.level = fit$conf.level,
        plot = FALSE
      )
    )
  } else {
    ci_study_measure_expr(
      ci_study_silent_confint(
        fit,
        parameter = "delta",
        method = delta_method,
        algorithm = algorithm,
        conf.level = fit$conf.level,
        plot = FALSE
      )
    )
  }

  ci_study_record_from_bounds(
    rep_index = rep_index,
    method = method,
    interval_target = interval_target,
    truth = truth,
    bounds = ci_study_confint_bounds(measure$value),
    fit_runtime_sec = fit_runtime_sec,
    interval_runtime_sec = measure$elapsed,
    measure_error = measure$error
  )
}

ci_study_joint_size <- function(x, y, z, threshold) {
  accepted <- is.finite(z) & z <= threshold
  dx <- if (length(x) >= 2L) stats::median(abs(diff(x))) else NA_real_
  dy <- if (length(y) >= 2L) stats::median(abs(diff(y))) else NA_real_
  area <- if (is.finite(dx) && is.finite(dy)) {
    sum(accepted) * dx * dy
  } else {
    NA_real_
  }

  indices <- which(accepted, arr.ind = TRUE)
  if (!nrow(indices)) {
    return(list(
      area = as.numeric(area),
      mu_range = NA_real_,
      log_or_range = NA_real_
    ))
  }

  list(
    area = as.numeric(area),
    mu_range = diff(range(x[indices[, 1L]], finite = TRUE)),
    log_or_range = diff(range(y[indices[, 2L]], finite = TRUE))
  )
}

ci_study_joint_bracket <- function(values, point) {
  if (point == max(values)) {
    return(length(values) - 1L)
  }

  index <- findInterval(point, values)
  as.integer(max(1L, min(index, length(values) - 1L)))
}

ci_study_joint_coverage <- function(surface, true_mu_delta, true_log_or_delta, conf_level) {
  threshold <- stats::qchisq(conf_level, df = 2)

  if (!(is.list(surface) && all(c("mu_delta", "log_or_delta", "surface") %in% names(surface)))) {
    return(list(
      covered = NA,
      outside_grid = NA,
      nonfinite_endpoint = NA,
      joint_area = NA_real_,
      joint_mu_range = NA_real_,
      joint_log_or_range = NA_real_,
      error_class = "invalid_surface",
      error_message = "Joint confidence result did not contain the expected surface fields."
    ))
  }

  x <- as.numeric(surface$mu_delta)
  y <- as.numeric(surface$log_or_delta)
  z <- as.matrix(surface$surface)

  if (length(x) < 2L || length(y) < 2L || !all(dim(z) == c(length(x), length(y))) ||
      anyDuplicated(x) || anyDuplicated(y) || !all(is.finite(c(x, y)))) {
    return(list(
      covered = NA,
      outside_grid = NA,
      nonfinite_endpoint = TRUE,
      joint_area = NA_real_,
      joint_mu_range = NA_real_,
      joint_log_or_range = NA_real_,
      error_class = "invalid_surface",
      error_message = "Joint confidence surface axes were invalid."
    ))
  }

  x_order <- order(x)
  y_order <- order(y)
  x <- x[x_order]
  y <- y[y_order]
  z <- z[x_order, y_order, drop = FALSE]
  size <- ci_study_joint_size(x, y, z, threshold)

  if (!all(is.finite(c(true_mu_delta, true_log_or_delta)))) {
    return(list(
      covered = NA,
      outside_grid = NA,
      nonfinite_endpoint = NA,
      joint_area = size$area,
      joint_mu_range = size$mu_range,
      joint_log_or_range = size$log_or_range,
      error_class = "invalid_truth",
      error_message = "The joint target truth was not finite."
    ))
  }

  outside <- true_mu_delta < min(x) || true_mu_delta > max(x) ||
    true_log_or_delta < min(y) || true_log_or_delta > max(y)
  if (outside) {
    return(list(
      covered = NA,
      outside_grid = TRUE,
      nonfinite_endpoint = FALSE,
      joint_area = size$area,
      joint_mu_range = size$mu_range,
      joint_log_or_range = size$log_or_range,
      error_class = "outside_grid",
      error_message = "The true joint parameter lay outside the evaluated surface grid."
    ))
  }

  x_index <- ci_study_joint_bracket(x, true_mu_delta)
  y_index <- ci_study_joint_bracket(y, true_log_or_delta)
  corners <- z[x_index + 0:1, y_index + 0:1, drop = FALSE]
  if (!all(is.finite(corners))) {
    return(list(
      covered = NA,
      outside_grid = FALSE,
      nonfinite_endpoint = TRUE,
      joint_area = size$area,
      joint_mu_range = size$mu_range,
      joint_log_or_range = size$log_or_range,
      error_class = "nonfinite_surface",
      error_message = "The interpolation neighborhood contained non-finite surface values."
    ))
  }

  x1 <- x[[x_index]]
  x2 <- x[[x_index + 1L]]
  y1 <- y[[y_index]]
  y2 <- y[[y_index + 1L]]
  wx <- if (x2 == x1) 0 else (true_mu_delta - x1) / (x2 - x1)
  wy <- if (y2 == y1) 0 else (true_log_or_delta - y1) / (y2 - y1)
  interpolated <- (1 - wx) * (1 - wy) * corners[1, 1] +
    wx * (1 - wy) * corners[2, 1] +
    (1 - wx) * wy * corners[1, 2] +
    wx * wy * corners[2, 2]

  list(
    covered = if (is.finite(interpolated)) interpolated <= threshold else NA,
    outside_grid = FALSE,
    nonfinite_endpoint = !is.finite(interpolated),
    joint_area = size$area,
    joint_mu_range = size$mu_range,
    joint_log_or_range = size$log_or_range,
    error_class = if (is.finite(interpolated)) "none" else "nonfinite_surface",
    error_message = if (is.finite(interpolated)) NA_character_ else "The interpolated joint statistic was non-finite."
  )
}

ci_study_joint_record <- function(fit, method, truths, fit_runtime_sec, rep_index, joint_resolution) {
  measure <- ci_study_measure_expr(
    ci_study_silent_confint(
      fit,
      parameter = "joint",
      conf.level = fit$conf.level,
      plot = FALSE,
      resolution = joint_resolution
    )
  )

  if (!is.null(measure$error)) {
    return(ci_study_make_record(
      rep_index = rep_index,
      method = method,
      interval_target = "joint_region",
      fit_success = TRUE,
      fit_runtime_sec = fit_runtime_sec,
      interval_success = FALSE,
      interval_runtime_sec = measure$elapsed,
      error_class = "interval_error",
      error_message = ci_study_condition_message(measure$error)
    ))
  }

  coverage <- ci_study_joint_coverage(
    measure$value,
    truths$true_mu_delta,
    truths$true_log_or_delta,
    fit$conf.level
  )

  ci_study_make_record(
    rep_index = rep_index,
    method = method,
    interval_target = "joint_region",
    fit_success = TRUE,
    fit_runtime_sec = fit_runtime_sec,
    interval_success = !identical(coverage$error_class, "invalid_surface"),
    covered = coverage$covered,
    nonfinite_endpoint = coverage$nonfinite_endpoint,
    outside_grid = coverage$outside_grid,
    width = coverage$joint_area,
    interval_runtime_sec = measure$elapsed,
    joint_area = coverage$joint_area,
    joint_mu_range = coverage$joint_mu_range,
    joint_log_or_range = coverage$joint_log_or_range,
    error_class = coverage$error_class,
    error_message = coverage$error_message
  )
}

ci_study_fit_failure_records <- function(rep_index,
                                         method,
                                         interval_targets,
                                         fit_runtime_sec,
                                         error_class,
                                         error_message) {
  lapply(interval_targets, function(interval_target) {
    ci_study_make_record(
      rep_index = rep_index,
      method = method,
      interval_target = interval_target,
      fit_success = FALSE,
      fit_runtime_sec = fit_runtime_sec,
      interval_success = FALSE,
      error_class = error_class,
      error_message = error_message
    )
  })
}

ci_study_evaluate_fit <- function(data, method, atom, conf_level, truths, rep_index, joint_resolution) {
  interval_targets <- ci_study_interval_targets()$interval_target
  fit_measure <- ci_study_fit_method(data, method, atom, conf_level)
  fit <- fit_measure$value

  fit_success <- is.null(fit_measure$error) &&
    inherits(fit, "trunc_comp_fit") &&
    isTRUE(fit$success)

  if (!fit_success) {
    error_class <- if (!is.null(fit_measure$error)) "fit_error" else "fit_failed"
    error_message <- if (!is.null(fit_measure$error)) {
      ci_study_condition_message(fit_measure$error)
    } else if (!is.null(fit$error)) {
      fit$error
    } else {
      "Fit did not return a successful trunc_comp_fit object."
    }

    return(ci_study_fit_failure_records(
      rep_index = rep_index,
      method = method,
      interval_targets = interval_targets,
      fit_runtime_sec = fit_measure$elapsed,
      error_class = error_class,
      error_message = error_message
    ))
  }

  c(
    ci_study_component_records(fit, method, truths, fit_measure$elapsed, rep_index),
    list(ci_study_joint_record(fit, method, truths, fit_measure$elapsed, rep_index, joint_resolution)),
    list(ci_study_delta_record(
      fit = fit,
      method = method,
      interval_target = "delta_welch",
      delta_method = "welch",
      truth = truths$true_delta,
      fit_runtime_sec = fit_measure$elapsed,
      rep_index = rep_index
    )),
    list(ci_study_delta_record(
      fit = fit,
      method = method,
      interval_target = "delta_projected_optimize",
      delta_method = "projected",
      truth = truths$true_delta,
      fit_runtime_sec = fit_measure$elapsed,
      rep_index = rep_index,
      algorithm = "optimize"
    )),
    list(ci_study_delta_record(
      fit = fit,
      method = method,
      interval_target = "delta_profile_optimize",
      delta_method = "profile",
      truth = truths$true_delta,
      fit_runtime_sec = fit_measure$elapsed,
      rep_index = rep_index,
      algorithm = "optimize"
    ))
  )
}

ci_study_true_count <- function(x) {
  sum(!is.na(x) & x)
}

ci_study_finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

ci_study_finite_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

ci_study_finite_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x)) {
    as.numeric(stats::quantile(x, probs = prob, names = FALSE, type = 7))
  } else {
    NA_real_
  }
}

ci_study_aggregate_interval_observations <- function(observations) {
  if (!nrow(observations)) {
    return(data.frame())
  }

  target_info <- ci_study_interval_targets()
  method_labels <- ci_study_method_labels()
  group_keys <- interaction(
    observations$cell_id,
    observations$method,
    observations$interval_target,
    drop = TRUE,
    lex.order = TRUE
  )
  groups <- split(seq_len(nrow(observations)), group_keys)

  rows <- lapply(groups, function(index) {
    x <- observations[index, , drop = FALSE]
    coverage_valid <- !is.na(x$covered)
    coverage_reps <- sum(coverage_valid)
    coverage_rate <- if (coverage_reps) {
      mean(x$covered[coverage_valid])
    } else {
      NA_real_
    }

    data.frame(
      cell_id = x$cell_id[[1]],
      method = x$method[[1]],
      method_label = unname(method_labels[[x$method[[1]]]]),
      interval_target = x$interval_target[[1]],
      interval_label = target_info$interval_label[match(x$interval_target[[1]], target_info$interval_target)],
      interval_order = target_info$interval_order[match(x$interval_target[[1]], target_info$interval_target)],
      attempted_reps = nrow(x),
      fit_successful_reps = sum(x$fit_success),
      successful_reps = coverage_reps,
      covered_reps = ci_study_true_count(x$covered),
      coverage_rate = coverage_rate,
      coverage_mcse = if (coverage_reps) sqrt(coverage_rate * (1 - coverage_rate) / coverage_reps) else NA_real_,
      fit_failure_count = sum(!x$fit_success),
      fit_failure_rate = mean(!x$fit_success),
      interval_failure_count = sum(x$fit_success & is.na(x$covered)),
      interval_failure_rate = mean(x$fit_success & is.na(x$covered)),
      interval_error_count = sum(
        x$fit_success &
          x$error_class %in% c("interval_error", "invalid_interval", "invalid_endpoint_order", "invalid_surface")
      ),
      nonfinite_endpoint_count = ci_study_true_count(x$nonfinite_endpoint),
      nonfinite_endpoint_rate = mean(!is.na(x$nonfinite_endpoint) & x$nonfinite_endpoint),
      outside_grid_count = ci_study_true_count(x$outside_grid),
      outside_grid_rate = mean(!is.na(x$outside_grid) & x$outside_grid),
      mean_width = ci_study_finite_mean(x$width),
      median_width = ci_study_finite_median(x$width),
      q25_width = ci_study_finite_quantile(x$width, 0.25),
      q75_width = ci_study_finite_quantile(x$width, 0.75),
      mean_runtime_sec = ci_study_finite_mean(x$interval_runtime_sec),
      median_runtime_sec = ci_study_finite_median(x$interval_runtime_sec),
      mean_fit_runtime_sec = ci_study_finite_mean(x$fit_runtime_sec),
      median_fit_runtime_sec = ci_study_finite_median(x$fit_runtime_sec),
      mean_joint_area = ci_study_finite_mean(x$joint_area),
      median_joint_area = ci_study_finite_median(x$joint_area),
      mean_joint_mu_range = ci_study_finite_mean(x$joint_mu_range),
      mean_joint_log_or_range = ci_study_finite_mean(x$joint_log_or_range),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$cell_id, out$method, out$interval_order), , drop = FALSE]
}

ci_study_aggregate_diagnostics <- function(observations) {
  diagnostics <- observations[
    !is.na(observations$error_class) & observations$error_class != "none",
    c("cell_id", "method", "interval_target", "error_class"),
    drop = FALSE
  ]
  if (!nrow(diagnostics)) {
    return(data.frame())
  }

  diagnostics$count <- 1L
  out <- stats::aggregate(
    count ~ cell_id + method + interval_target + error_class,
    data = diagnostics,
    FUN = sum
  )
  out[order(out$cell_id, out$method, out$interval_target, out$error_class), , drop = FALSE]
}

ci_study_run_chunk <- function(chunk, output_dir) {
  chunk <- chunk[1, , drop = FALSE]
  scenario <- simulation_study_scenarios()[[chunk$scenario_id]]
  params <- scenario$effect_parameters(chunk$h)
  methods <- ci_study_methods()
  global_rep_indices <- seq.int(chunk$rep_start, chunk$rep_end)
  records <- vector(
    "list",
    length(global_rep_indices) * length(methods) * nrow(ci_study_interval_targets())
  )
  record_index <- 0L
  combined_deltas <- numeric(length(global_rep_indices))

  truths <- list(
    true_mu_delta = chunk$true_mu_delta,
    true_log_or_delta = chunk$true_log_or_delta,
    true_alpha_delta = chunk$true_alpha_delta,
    true_delta = chunk$true_delta
  )

  for (local_index in seq_along(global_rep_indices)) {
    rep_index <- global_rep_indices[[local_index]]
    set.seed(as.integer(chunk$seed + rep_index - 1L))
    data <- TruncComp2::simulate_truncated_data(
      n = chunk$n,
      f0 = params$f0,
      f1 = params$f1,
      pi0 = params$pi0,
      pi1 = params$pi1,
      atom = chunk$atom
    )

    combined_deltas[[local_index]] <- mean(data$Y[data$R == 1]) - mean(data$Y[data$R == 0])

    for (method in methods) {
      method_records <- ci_study_evaluate_fit(
        data = data,
        method = method,
        atom = chunk$atom,
        conf_level = chunk$conf_level,
        truths = truths,
        rep_index = rep_index,
        joint_resolution = chunk$joint_resolution
      )

      for (record in method_records) {
        record_index <- record_index + 1L
        records[[record_index]] <- record
      }
    }
  }

  interval_observations <- do.call(rbind, records[seq_len(record_index)])
  interval_observations$cell_id <- chunk$cell_id
  interval_observations$chunk_id <- chunk$chunk_id
  interval_observations <- interval_observations[, c(
    "cell_id", "chunk_id", "rep_index", "method", "interval_target",
    "fit_success", "interval_success", "covered", "nonfinite_endpoint",
    "outside_grid", "width", "interval_runtime_sec", "fit_runtime_sec",
    "joint_area", "joint_mu_range", "joint_log_or_range",
    "error_class", "error_message"
  )]

  chunk_summary <- data.frame(
    cell_id = chunk$cell_id,
    scenario_id = chunk$scenario_id,
    scenario_order = chunk$scenario_order,
    scenario_label = chunk$scenario_label,
    short_label = chunk$short_label,
    purpose = chunk$purpose,
    n = chunk$n,
    h = chunk$h,
    is_null = chunk$is_null,
    reps = chunk$reps,
    chunk_reps = chunk$chunk_reps,
    expected_chunks = chunk$expected_chunks,
    chunk_id = chunk$chunk_id,
    rep_start = chunk$rep_start,
    rep_end = chunk$rep_end,
    chunk_reps_actual = chunk$chunk_reps_actual,
    attempted_reps = length(global_rep_indices),
    conf_level = chunk$conf_level,
    atom = chunk$atom,
    base_seed = chunk$base_seed,
    seed = chunk$seed,
    joint_resolution = chunk$joint_resolution,
    effect_text = chunk$effect_text,
    pi0 = chunk$pi0,
    pi1 = chunk$pi1,
    survivor_mean0 = chunk$survivor_mean0,
    survivor_mean1 = chunk$survivor_mean1,
    true_mu_delta = chunk$true_mu_delta,
    true_log_or_delta = chunk$true_log_or_delta,
    true_alpha_delta = chunk$true_alpha_delta,
    true_delta = chunk$true_delta,
    empirical_combined_delta_sum = sum(combined_deltas),
    empirical_combined_delta_sumsq = sum(combined_deltas ^ 2),
    empirical_combined_delta = mean(combined_deltas),
    empirical_combined_delta_sd = if (length(combined_deltas) > 1L) stats::sd(combined_deltas) else NA_real_,
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE
  )

  result <- list(
    version = "trunccomp2-ci-study-cell-v1",
    chunk_summary = chunk_summary,
    interval_metrics = ci_study_aggregate_interval_observations(interval_observations),
    replicate_diagnostics = ci_study_aggregate_diagnostics(interval_observations),
    interval_observations = interval_observations
  )

  ci_study_write_rds_atomic(result, ci_study_chunk_path(output_dir, chunk))
  result
}

ci_study_infer_config_from_chunks <- function(output_dir) {
  chunk_paths <- ci_study_existing_chunk_paths(output_dir)
  if (!length(chunk_paths)) {
    stop("No CI-study chunk files were found. Run the CI-study driver first.", call. = FALSE)
  }

  chunk_results <- lapply(chunk_paths, readRDS)
  chunk_summaries <- do.call(rbind, lapply(chunk_results, function(x) x$chunk_summary))
  rownames(chunk_summaries) <- NULL

  singleton <- function(field) {
    values <- unique(chunk_summaries[[field]])
    if (length(values) != 1L) {
      stop("Unable to infer a single CI-study configuration from field ", field, ".", call. = FALSE)
    }
    values[[1]]
  }

  ci_study_default_config(
    reps = as.integer(singleton("reps")),
    chunk_reps = as.integer(singleton("chunk_reps")),
    scenario_ids = unique(chunk_summaries$scenario_id[order(chunk_summaries$scenario_order)]),
    n_seq = sort(unique(as.integer(chunk_summaries$n))),
    effect_levels = sort(unique(as.integer(chunk_summaries$h))),
    conf_level = as.numeric(singleton("conf_level")),
    atom = as.numeric(singleton("atom")),
    base_seed = as.integer(singleton("base_seed")),
    joint_resolution = as.integer(singleton("joint_resolution"))
  )
}

ci_study_aggregate_cell_summaries <- function(chunk_summaries, design) {
  expected_chunks <- data.frame(
    cell_id = design$cell_id,
    expected_chunks = design$expected_chunks,
    stringsAsFactors = FALSE
  )

  rows <- lapply(seq_len(nrow(design)), function(index) {
    cell <- design[index, , drop = FALSE]
    chunks <- chunk_summaries[chunk_summaries$cell_id == cell$cell_id, , drop = FALSE]
    completed_chunks <- nrow(chunks)
    attempted_reps <- if (completed_chunks) sum(chunks$attempted_reps) else 0L
    empirical_sum <- if (completed_chunks) sum(chunks$empirical_combined_delta_sum) else NA_real_
    empirical_sumsq <- if (completed_chunks) sum(chunks$empirical_combined_delta_sumsq) else NA_real_
    empirical_mean <- if (attempted_reps > 0L) empirical_sum / attempted_reps else NA_real_
    empirical_sd <- if (attempted_reps > 1L) {
      sqrt(max(0, (empirical_sumsq - empirical_sum ^ 2 / attempted_reps) / (attempted_reps - 1L)))
    } else {
      NA_real_
    }

    data.frame(
      cell,
      completed_chunks = completed_chunks,
      cell_complete = completed_chunks == expected_chunks$expected_chunks[expected_chunks$cell_id == cell$cell_id],
      attempted_reps = attempted_reps,
      empirical_combined_delta = empirical_mean,
      empirical_combined_delta_sd = empirical_sd,
      completed_at = if (completed_chunks) max(chunks$completed_at) else NA_character_,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$scenario_order, out$h, out$n), , drop = FALSE]
}

aggregate_ci_study_results <- function(output_dir, config = NULL) {
  if (is.null(config)) {
    config_path <- ci_study_config_path(output_dir)
    if (file.exists(config_path)) {
      config <- readRDS(config_path)
    } else {
      config <- ci_study_infer_config_from_chunks(output_dir)
      saveRDS(config, file = config_path)
    }
  }

  design <- ci_study_design(config)
  manifest <- ci_study_chunk_manifest(design, config)
  chunk_paths <- vapply(
    split(manifest, seq_len(nrow(manifest))),
    function(chunk) ci_study_chunk_path(output_dir, chunk),
    character(1)
  )
  existing <- file.exists(chunk_paths)

  if (!any(existing)) {
    stop(
      "No CI-study chunk files were found. Copy the raw chunk outputs into simulation-study/results/.../cells/ or run the CI-study driver first.",
      call. = FALSE
    )
  }

  chunk_results <- lapply(chunk_paths[existing], readRDS)
  chunk_summaries <- do.call(rbind, lapply(chunk_results, function(x) x$chunk_summary))
  interval_observations <- do.call(rbind, lapply(chunk_results, function(x) x$interval_observations))
  rownames(chunk_summaries) <- NULL
  rownames(interval_observations) <- NULL

  cell_metrics <- ci_study_aggregate_cell_summaries(chunk_summaries, design)
  interval_metrics <- ci_study_aggregate_interval_observations(interval_observations)
  interval_metrics <- merge(
    cell_metrics[, c(
      "cell_id", "scenario_id", "scenario_order", "scenario_label", "short_label",
      "purpose", "n", "h", "is_null", "reps", "chunk_reps", "conf_level",
      "atom", "seed", "joint_resolution", "effect_text", "pi0", "pi1",
      "survivor_mean0", "survivor_mean1", "true_mu_delta", "true_log_or_delta",
      "true_alpha_delta", "true_delta", "empirical_combined_delta",
      "empirical_combined_delta_sd", "cell_complete"
    )],
    interval_metrics,
    by = "cell_id",
    sort = FALSE
  )
  interval_metrics <- interval_metrics[
    order(
      interval_metrics$scenario_order,
      interval_metrics$h,
      interval_metrics$n,
      interval_metrics$method,
      interval_metrics$interval_order
    ),
    ,
    drop = FALSE
  ]

  diagnostics <- ci_study_aggregate_diagnostics(interval_observations)
  if (nrow(diagnostics)) {
    diagnostics <- merge(
      cell_metrics[, c("cell_id", "scenario_id", "scenario_order", "n", "h", "is_null")],
      diagnostics,
      by = "cell_id",
      sort = FALSE
    )
  }

  results <- list(
    version = config$version,
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    config = config,
    complete = all(existing),
    expected_cells = nrow(design),
    completed_cells = sum(cell_metrics$cell_complete),
    expected_chunks = nrow(manifest),
    completed_chunks = sum(existing),
    scenarios = ci_study_scenario_table(config),
    design = design,
    manifest = manifest,
    cell_metrics = cell_metrics,
    interval_metrics = interval_metrics,
    diagnostics = diagnostics
  )

  saveRDS(results, file = file.path(output_dir, "ci-study.rds"))
  utils::write.csv(
    results$cell_metrics,
    file = file.path(output_dir, "ci-study-cell-metrics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    results$interval_metrics,
    file = file.path(output_dir, "ci-study-interval-metrics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    results$diagnostics,
    file = file.path(output_dir, "ci-study-diagnostics.csv"),
    row.names = FALSE
  )

  results
}

run_ci_study <- function(repo_root,
                         output_dir,
                         reps = 10000L,
                         chunk_reps = 500L,
                         workers = 1L,
                         scenario_ids = NULL,
                         n_seq = NULL,
                         effect_levels = NULL,
                         conf_level = 0.95,
                         joint_resolution = 35L,
                         overwrite = FALSE) {
  prepared <- ci_study_prepare_run(
    output_dir = output_dir,
    reps = reps,
    chunk_reps = chunk_reps,
    scenario_ids = scenario_ids,
    n_seq = n_seq,
    effect_levels = effect_levels,
    conf_level = conf_level,
    joint_resolution = joint_resolution,
    overwrite = overwrite
  )

  load_local_trunccomp2(repo_root)

  if (nrow(prepared$pending_manifest)) {
    chunk_tasks <- split(prepared$pending_manifest, seq_len(nrow(prepared$pending_manifest)))
    worker_fun <- function(chunk) {
      message(sprintf(
        "Running %s, h = %d, n = %d, chunk = %d/%d",
        chunk$scenario_id,
        chunk$h,
        chunk$n,
        chunk$chunk_id,
        chunk$expected_chunks
      ))
      ci_study_run_chunk(chunk, output_dir = prepared$output_dir)
    }

    if (.Platform$OS.type != "windows" && workers > 1L) {
      parallel::mclapply(chunk_tasks, worker_fun, mc.cores = workers, mc.preschedule = FALSE)
    } else {
      lapply(chunk_tasks, worker_fun)
    }
  }

  aggregate_ci_study_results(prepared$output_dir, config = prepared$config)
}

ci_study_finalize_results <- function(output_dir, ci_results) {
  message(
    sprintf(
      "CI study aggregation complete: %d/%d chunks and %d/%d cells available. Results written to %s",
      ci_results$completed_chunks,
      ci_results$expected_chunks,
      ci_results$completed_cells,
      ci_results$expected_cells,
      file.path(output_dir, "ci-study.rds")
    )
  )

  invisible(ci_results)
}

ci_study_submit_slurm_jobs <- function(repo_root,
                                       output_dir,
                                       manifest_path,
                                       chunk_count,
                                       partition = "standard",
                                       time = "04:00:00",
                                       cpus_per_task = 1L,
                                       mem = "4G",
                                       array_parallelism = 80L,
                                       collector_time = "00:30:00",
                                       collector_mem = "8G",
                                       job_name_prefix = "trunccomp2_ci",
                                       dry_run = FALSE) {
  if (!(length(cpus_per_task) == 1L && is.numeric(cpus_per_task) && is.finite(cpus_per_task) &&
        cpus_per_task >= 1)) {
    stop("cpus_per_task must be a single positive integer.", call. = FALSE)
  }

  study_dir <- ci_study_root(repo_root)
  chunk_script <- file.path(study_dir, "slurm", "run-ci-study-cell.sbatch")
  collect_script <- file.path(study_dir, "slurm", "collect-ci-study-results.sbatch")
  if (!file.exists(chunk_script)) {
    stop("Missing Slurm CI-study chunk script at ", chunk_script, call. = FALSE)
  }
  if (!file.exists(collect_script)) {
    stop("Missing Slurm CI-study collector script at ", collect_script, call. = FALSE)
  }

  slurm_dir <- ci_study_slurm_dir(output_dir)
  array_spec <- .simulation_study_array_spec(chunk_count, array_parallelism = array_parallelism)

  chunk_args <- c(
    "--parsable",
    paste0("--job-name=", job_name_prefix, "_chunks"),
    paste0("--partition=", partition),
    paste0("--time=", time),
    paste0("--cpus-per-task=", as.integer(cpus_per_task)),
    paste0("--array=", array_spec),
    paste0("--output=", file.path(slurm_dir, "chunk-%A_%a.out")),
    paste0("--error=", file.path(slurm_dir, "chunk-%A_%a.err")),
    paste0("--chdir=", repo_root)
  )
  if (!is.null(mem) && nzchar(mem)) {
    chunk_args <- c(chunk_args, paste0("--mem=", mem))
  }
  chunk_args <- c(chunk_args, chunk_script, repo_root, output_dir, manifest_path)

  chunk_submission <- .simulation_study_submit_sbatch(chunk_args, dry_run = dry_run)
  dependency_job_id <- if (dry_run) "<array-job-id>" else chunk_submission$job_id

  collect_args <- c(
    "--parsable",
    paste0("--job-name=", job_name_prefix, "_collect"),
    paste0("--partition=", partition),
    paste0("--time=", collector_time),
    "--cpus-per-task=1",
    paste0("--dependency=afterany:", dependency_job_id),
    paste0("--output=", file.path(slurm_dir, "collect-%j.out")),
    paste0("--error=", file.path(slurm_dir, "collect-%j.err")),
    paste0("--chdir=", repo_root)
  )
  if (!is.null(collector_mem) && nzchar(collector_mem)) {
    collect_args <- c(collect_args, paste0("--mem=", collector_mem))
  }
  collect_args <- c(collect_args, collect_script, repo_root, output_dir)

  collect_submission <- .simulation_study_submit_sbatch(collect_args, dry_run = dry_run)

  metadata <- c(
    sprintf("generated_at=%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("manifest=%s", manifest_path),
    sprintf("chunk_count=%d", as.integer(chunk_count)),
    sprintf("array_job_id=%s", if (dry_run) "<dry-run>" else chunk_submission$job_id),
    sprintf("collector_job_id=%s", if (dry_run) "<dry-run>" else collect_submission$job_id),
    sprintf("array_command=%s", chunk_submission$command),
    sprintf("collector_command=%s", collect_submission$command)
  )
  writeLines(metadata, con = ci_study_slurm_submission_path(output_dir), useBytes = TRUE)

  list(
    manifest_path = manifest_path,
    chunk_count = as.integer(chunk_count),
    array_job_id = chunk_submission$job_id,
    collector_job_id = collect_submission$job_id,
    array_command = chunk_submission$command,
    collector_command = collect_submission$command,
    dry_run = dry_run
  )
}
