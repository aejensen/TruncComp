args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) default else sub(prefix, "", hit[[1]])
}
arg_numeric_csv <- function(name, default) {
  value <- arg_value(name, default)
  out <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1]]))
  if (!length(out) || any(!is.finite(out)) || any(out <= 0)) {
    stop("--", name, " must be a comma-separated list of positive finite numbers.", call. = FALSE)
  }
  sort(unique(out))
}
parse_numeric_pair <- function(value, arg_name) {
  out <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1]]))
  if (length(out) != 2L || any(!is.finite(out)) || any(out <= 0)) {
    stop("--", arg_name, " must be two comma-separated positive finite numbers.", call. = FALSE)
  }
  out
}
arg_numeric_pair <- function(name, default) {
  parse_numeric_pair(arg_value(name, default), name)
}
format_tag_number <- function(x) {
  gsub("\\.", "p", format(x, trim = TRUE, scientific = FALSE))
}
rbind_fill <- function(rows) {
  if (!length(rows)) return(data.frame())
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(all_names, names(row))
    for (name in missing) {
      row[[name]] <- NA
    }
    row[all_names]
  })
  do.call(rbind, rows)
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) {
  stop("Unable to determine script path.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
source(file.path(script_dir, "ist3_bayes_speed_utils.R"), local = globalenv())

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("cmdstanr is required for this experiment.", call. = FALSE)
}
if (!requireNamespace("posterior", quietly = TRUE)) {
  stop("posterior is required for this experiment.", call. = FALSE)
}

mode <- arg_value("mode", "all")
if (!mode %in% c("validate", "screen", "fit", "all")) {
  stop("--mode must be validate, screen, fit, or all.", call. = FALSE)
}

seed <- as.integer(arg_value("seed", "20260419"))
target_final_seconds <- as.numeric(arg_value("target-final-seconds", "1800"))
fit_adapt_delta_override <- arg_value("fit-adapt-delta", NULL)
fit_adapt_delta_override <- if (is.null(fit_adapt_delta_override)) {
  NA_real_
} else {
  as.numeric(fit_adapt_delta_override)
}
if (!is.na(fit_adapt_delta_override) && (!is.finite(fit_adapt_delta_override) || fit_adapt_delta_override <= 0 || fit_adapt_delta_override >= 1)) {
  stop("--fit-adapt-delta must be a number strictly between 0 and 1.", call. = FALSE)
}
stop_after_medium <- tolower(arg_value("stop-after-medium", "false")) %in% c("1", "true", "yes")
threads_per_chain <- as.integer(arg_value("threads-per-chain", "1"))
grainsize <- as.integer(arg_value("grainsize", if (threads_per_chain > 1L) "16" else "128"))
screen_warmup <- as.integer(arg_value("screen-warmup", "50"))
screen_sampling <- as.integer(arg_value("screen-sampling", "50"))
medium_warmup <- as.integer(arg_value("medium-warmup", "150"))
medium_sampling <- as.integer(arg_value("medium-sampling", "150"))
final_warmup <- as.integer(arg_value("final-warmup", "300"))
final_sampling <- as.integer(arg_value("final-sampling", "300"))
model_variant <- arg_value("model", "h2_ordered")
if (!model_variant %in% c("h2_ordered", "h2_ordered_logit", "h2_ordered_logit_ncstick", "h2_ordered_logit_marginal_stick", "h1_single", "h1_logitnormal", "h2_logitnormal", "h2_logitnormal_bounded_mu", "h2_logitnormal_bounded_mu_collapsed_rho")) {
  stop("--model must be h2_ordered, h2_ordered_logit, h2_ordered_logit_ncstick, h2_ordered_logit_marginal_stick, h1_single, h1_logitnormal, h2_logitnormal, h2_logitnormal_bounded_mu, or h2_logitnormal_bounded_mu_collapsed_rho.", call. = FALSE)
}
model_H <- if (model_variant %in% c("h1_single", "h1_logitnormal")) 1L else 2L
uses_stick_quadrature <- model_variant %in% c(
  "h2_ordered_logit_marginal_stick"
)
heaping_grids <- arg_numeric_csv("heaping-grids", "1,5,10")
stick_quad_n <- as.integer(arg_value("stick-quad-n", "64"))
init_strategy <- arg_value("init-strategy", "empirical")
if (!init_strategy %in% c("empirical", "jittered", "control_low_main", "control_high_main")) {
  stop("--init-strategy must be empirical, jittered, control_low_main, or control_high_main.", call. = FALSE)
}
if (identical(init_strategy, "jittered")) {
  init_strategy <- "empirical"
}
gap_prior <- arg_numeric_pair("gap-prior", "1,1")
fixed_v_prior_arg <- arg_value("fixed-v-prior", NULL)
use_fixed_v_prior <- !is.null(fixed_v_prior_arg)
fixed_v_prior <- if (use_fixed_v_prior) {
  parse_numeric_pair(fixed_v_prior_arg, "fixed-v-prior")
} else {
  c(2, 2)
}
if (use_fixed_v_prior) {
  stop("--fixed-v-prior is disabled for these experiments; stick weights stay beta(1, alpha).", call. = FALSE)
}
candidate_labels <- arg_value("candidates", NULL)
heaping_label <- if (identical(heaping_grids, 1)) {
  "no_heaping"
} else {
  paste0("heap", paste(format(heaping_grids, trim = TRUE, scientific = FALSE), collapse = "_"))
}
prior_control_label <- paste0(
  "gap", format_tag_number(gap_prior[[1]]), "_", format_tag_number(gap_prior[[2]]),
  if (use_fixed_v_prior) {
    paste0("-fixedv", format_tag_number(fixed_v_prior[[1]]), "_", format_tag_number(fixed_v_prior[[2]]))
  } else {
    "-stickv"
  }
)
artifact_tag <- paste0("-", model_variant, "-", heaping_label, "-", init_strategy, "-", prior_control_label)
build_dir <- ensure_dir(file.path(repo_root, "manuscript", "build", "ist3-bayes-speed"))
draws_dir <- ensure_dir(file.path(build_dir, "draws"))
model_dir <- ensure_dir(file.path(build_dir, "cmdstan"))

message("IST-3 bounded-score speed experiment")
message("Repository: ", repo_root)
message("Build directory: ", build_dir)
message("Model variant: ", model_variant, "; H = ", model_H)
message("Heaping grids: ", paste(heaping_grids, collapse = ", "))
message("Init strategy: ", init_strategy)
message("Gap prior: beta(", gap_prior[[1]], ", ", gap_prior[[2]], ")")
if (use_fixed_v_prior) {
  message("Stick weight prior: fixed beta(", fixed_v_prior[[1]], ", ", fixed_v_prior[[2]], ")")
} else {
  message("Stick weight prior: beta(1, alpha)")
}
message("Threads per chain: ", threads_per_chain, "; reduce_sum grainsize: ", grainsize)
if (uses_stick_quadrature) {
  message("Stick quadrature nodes: ", stick_quad_n)
}
if (is.finite(fit_adapt_delta_override)) {
  message("Fit adapt_delta override: ", fit_adapt_delta_override)
}
if (stop_after_medium) {
  message("Will stop after medium validation.")
}
message(
  "Iterations: screen ", screen_warmup, "/", screen_sampling,
  "; medium ", medium_warmup, "/", medium_sampling,
  "; final ", final_warmup, "/", final_sampling
)

d <- load_ist3_experiment_data(repo_root)
candidates <- ist3_prior_candidates()
candidates <- lapply(candidates, function(prior) {
  prior$gap_prior_alpha <- gap_prior[[1]]
  prior$gap_prior_beta <- gap_prior[[2]]
  prior$use_fixed_v_prior <- use_fixed_v_prior
  prior$v_prior_alpha <- fixed_v_prior[[1]]
  prior$v_prior_beta <- fixed_v_prior[[2]]
  prior$kernel <- if (model_variant %in% c("h1_logitnormal", "h2_logitnormal", "h2_logitnormal_bounded_mu", "h2_logitnormal_bounded_mu_collapsed_rho")) "logitnormal" else "beta"
  prior$mu_parameterization <- if (model_variant %in% c("h2_logitnormal_bounded_mu", "h2_logitnormal_bounded_mu_collapsed_rho")) "bounded_ordered" else "ordered"
  prior$rho_parameterization <- if (identical(model_variant, "h2_logitnormal_bounded_mu_collapsed_rho")) "collapsed" else "sampled"
  prior$m_parameterization <- if (model_variant %in% c("h2_ordered_logit", "h2_ordered_logit_ncstick", "h2_ordered_logit_marginal_stick")) {
    "ordered_logit"
  } else {
    "ordered_gap"
  }
  prior$stick_parameterization <- if (uses_stick_quadrature) {
    "marginal"
  } else if (identical(model_variant, "h2_ordered_logit_ncstick")) {
    "noncentered"
  } else {
    "centered"
  }
  prior
})
if (!is.null(candidate_labels)) {
  requested_candidates <- trimws(strsplit(candidate_labels, ",", fixed = TRUE)[[1]])
  requested_candidates <- requested_candidates[nzchar(requested_candidates)]
  unknown_candidates <- setdiff(requested_candidates, names(candidates))
  if (length(unknown_candidates)) {
    stop(
      "Unknown --candidates value(s): ", paste(unknown_candidates, collapse = ", "),
      ". Available: ", paste(names(candidates), collapse = ", "),
      call. = FALSE
    )
  }
  candidates <- candidates[requested_candidates]
}
message("Prior candidates: ", paste(names(candidates), collapse = ", "))
reference_standata <- build_ist3_standata(
  d,
  candidates[[1]],
  H = model_H,
  heaping_grids = heaping_grids,
  grainsize = grainsize
)
if (uses_stick_quadrature) {
  reference_standata <- add_stick_quadrature(reference_standata, stick_quad_n)
}
validation <- validate_aggregation(reference_standata)
validation_df <- as.data.frame(validation, stringsAsFactors = FALSE)
validation_df$work_reduction <- validation_df$old_beta_interval_work / validation_df$new_beta_interval_work
utils::write.csv(validation_df, file.path(build_dir, paste0("aggregation-validation", artifact_tag, ".csv")), row.names = FALSE)
message(
  "Aggregation: ", validation$n_observed, " survivor rows -> ",
  validation$n_score_cells, " score cells; work reduction ",
  round(validation_df$work_reduction, 2), "x."
)
if (abs(validation$continuous_loglik_diff) > 1e-8) {
  stop("Aggregated and observation-level survivor log likelihoods did not match.", call. = FALSE)
}

prior_rows <- do.call(rbind, lapply(candidates, function(prior) {
  prior_predictive_summary(reference_standata, prior)
}))
utils::write.csv(prior_rows, file.path(build_dir, paste0("prior-predictive-screen", artifact_tag, ".csv")), row.names = FALSE)

if (identical(mode, "validate")) {
  message("Validation mode complete.")
  quit(status = 0)
}

stan_file <- file.path(
  script_dir,
  if (identical(model_variant, "h1_logitnormal")) {
    "ist3_bounded_score_single_logitnormal.stan"
  } else if (identical(model_variant, "h1_single")) {
    "ist3_bounded_score_single_beta.stan"
  } else if (identical(model_variant, "h2_logitnormal")) {
    "ist3_bounded_score_logitnormal_aggregated.stan"
  } else if (identical(model_variant, "h2_logitnormal_bounded_mu")) {
    "ist3_bounded_score_logitnormal_bounded_mu.stan"
  } else if (identical(model_variant, "h2_logitnormal_bounded_mu_collapsed_rho")) {
    "ist3_bounded_score_logitnormal_bounded_mu_collapsed_rho.stan"
  } else if (identical(model_variant, "h2_ordered_logit_marginal_stick")) {
    "ist3_bounded_score_ordered_logit_marginal_stick.stan"
  } else if (identical(model_variant, "h2_ordered_logit_ncstick")) {
    "ist3_bounded_score_ordered_logit_ncstick.stan"
  } else if (identical(model_variant, "h2_ordered_logit")) {
    "ist3_bounded_score_ordered_logit.stan"
  } else {
    "ist3_bounded_score_aggregated.stan"
  }
)
message("Compiling external Stan model with cmdstanr and STAN_THREADS.")
model <- cmdstanr::cmdstan_model(
  stan_file,
  dir = model_dir,
  cpp_options = list(stan_threads = TRUE),
  force_recompile = FALSE
)

benchmark_path <- file.path(build_dir, paste0("benchmark-results", artifact_tag, ".csv"))
if (file.exists(benchmark_path)) {
  existing_benchmark <- utils::read.csv(benchmark_path, stringsAsFactors = FALSE)
  if (!"model_variant" %in% names(existing_benchmark)) {
    existing_benchmark$model_variant <- model_variant
  }
  benchmark_rows <- split(existing_benchmark, seq_len(nrow(existing_benchmark)))
} else {
  benchmark_rows <- list()
}
result_objects <- list()

append_benchmark <- function(result, status, projected_final_sec, artifact_prefix) {
  row <- result$metrics
  row$status <- status
  row$projected_final_sec <- projected_final_sec
  row$artifact_prefix <- artifact_prefix
  benchmark_rows[[length(benchmark_rows) + 1L]] <<- row
  utils::write.csv(rbind_fill(benchmark_rows), benchmark_path, row.names = FALSE)
}

run_and_record <- function(prior, stage, adapt_delta, iter_warmup, iter_sampling, chains, seed_offset) {
  label <- prior$label
  artifact_prefix <- paste(
    stage,
    paste0("w", iter_warmup, "s", iter_sampling),
    model_variant,
    heaping_label,
    init_strategy,
    prior_control_label,
    label,
    paste0("delta", gsub("\\.", "", as.character(adapt_delta))),
    sep = "-"
  )
  standata <- build_ist3_standata(
    d,
    prior,
    H = model_H,
    heaping_grids = heaping_grids,
    grainsize = grainsize
  )
  if (uses_stick_quadrature) {
    standata <- add_stick_quadrature(standata, stick_quad_n)
  }
  result <- run_cmdstan_fit(
    model = model,
    standata = standata,
    prior = prior,
    label = label,
    stage = stage,
    adapt_delta = adapt_delta,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    chains = chains,
    seed = seed + seed_offset,
    output_dir = ensure_dir(file.path(draws_dir, artifact_prefix)),
    threads_per_chain = threads_per_chain,
    init_strategy = init_strategy
  )
  result$metrics$model_variant <- model_variant
  result$metrics$heaping_grids <- paste(heaping_grids, collapse = ",")
  result$metrics$init_strategy <- init_strategy
  result$metrics$gap_prior_alpha <- gap_prior[[1]]
  result$metrics$gap_prior_beta <- gap_prior[[2]]
  result$metrics$use_fixed_v_prior <- use_fixed_v_prior
  result$metrics$v_prior_alpha <- fixed_v_prior[[1]]
  result$metrics$v_prior_beta <- fixed_v_prior[[2]]
  result$metrics$stick_quad_n <- if (uses_stick_quadrature) stick_quad_n else NA_integer_
  status <- diagnostic_acceptance(result$metrics, final = identical(stage, "final"))
  projected_final_sec <- result$metrics$elapsed_sec /
    (chains * (iter_warmup + iter_sampling)) * (4 * (final_warmup + final_sampling))
  append_benchmark(result, status, projected_final_sec, artifact_prefix)
  saveRDS(
    list(
      metrics = result$metrics,
      status = status,
      projected_final_sec = projected_final_sec,
      estimand_summary = result$estimand_summary,
      prior = prior,
      started = result$started,
      finished = result$finished,
      output_files = result$fit$output_files(),
      component_diagnostics = result$component_diagnostics
    ),
    file.path(build_dir, paste0(artifact_prefix, ".rds"))
  )
  result$status <- status
  result$projected_final_sec <- projected_final_sec
  result$artifact_prefix <- artifact_prefix
  result
}

if (mode %in% c("screen", "all")) {
  screen_results <- list()
  adapt_ladder <- c(0.90, 0.95, 0.99)
  candidate_index <- 0L
  for (prior in candidates) {
    candidate_index <- candidate_index + 1L
    for (delta in adapt_ladder) {
      message("Screening ", prior$label, " at adapt_delta=", delta)
      result <- run_and_record(
        prior = prior,
        stage = "screen",
        adapt_delta = delta,
        iter_warmup = screen_warmup,
        iter_sampling = screen_sampling,
        chains = 1,
        seed_offset = 1000L + candidate_index * 10L + which(adapt_ladder == delta)
      )
      screen_results[[length(screen_results) + 1L]] <- result
      if (result$metrics$divergences == 0L) break
    }
  }

  screen_table <- do.call(rbind, lapply(screen_results, function(x) {
    cbind(x$metrics, status = x$status, projected_final_sec = x$projected_final_sec)
  }))
  screen_table$screen_rank <- with(
    screen_table,
    rank(
      divergences * 1e6 +
        (max_treedepth >= 12) * 1e5 +
        ifelse(is.na(min_bfmi), 1e4, pmax(0, 0.3 - min_bfmi) * 1e4) +
        projected_final_sec,
      ties.method = "first"
    )
  )
  utils::write.csv(screen_table, file.path(build_dir, paste0("screen-ranked", artifact_tag, ".csv")), row.names = FALSE)
} else {
  screen_ranked_path <- file.path(build_dir, paste0("screen-ranked", artifact_tag, ".csv"))
  if (!file.exists(screen_ranked_path)) {
    stop("No screen-ranked.csv found. Run --mode=screen before --mode=fit.", call. = FALSE)
  }
  screen_table <- utils::read.csv(screen_ranked_path, stringsAsFactors = FALSE)
}

if (identical(mode, "screen")) {
  message("Screen mode complete.")
  quit(status = 0)
}

screen_table <- screen_table[order(screen_table$screen_rank), , drop = FALSE]
candidate_lookup <- candidates
medium_candidates <- unique(screen_table$label)
selected <- NULL
medium_results <- list()

for (label in medium_candidates) {
  best_screen <- screen_table[screen_table$label == label, , drop = FALSE][1, ]
  prior <- candidate_lookup[[label]]
  delta <- if (is.finite(fit_adapt_delta_override)) fit_adapt_delta_override else best_screen$adapt_delta
  message("Medium validation for ", label, " at adapt_delta=", delta)
  result <- run_and_record(
    prior = prior,
    stage = "medium",
    adapt_delta = delta,
    iter_warmup = medium_warmup,
    iter_sampling = medium_sampling,
    chains = 4,
    seed_offset = 3000L + length(medium_results)
  )
  medium_results[[length(medium_results) + 1L]] <- result
  medium_status <- diagnostic_acceptance(result$metrics, final = FALSE)
  rhat_ok <- is.finite(result$metrics$max_rhat) && result$metrics$max_rhat <= 1.05
  ess_ok <- is.finite(result$metrics$min_bulk_ess) && result$metrics$min_bulk_ess >= 50
  runtime_ok <- result$projected_final_sec <= target_final_seconds
  if (identical(medium_status, "accepted") && rhat_ok && ess_ok && runtime_ok) {
    selected <- list(prior = prior, adapt_delta = delta, medium = result)
    break
  }
}

if (is.null(selected)) {
  message("No medium candidate satisfied all criteria; selecting the best available medium result.")
  medium_table <- do.call(rbind, lapply(medium_results, function(x) {
    cbind(x$metrics, status = x$status, projected_final_sec = x$projected_final_sec)
  }))
  medium_table$rank_score <- with(
    medium_table,
    divergences * 1e7 +
      (max_treedepth >= 12) * 1e6 +
      ifelse(is.na(max_rhat), 1e5, pmax(0, max_rhat - 1.05) * 1e5) +
      ifelse(is.na(min_bulk_ess), 1e4, pmax(0, 50 - min_bulk_ess) * 1e3) +
      projected_final_sec
  )
  medium_table <- medium_table[order(medium_table$rank_score), , drop = FALSE]
  selected_label <- medium_table$label[[1]]
  selected <- list(
    prior = candidate_lookup[[selected_label]],
    adapt_delta = medium_table$adapt_delta[[1]],
    medium = medium_results[[which(vapply(medium_results, function(x) x$metrics$label[[1]] == selected_label, logical(1)))[[1]]]]
  )
}

if (stop_after_medium) {
  message("Stop-after-medium requested; skipping final run.")
  quit(status = 0)
}

message("Final run for ", selected$prior$label, " at adapt_delta=", selected$adapt_delta)
final_result <- run_and_record(
  prior = selected$prior,
  stage = "final",
  adapt_delta = selected$adapt_delta,
  iter_warmup = final_warmup,
  iter_sampling = final_sampling,
  chains = 4,
  seed_offset = 5000L
)
final_status <- diagnostic_acceptance(final_result$metrics, final = TRUE)

summary_path <- file.path(
  build_dir,
  paste0(
    "selected-fit-diagnostics-w",
    final_warmup,
    "s",
    final_sampling,
    artifact_tag,
    ".txt"
  )
)
sink(summary_path)
cat("IST-3 bounded-score speed experiment\n")
cat("Selected prior:", selected$prior$label, "\n")
cat("Adapt delta:", selected$adapt_delta, "\n")
cat("Final status:", final_status, "\n\n")
print(final_result$metrics)
cat("\nEstimand summary\n")
print(final_result$estimand_summary)
cat("\nComponent diagnostics by chain\n")
print(final_result$component_diagnostics)
cat("\nPrior\n")
print(selected$prior)
sink()

if (!identical(final_status, "accepted")) {
  warning("Final run did not satisfy all target criteria: ", final_status, call. = FALSE)
}
message("Final status: ", final_status)
message("Diagnostics written to ", summary_path)
