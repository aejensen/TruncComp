command <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command, value = TRUE)
if (length(file_argument) != 1L) {
  stop("Run this test with Rscript.")
}

test_file <- normalizePath(sub("^--file=", "", file_argument))
manuscript_dir <- dirname(dirname(test_file))
repo_root <- dirname(manuscript_dir)

source(file.path(manuscript_dir, "R", "utils.R"))
source(file.path(manuscript_dir, "R", "data.R"))
source(file.path(manuscript_dir, "R", "analysis.R"))
source(file.path(manuscript_dir, "R", "build.R"))

cache_path <- application_bayes_cache_path(manuscript_dir)
cache_sha256_before <- .application_sha256_file(cache_path)
cached <- readRDS(cache_path)
data <- utils::read.csv(
  file.path(manuscript_dir, "application-data", "ist3-standardized.csv"),
  stringsAsFactors = FALSE
)
contract <- application_bayes_contract(manuscript_dir)

validation <- application_bayes_validate_cache(
  cached = cached,
  d = data,
  cache_path = cache_path,
  manuscript_dir = manuscript_dir
)
stopifnot(validation$valid)
stopifnot(length(validation$errors) == 0L)
stopifnot(identical(cache_sha256_before, contract$locked_cache_sha256))
stopifnot(identical(
  .application_sha256_object(application_bayes_settings()),
  contract$settings_sha256
))

application_source <- paste(
  readLines(application_bayes_model_path(manuscript_dir), warn = FALSE),
  collapse = "\n"
)
generic_score_source <- paste(
  readLines(
    file.path(
      repo_root,
      "inst",
      "stan",
      "trunc_comp_bayes_bounded_score_logit_normal.stan"
    ),
    warn = FALSE
  ),
  collapse = "\n"
)
generic_continuous_source <- paste(
  readLines(
    file.path(
      repo_root,
      "inst",
      "stan",
      "trunc_comp_bayes_bounded_continuous_logit_normal.stan"
    ),
    warn = FALSE
  ),
  collapse = "\n"
)
generic_score_header <- paste(
  readLines(
    file.path(
      repo_root,
      "src",
      "stanExports_trunc_comp_bayes_bounded_score_logit_normal.h"
    ),
    warn = FALSE
  ),
  collapse = "\n"
)
generic_continuous_header <- paste(
  readLines(
    file.path(
      repo_root,
      "src",
      "stanExports_trunc_comp_bayes_bounded_continuous_logit_normal.h"
    ),
    warn = FALSE
  ),
  collapse = "\n"
)
stopifnot(grepl("ordered[H] mu_logit_comp[2]", application_source, fixed = TRUE))
stopifnot(grepl("vector[H] mu_logit_comp[2]", generic_score_source, fixed = TRUE))
stopifnot(grepl("vector[H] mu_logit_comp[2]", generic_continuous_source, fixed = TRUE))
stopifnot(!grepl("read_constrain_ordered", generic_score_header, fixed = TRUE))
stopifnot(!grepl("read_constrain_ordered", generic_continuous_header, fixed = TRUE))
stopifnot(identical(
  application_bayes_model_code_fingerprint(application_source),
  contract$model_code_sha256
))

sampled_parameter_diagnostics <- application_bayes_sampled_parameter_diagnostics(cached$fit)
stopifnot(identical(
  sampled_parameter_diagnostics$parameter_table$variable,
  application_bayes_sampled_parameter_names()
))
stopifnot(nrow(sampled_parameter_diagnostics$parameter_table) == 17L)
stopifnot(identical(format_fixed(sampled_parameter_diagnostics$max_rhat, 4), "1.0018"))
stopifnot(identical(format_count(sampled_parameter_diagnostics$min_bulk_ess), "3,584"))
stopifnot(identical(format_count(sampled_parameter_diagnostics$min_tail_ess), "4,169"))

diagnostic_macros <- .application_bayes_diagnostic_macros(cached)
stopifnot(identical(diagnostic_macros$ApplicationBayesMaxRhat, "1.0018"))
stopifnot(identical(diagnostic_macros$ApplicationBayesMinBulkEss, "3,584"))
stopifnot(identical(diagnostic_macros$ApplicationBayesMinTailEss, "4,169"))

tampered_settings <- cached
tampered_settings$settings$chains <- 3L
settings_validation <- application_bayes_validate_cache(
  cached = tampered_settings,
  d = data,
  cache_path = cache_path,
  manuscript_dir = manuscript_dir
)
stopifnot(!settings_validation$valid)
stopifnot(any(grepl("cached application settings", settings_validation$errors, fixed = TRUE)))

tampered_data <- data
tampered_data$Y[[1L]] <- tampered_data$Y[[1L]] + 1
data_validation <- application_bayes_validate_cache(
  cached = cached,
  d = tampered_data,
  cache_path = cache_path,
  manuscript_dir = manuscript_dir
)
stopifnot(!data_validation$valid)
stopifnot(any(grepl("analysis data fingerprint", data_validation$errors, fixed = TRUE)))

old_refresh <- Sys.getenv("TRUNCCOMP_REFRESH_IST3_BAYES", unset = NA_character_)
old_skip <- Sys.getenv("TRUNCCOMP_SKIP_IST3_BAYES", unset = NA_character_)
on.exit({
  if (is.na(old_refresh)) Sys.unsetenv("TRUNCCOMP_REFRESH_IST3_BAYES") else
    Sys.setenv(TRUNCCOMP_REFRESH_IST3_BAYES = old_refresh)
  if (is.na(old_skip)) Sys.unsetenv("TRUNCCOMP_SKIP_IST3_BAYES") else
    Sys.setenv(TRUNCCOMP_SKIP_IST3_BAYES = old_skip)
}, add = TRUE)
Sys.unsetenv(c("TRUNCCOMP_REFRESH_IST3_BAYES", "TRUNCCOMP_SKIP_IST3_BAYES"))

loaded <- fit_or_load_application_bayes(
  d = data,
  metadata = list(atom = -1),
  manuscript_dir = manuscript_dir
)
stopifnot(isTRUE(loaded$from_cache))
stopifnot(identical(loaded$model_id, contract$model_id))
stopifnot(!file.exists(application_bayes_candidate_cache_path(manuscript_dir)))

valid_candidate <- loaded
valid_candidate$from_cache <- FALSE
valid_candidate$cache_path <- application_bayes_candidate_cache_path(manuscript_dir)
valid_candidate$model_id <- contract$model_id
valid_candidate$contract <- contract
valid_candidate$fit$settings$model_name <- contract$model_id
candidate_validation <- application_bayes_validate_candidate(
  candidate = valid_candidate,
  d = data,
  manuscript_dir = manuscript_dir
)
stopifnot(candidate_validation$valid)

candidate_roundtrip_path <- tempfile("ist3-valid-candidate-", fileext = ".rds")
saveRDS(valid_candidate, candidate_roundtrip_path)
roundtrip_validation <- application_bayes_validate_candidate(
  candidate = readRDS(candidate_roundtrip_path),
  d = data,
  manuscript_dir = manuscript_dir
)
unlink(candidate_roundtrip_path)
stopifnot(roundtrip_validation$valid)

invalid_injected_fit <- valid_candidate
invalid_injected_fit$fit$settings$model_name <- "trunc_comp_bayes_bounded_score_logit_normal"
invalid_candidate_validation <- application_bayes_validate_candidate(
  candidate = invalid_injected_fit,
  d = data,
  manuscript_dir = manuscript_dir
)
stopifnot(!invalid_candidate_validation$valid)
stopifnot(any(grepl(
  "fitted model identifier",
  invalid_candidate_validation$errors,
  fixed = TRUE
)))
stopifnot(!file.exists(application_bayes_candidate_cache_path(manuscript_dir)))

missing_cache_manuscript <- tempfile("ist3-contract-missing-cache-")
dir.create(file.path(missing_cache_manuscript, "application-data"), recursive = TRUE)
missing_error <- tryCatch(
  fit_or_load_application_bayes(
    d = data,
    metadata = list(atom = -1),
    manuscript_dir = missing_cache_manuscript
  ),
  error = identity
)
stopifnot(inherits(missing_error, "error"))
stopifnot(grepl("No Bayesian fit was started", conditionMessage(missing_error), fixed = TRUE))

stopifnot(identical(
  .application_sha256_file(cache_path),
  cache_sha256_before
))
cat("IST-3 Bayesian application contract tests passed.\n")
