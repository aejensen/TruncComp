args <- commandArgs(trailingOnly = TRUE)

.parse_flag <- function(arguments, flag, default = NULL) {
  match <- grep(paste0("^", flag, "="), arguments, value = TRUE)
  if (!length(match)) {
    return(default)
  }

  sub(paste0("^", flag, "="), "", match[[1]])
}

.parse_csv_integers <- function(value) {
  if (is.null(value) || !nzchar(value)) {
    return(NULL)
  }

  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

.parse_csv_strings <- function(value) {
  if (is.null(value) || !nzchar(value)) {
    return(NULL)
  }

  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) {
  stop("Unable to determine the script path.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
study_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
repo_root <- normalizePath(file.path(study_dir, ".."), mustWork = TRUE)

source(file.path(repo_root, "manuscript", "R", "utils.R"), local = globalenv())
source(file.path(study_dir, "R", "simulation-study.R"), local = globalenv())
source(file.path(study_dir, "R", "ci-study.R"), local = globalenv())
source(file.path(study_dir, "R", "slurm.R"), local = globalenv())

output_dir <- .parse_flag(args, "--output-dir", ci_study_output_dir(repo_root))
reps <- as.integer(.parse_flag(args, "--reps", "10000"))
chunk_reps <- as.integer(.parse_flag(args, "--chunk-reps", "500"))
scenario_ids <- .parse_csv_strings(.parse_flag(args, "--scenarios", NULL))
n_seq <- .parse_csv_integers(.parse_flag(args, "--n", NULL))
effect_levels <- .parse_csv_integers(.parse_flag(args, "--h", NULL))
conf_level <- as.numeric(.parse_flag(args, "--conf-level", "0.95"))
joint_resolution <- as.integer(.parse_flag(args, "--joint-resolution", "35"))
overwrite <- "--overwrite" %in% args
dry_run <- "--dry-run" %in% args

partition <- .parse_flag(args, "--partition", "standard")
time_limit <- .parse_flag(args, "--time", "04:00:00")
cpus_per_task <- as.integer(.parse_flag(args, "--cpus-per-task", "1"))
mem <- .parse_flag(args, "--mem", "4G")
array_parallelism <- as.integer(.parse_flag(args, "--array-parallelism", "80"))
collector_time <- .parse_flag(args, "--collector-time", "00:30:00")
collector_mem <- .parse_flag(args, "--collector-mem", "8G")
job_name_prefix <- .parse_flag(args, "--job-name-prefix", "trunccomp2_ci")

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

if (!nrow(prepared$pending_manifest)) {
  message("No pending CI-study chunks matched the requested filters.")
  if (!dry_run) {
    results <- aggregate_ci_study_results(prepared$output_dir, config = prepared$config)
    ci_study_finalize_results(
      output_dir = prepared$output_dir,
      ci_results = results
    )
  }
  quit(save = "no", status = 0L)
}

manifest_path <- ci_study_write_manifest(prepared$pending_manifest, prepared$output_dir)
submission <- ci_study_submit_slurm_jobs(
  repo_root = repo_root,
  output_dir = prepared$output_dir,
  manifest_path = manifest_path,
  chunk_count = nrow(prepared$pending_manifest),
  partition = partition,
  time = time_limit,
  cpus_per_task = cpus_per_task,
  mem = mem,
  array_parallelism = array_parallelism,
  collector_time = collector_time,
  collector_mem = collector_mem,
  job_name_prefix = job_name_prefix,
  dry_run = dry_run
)

if (dry_run) {
  message(sprintf("Dry run only. Prepared %d pending CI-study chunk(s).", submission$chunk_count))
  message(sprintf("Manifest written to %s", submission$manifest_path))
  message(submission$array_command)
  message(submission$collector_command)
} else {
  message(sprintf("Submitted %d CI-study chunk(s) via Slurm.", submission$chunk_count))
  message(sprintf("Array job id: %s", submission$array_job_id))
  message(sprintf("Collector job id: %s", submission$collector_job_id))
  message(sprintf("Manifest written to %s", submission$manifest_path))
  message(sprintf("Submission details written to %s", ci_study_slurm_submission_path(prepared$output_dir)))
}
