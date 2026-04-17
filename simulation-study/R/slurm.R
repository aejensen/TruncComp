simulation_study_slurm_submission_path <- function(output_dir) {
  file.path(simulation_study_slurm_dir(output_dir), "latest-submission.txt")
}

.simulation_study_shell_command <- function(command, args) {
  paste(c(command, shQuote(args)), collapse = " ")
}

.simulation_study_submit_sbatch <- function(args, dry_run = FALSE) {
  rendered <- .simulation_study_shell_command("sbatch", args)
  if (dry_run) {
    return(list(job_id = NA_character_, command = rendered, stdout = character()))
  }

  output <- system2("sbatch", args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(paste(output, collapse = "\n"), call. = FALSE)
  }

  if (!length(output)) {
    stop("sbatch returned no output.", call. = FALSE)
  }

  job_id <- sub(";.*$", "", trimws(output[[length(output)]]))
  if (!nzchar(job_id)) {
    stop("Unable to parse the sbatch job id.", call. = FALSE)
  }

  list(job_id = job_id, command = rendered, stdout = output)
}

.simulation_study_array_spec <- function(cell_count, array_parallelism = NULL) {
  if (!(length(cell_count) == 1L && is.numeric(cell_count) && is.finite(cell_count) && cell_count >= 1)) {
    stop("cell_count must be a single positive integer.", call. = FALSE)
  }

  cell_count <- as.integer(cell_count)
  if (is.null(array_parallelism) || !nzchar(as.character(array_parallelism))) {
    return(sprintf("1-%d", cell_count))
  }

  if (!(length(array_parallelism) == 1L && is.numeric(array_parallelism) && is.finite(array_parallelism) &&
        array_parallelism >= 1)) {
    stop("array_parallelism must be a single positive integer when provided.", call. = FALSE)
  }

  sprintf("1-%d%%%d", cell_count, as.integer(array_parallelism))
}

simulation_study_submit_slurm_jobs <- function(repo_root,
                                               output_dir,
                                               manifest_path,
                                               cell_count,
                                               partition = "standard",
                                               time = "04:00:00",
                                               cpus_per_task = 1L,
                                               mem = "4G",
                                               array_parallelism = 4L,
                                               collector_time = "00:15:00",
                                               collector_mem = "1G",
                                               job_name_prefix = "trunccomp2_sim",
                                               dry_run = FALSE) {
  if (!(length(cpus_per_task) == 1L && is.numeric(cpus_per_task) && is.finite(cpus_per_task) &&
        cpus_per_task >= 1)) {
    stop("cpus_per_task must be a single positive integer.", call. = FALSE)
  }

  study_dir <- simulation_study_root(repo_root)
  cell_script <- file.path(study_dir, "slurm", "run-simulation-study-cell.sbatch")
  collect_script <- file.path(study_dir, "slurm", "collect-simulation-study-results.sbatch")
  if (!file.exists(cell_script)) {
    stop("Missing Slurm cell script at ", cell_script, call. = FALSE)
  }
  if (!file.exists(collect_script)) {
    stop("Missing Slurm collector script at ", collect_script, call. = FALSE)
  }

  slurm_dir <- simulation_study_slurm_dir(output_dir)
  array_spec <- .simulation_study_array_spec(cell_count, array_parallelism = array_parallelism)

  cell_args <- c(
    "--parsable",
    paste0("--job-name=", job_name_prefix, "_cells"),
    paste0("--partition=", partition),
    paste0("--time=", time),
    paste0("--cpus-per-task=", as.integer(cpus_per_task)),
    paste0("--array=", array_spec),
    paste0("--output=", file.path(slurm_dir, "cell-%A_%a.out")),
    paste0("--error=", file.path(slurm_dir, "cell-%A_%a.err")),
    paste0("--chdir=", repo_root)
  )
  if (!is.null(mem) && nzchar(mem)) {
    cell_args <- c(cell_args, paste0("--mem=", mem))
  }
  cell_args <- c(cell_args, cell_script, repo_root, output_dir, manifest_path)

  cell_submission <- .simulation_study_submit_sbatch(cell_args, dry_run = dry_run)
  dependency_job_id <- if (dry_run) "<array-job-id>" else cell_submission$job_id

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
    sprintf("cell_count=%d", as.integer(cell_count)),
    sprintf("array_job_id=%s", if (dry_run) "<dry-run>" else cell_submission$job_id),
    sprintf("collector_job_id=%s", if (dry_run) "<dry-run>" else collect_submission$job_id),
    sprintf("array_command=%s", cell_submission$command),
    sprintf("collector_command=%s", collect_submission$command)
  )
  writeLines(metadata, con = simulation_study_slurm_submission_path(output_dir), useBytes = TRUE)

  list(
    manifest_path = manifest_path,
    cell_count = as.integer(cell_count),
    array_job_id = cell_submission$job_id,
    collector_job_id = collect_submission$job_id,
    array_command = cell_submission$command,
    collector_command = collect_submission$command,
    dry_run = dry_run
  )
}
