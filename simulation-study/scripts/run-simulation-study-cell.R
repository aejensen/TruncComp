args <- commandArgs(trailingOnly = TRUE)

.parse_flag <- function(arguments, flag, default = NULL) {
  match <- grep(paste0("^", flag, "="), arguments, value = TRUE)
  if (!length(match)) {
    return(default)
  }

  sub(paste0("^", flag, "="), "", match[[1]])
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

output_dir <- .parse_flag(args, "--output-dir", simulation_study_output_dir(repo_root))
manifest_path <- .parse_flag(args, "--manifest", simulation_study_manifest_path(output_dir))
task_id <- .parse_flag(args, "--task-id", Sys.getenv("SLURM_ARRAY_TASK_ID", unset = ""))

if (!nzchar(task_id)) {
  stop("A task id is required via --task-id or SLURM_ARRAY_TASK_ID.", call. = FALSE)
}

task_id <- as.integer(task_id)
if (!is.finite(task_id) || task_id < 1L) {
  stop("task id must be a positive integer.", call. = FALSE)
}

if (!file.exists(manifest_path)) {
  stop("No simulation-study manifest found at ", manifest_path, call. = FALSE)
}

manifest <- readRDS(manifest_path)
if (!(is.data.frame(manifest) && nrow(manifest) >= task_id)) {
  stop(
    sprintf("Task id %d is outside the manifest bounds (n = %d).", task_id, if (is.data.frame(manifest)) nrow(manifest) else 0L),
    call. = FALSE
  )
}

cell <- manifest[task_id, , drop = FALSE]

load_local_trunccomp2(repo_root)
message(sprintf("Running manifest task %d: %s, h = %d, n = %d", task_id, cell$scenario_id, cell$h, cell$n))
simulation_study_run_cell(cell, output_dir = output_dir)
