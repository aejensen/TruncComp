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

output_dir <- .parse_flag(args, "--output-dir", simulation_study_output_dir(repo_root))
reps <- as.integer(.parse_flag(args, "--reps", "10000"))
workers <- as.integer(.parse_flag(args, "--workers", "1"))
scenario_ids <- .parse_csv_strings(.parse_flag(args, "--scenarios", NULL))
n_seq <- .parse_csv_integers(.parse_flag(args, "--n", NULL))
effect_levels <- .parse_csv_integers(.parse_flag(args, "--h", NULL))
overwrite <- "--overwrite" %in% args

results <- run_simulation_study(
  repo_root = repo_root,
  output_dir = output_dir,
  reps = reps,
  workers = workers,
  scenario_ids = scenario_ids,
  n_seq = n_seq,
  effect_levels = effect_levels,
  overwrite = overwrite
)

message(
  sprintf(
    "Simulation study complete: %d/%d cells available. Aggregated results written to %s",
    results$completed_cells,
    results$expected_cells,
    file.path(output_dir, "simulation-study.rds")
  )
)
