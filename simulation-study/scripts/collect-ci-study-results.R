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
source(file.path(study_dir, "R", "ci-study.R"), local = globalenv())

output_dir <- .parse_flag(args, "--output-dir", ci_study_output_dir(repo_root))
results <- aggregate_ci_study_results(output_dir)

ci_study_finalize_results(
  output_dir = output_dir,
  ci_results = results
)
