args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2 || !identical(args[[1]], "--target")) {
  stop("Usage: Rscript manuscript/scripts/build-assets.R --target manuscript", call. = FALSE)
}

target <- args[[2]]
if (!identical(target, "manuscript")) {
  stop("Only --target manuscript is currently supported.", call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) {
  stop("Unable to determine the build script path.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
manuscript_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
repo_root <- normalizePath(file.path(manuscript_dir, ".."), mustWork = TRUE)

source(file.path(manuscript_dir, "R", "utils.R"), local = globalenv())
for (file_name in c("data.R", "analysis.R", "figures.R", "tables.R", "build.R")) {
  source(file.path(manuscript_dir, "R", file_name), local = globalenv())
}

build_manuscript_assets(file.path(manuscript_dir, "build"), repo_root = repo_root)
