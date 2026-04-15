script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if(length(script_arg) > 0) {
  normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath("tools/check-package.R", winslash = "/", mustWork = FALSE)
}
package_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

description <- read.dcf(file.path(package_root, "DESCRIPTION"))
package_name <- description[1, "Package"]
r_bin <- file.path(R.home("bin"), "R")
library_dir <- tempfile("trunccomp2-lib-")
dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))

ensure_installed <- function(pkg) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", lib = .libPaths()[1])
  }
}

run_command <- function(args, wd = package_root) {
  old_wd <- getwd()
  setwd(wd)
  on.exit(setwd(old_wd), add = TRUE)

  status <- system2(r_bin, args = args)
  if(status != 0) {
    stop("Command failed: ", paste(c(r_bin, args), collapse = " "), call. = FALSE)
  }
}

cat("Preparing package-local verification for ", package_name, "\n", sep = "")

ensure_installed("remotes")
remotes::install_deps(
  package_root,
  dependencies = c("Depends", "Imports", "LinkingTo", "Suggests"),
  upgrade = "never",
  quiet = TRUE
)

run_command(c("CMD", "INSTALL", "-l", library_dir, package_root))

library(testthat)
library(package_name, character.only = TRUE)

cat("Running package-local testthat suite\n")
testthat::test_dir(file.path(package_root, "tests", "testthat"), reporter = "summary")

build_dir <- tempfile("trunccomp2-build-")
dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)

cat("Building source tarball\n")
run_command(c("CMD", "build", package_root), wd = build_dir)

tarball <- list.files(
  build_dir,
  pattern = paste0("^", package_name, "_.*[.]tar[.]gz$"),
  full.names = TRUE
)

if(length(tarball) != 1L) {
  stop("Expected one built tarball in ", build_dir, call. = FALSE)
}

cat("Running R CMD check\n")
run_command(c("CMD", "check", "--no-manual", "--no-build-vignettes", tarball), wd = build_dir)

cat("Verification completed successfully for ", package_name, "\n", sep = "")
