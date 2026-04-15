script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if(length(script_arg) > 0) {
  normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath("tools/generate-el-fixture.R", winslash = "/", mustWork = FALSE)
}
package_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
output_path <- if(length(args) > 0) args[1] else {
  file.path(package_root, "tests", "testthat", "fixtures", "el_means_reference.rds")
}
output_path <- normalizePath(output_path, winslash = "/", mustWork = FALSE)

old_wd <- getwd()
setwd(package_root)
on.exit(setwd(old_wd), add = TRUE)

lib_dir <- tempfile("el-lib-")
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
install.packages("EL", repos = "https://cran.r-project.org", lib = lib_dir, quiet = TRUE)

old_lib_paths <- .libPaths(c(lib_dir, .libPaths()))
on.exit(.libPaths(old_lib_paths), add = TRUE)

TruncComp2Example <- readRDS(file.path("inst", "extdata", "TruncComp2Example.rds"))
library(EL)

cases <- list(
  moderate_shift = list(
    x = c(0.138, 0.926, -0.555, 0.055, 1.813, 0.697, 0.597, -0.052),
    y = c(-0.476, 0.780, 0.056, 2.265, 0.953, 0.000, -0.923)
  ),
  equal_location = list(
    x = c(-1.1, -0.5, -0.2, 0.1, 0.4, 0.9, 1.2),
    y = c(-1.3, -0.8, -0.1, 0.0, 0.2, 0.8, 1.0, 1.4)
  ),
  heavy_ties = list(
    x = c(0.0, 0.0, 2.0, 1.0, 1.5, 2.5, 2.5, 2.5),
    y = c(0.0, 0.5, 0.0, 2.5, 0.0, 0.0, 1.0, 0.0)
  ),
  small_sample = list(
    x = c(-0.033, 0.692, 0.069, 0.213),
    y = c(-0.675, 0.941, -0.619)
  ),
  identical_constants = list(
    x = c(1, 1, 1),
    y = c(1, 1, 1)
  ),
  separated_constants = list(
    x = c(1, 1, 1),
    y = c(2, 2, 2)
  ),
  trunccomp_example_alive = with(TruncComp2Example, list(
    x = Y[R == 1 & Y != 0],
    y = Y[R == 0 & Y != 0]
  ))
)

fixture <- list(
  generated_on = as.Date("2026-04-15"),
  upstream_package = "EL",
  upstream_version = as.character(utils::packageVersion("EL")),
  cases = lapply(names(cases), function(name) {
    x <- cases[[name]]$x
    y <- cases[[name]]$y
    fit <- suppressWarnings(EL.means(x, y))
    list(
      name = name,
      x = x,
      y = y,
      reference = list(
        estimate = as.numeric(fit$estimate),
        conf.int = as.numeric(fit$conf.int),
        statistic = as.numeric(fit$statistic),
        p.value = as.numeric(fit$p.value)
      )
    )
  })
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(fixture, output_path, version = 2)
