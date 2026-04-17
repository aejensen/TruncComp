script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if(length(script_arg) > 0) {
  normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath("tools/build-package-data.R", winslash = "/", mustWork = FALSE)
}
package_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

generate_adjusted_example <- function() {
  set.seed(266)

  n_per_arm <- 25L
  levels_L <- c("low", "mid", "high")
  pL0 <- c(0.55, 0.30, 0.15)
  pL1 <- c(0.20, 0.35, 0.45)
  pA <- c(0.30, 0.55, 0.75)
  muL <- c(2.0, 3.0, 4.0)
  sigma <- 0.9

  R <- c(rep.int(0, n_per_arm), rep.int(1, n_per_arm))
  L0 <- sample(levels_L, n_per_arm, replace = TRUE, prob = pL0)
  L1 <- sample(levels_L, n_per_arm, replace = TRUE, prob = pL1)
  L <- factor(c(L0, L1), levels = levels_L)

  alive <- c(
    stats::rbinom(n_per_arm, 1, pA[match(L0, levels_L)]),
    stats::rbinom(n_per_arm, 1, pA[match(L1, levels_L)])
  )

  mu <- muL[match(L, levels_L)]
  Y <- numeric(2L * n_per_arm)
  Y[alive == 1] <- stats::rnorm(sum(alive == 1), mean = mu[alive == 1], sd = sigma)

  data.frame(R = R, L = L, Y = Y)
}

trunc_comp_example <- readRDS(file.path(package_root, "inst", "extdata", "TruncComp2Example.rds"))
trunc_comp_adjusted_example <- generate_adjusted_example()

data_dir <- file.path(package_root, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

save(trunc_comp_example, file = file.path(data_dir, "trunc_comp_example.rda"), version = 2)
save(
  trunc_comp_adjusted_example,
  file = file.path(data_dir, "trunc_comp_adjusted_example.rda"),
  version = 2
)

cat("Saved package datasets in", data_dir, "\n")
