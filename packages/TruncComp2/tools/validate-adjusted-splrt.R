#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pkgload)
})

pkgload::load_all(".", quiet = TRUE)

generate_case <- function(seed, n = 20, treatment_effect = 0, logit_treatment = 0) {
  set.seed(seed)
  total_n <- 2 * n
  L1 <- stats::rnorm(total_n)
  L2 <- factor(sample(c("a", "b", "c"),
                      total_n,
                      replace = TRUE,
                      prob = c(0.45, 0.35, 0.20)))
  R <- c(rep(0, n), rep(1, n))
  A <- stats::rbinom(total_n,
                     1,
                     stats::plogis(-0.7 + logit_treatment * R + 0.6 * L1 - 0.5 * (L2 == "b") + 0.8 * (L2 == "c")))
  Y <- numeric(total_n)
  mu <- 1.5 + treatment_effect * R + 0.5 * L1 + 0.4 * (L2 == "b") + 0.9 * (L2 == "c")
  Y[A == 1] <- stats::rnorm(sum(A), mu[A == 1], 0.75)
  data.frame(Y = Y, R = R, L1 = L1, L2 = L2)
}

run_scenario <- function(label, seeds, treatment_effect, logit_treatment) {
  fits <- lapply(seeds, function(seed) {
    data <- generate_case(seed, treatment_effect = treatment_effect, logit_treatment = logit_treatment)
    truncComp(Y ~ R, atom = 0, data = data, method = "SPLRT", adjust = ~ L1 + L2)
  })

  successful <- vapply(fits, function(fit) isTRUE(fit$success), logical(1))
  p_values <- vapply(fits[successful], function(fit) fit$p, numeric(1))

  cat("\nScenario:", label, "\n")
  cat("Successful fits:", sum(successful), "of", length(fits), "\n")
  if(length(p_values) > 0) {
    cat("Mean p-value:", mean(p_values), "\n")
    cat("Rejection rate at 5%:", mean(p_values < 0.05), "\n")
  }
}

null_seeds <- 5001:5100
alt_seeds <- 6001:6100

cat("Validating adjusted SPLRT with fixed seeds.\n")
run_scenario("Approximate null calibration", null_seeds, treatment_effect = 0, logit_treatment = 0)
run_scenario("Basic power under conditional effects", alt_seeds, treatment_effect = 0.45, logit_treatment = 0.35)
