script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if(length(script_arg) > 0) {
  normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath("tools/validate-bayesian-dp.R", winslash = "/", mustWork = FALSE)
}
package_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

library(TruncComp2)

simulate_validation_case <- function(seed, n, f0, f1, pi0, pi1, atom = 0) {
  for(offset in 0:1000) {
    set.seed(seed + offset)
    simulated <- simulate_truncated_data(n, f0 = f0, f1 = f1, pi0 = pi0, pi1 = pi1, atom = atom)
    if(sum(simulated$A[simulated$R == 0]) >= 3 &&
       sum(simulated$A[simulated$R == 1]) >= 3 &&
       any(simulated$A == 0)) {
      return(simulated)
    }
  }

  stop("Unable to generate validation data.")
}

fit_validation_case <- function(data, seed) {
  bayes_fit <- suppressWarnings(trunc_comp_bayes(
    Y ~ R,
    atom = 0,
    data = data,
    mixture_components = 5,
    chains = 4,
    iter_warmup = 400,
    iter_sampling = 400,
    seed = seed,
    refresh = 0,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    cores = 1
  ))

  lrt_fit <- trunc_comp(Y ~ R, atom = 0, data = data, method = "lrt")
  splrt_fit <- trunc_comp(Y ~ R, atom = 0, data = data, method = "splrt")

  list(
    bayes = bayes_fit,
    lrt = lrt_fit,
    splrt = splrt_fit
  )
}

scenarios <- list(
  null = function() simulate_validation_case(1, 30, function(n) stats::rnorm(n, 1.4, 0.5), function(n) stats::rnorm(n, 1.4, 0.5), 0.7, 0.7),
  atom_only = function() simulate_validation_case(2, 30, function(n) stats::rnorm(n, 1.4, 0.5), function(n) stats::rnorm(n, 1.4, 0.5), 0.5, 0.8),
  continuous_only = function() simulate_validation_case(3, 30, function(n) stats::rnorm(n, 1.2, 0.4), function(n) stats::rnorm(n, 1.8, 0.4), 0.7, 0.7),
  opposing = function() simulate_validation_case(4, 30, function(n) stats::rnorm(n, 1.7, 0.4), function(n) stats::rnorm(n, 1.3, 0.4), 0.55, 0.8),
  bimodal = function() simulate_validation_case(5, 30, function(n) c(stats::rnorm(ceiling(n / 2), 0.8, 0.2), stats::rnorm(floor(n / 2), 1.8, 0.2)), function(n) c(stats::rnorm(ceiling(n / 2), 1.2, 0.2), stats::rnorm(floor(n / 2), 2.2, 0.2)), 0.7, 0.7),
  skewed = function() simulate_validation_case(6, 30, function(n) stats::rgamma(n, shape = 4, rate = 3), function(n) stats::rgamma(n, shape = 6, rate = 3), 0.7, 0.7),
  small_sample = function() simulate_validation_case(7, 12, function(n) stats::rnorm(n, 1.3, 0.5), function(n) stats::rnorm(n, 1.8, 0.5), 0.7, 0.65),
  large_sample = function() simulate_validation_case(8, 60, function(n) stats::rnorm(n, 1.3, 0.4), function(n) stats::rnorm(n, 1.9, 0.4), 0.65, 0.8)
)

results <- lapply(seq_along(scenarios), function(i) {
  name <- names(scenarios)[[i]]
  data <- scenarios[[i]]()
  fits <- fit_validation_case(data, seed = 900 + i)

  bayes_coef <- coef(fits$bayes)
  lrt_coef <- coef(fits$lrt)
  splrt_coef <- coef(fits$splrt)

  cat("\nScenario:", name, "\n")
  cat("  Bayesian delta:", format(bayes_coef["delta"], digits = 4), "\n")
  cat("  Bayesian mu_delta:", format(bayes_coef["mu_delta"], digits = 4), "\n")
  cat("  LRT delta:", format(lrt_coef["delta"], digits = 4), "\n")
  cat("  SPLRT delta:", format(splrt_coef["delta"], digits = 4), "\n")
  cat("  Divergences:", fits$bayes$diagnostics$divergences, "\n")
  cat("  Max Rhat:", format(fits$bayes$diagnostics$max_rhat, digits = 4), "\n")

  list(name = name, data = data, fits = fits)
})

reference_results <- Filter(function(x) x$name %in% c("null", "atom_only", "continuous_only", "opposing", "large_sample"), results)

for(result in reference_results) {
  diagnostics <- result$fits$bayes$diagnostics
  if(diagnostics$divergences > 0L) {
    stop("Validation failed: divergences detected in scenario ", result$name, ".", call. = FALSE)
  }
  if(!is.finite(diagnostics$max_rhat) || diagnostics$max_rhat > 1.01) {
    stop("Validation failed: max Rhat exceeds 1.01 in scenario ", result$name, ".", call. = FALSE)
  }
}

cat("\nBayesian DP validation completed successfully.\n")
