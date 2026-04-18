bayes_test_data <- function(seed = 20260418, n = 10) {
  for(offset in 0:1000) {
    set.seed(seed + offset)
    simulated <- simulate_truncated_data(
      n = n,
      f0 = function(m) stats::rnorm(m, 1.2, 0.35),
      f1 = function(m) stats::rnorm(m, 1.7, 0.35),
      pi0 = 0.75,
      pi1 = 0.65,
      atom = 0
    )

    if(sum(simulated$A[simulated$R == 0]) >= 3 &&
       sum(simulated$A[simulated$R == 1]) >= 3 &&
       any(simulated$A == 0)) {
      return(simulated)
    }
  }

  stop("Unable to generate a stable Bayesian test dataset.")
}

bayes_positive_test_data <- function(seed = 20260418, n = 10) {
  for(offset in 0:1000) {
    set.seed(seed + offset)
    simulated <- simulate_truncated_data(
      n = n,
      f0 = function(m) stats::rgamma(m, shape = 4, rate = 2.5),
      f1 = function(m) stats::rgamma(m, shape = 5, rate = 2.2),
      pi0 = 0.75,
      pi1 = 0.65,
      atom = 0
    )

    if(sum(simulated$A[simulated$R == 0]) >= 3 &&
       sum(simulated$A[simulated$R == 1]) >= 3 &&
       any(simulated$A == 0) &&
       all(simulated$Y[simulated$A == 1] > 0)) {
      return(simulated)
    }
  }

  stop("Unable to generate a stable positive-support Bayesian test dataset.")
}

bayes_test_fit_args <- function(seed = 42L) {
  list(
    atom = 0,
    mixture_components = 3,
    auto_select_mixture_components = FALSE,
    chains = 2,
    iter_warmup = 60,
    iter_sampling = 60,
    seed = seed,
    refresh = 0,
    control = list(adapt_delta = 0.9, max_treedepth = 10),
    cores = 1
  )
}

bayes_formula_fit <- function(data = bayes_test_data(),
                              seed = 42L,
                              continuous_support = "real_line",
                              ...) {
  suppressWarnings(do.call(
    trunc_comp_bayes,
    c(
      list(Y ~ R, data = data, continuous_support = continuous_support),
      bayes_test_fit_args(seed = seed),
      list(...)
    )
  ))
}

bayes_default_fit <- function(data = bayes_test_data(),
                              seed = 42L,
                              continuous_support = "real_line",
                              ...) {
  suppressWarnings(do.call(
    trunc_comp_bayes,
    c(
      list(data$Y, data$A, data$R, continuous_support = continuous_support),
      bayes_test_fit_args(seed = seed),
      list(...)
    )
  ))
}

bayes_positive_formula_fit <- function(data = bayes_positive_test_data(),
                                       seed = 42L,
                                       ...) {
  bayes_formula_fit(
    data = data,
    seed = seed,
    continuous_support = "positive_real",
    ...
  )
}
