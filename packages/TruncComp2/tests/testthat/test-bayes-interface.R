test_that("Bayesian formula and default interfaces agree under a fixed seed", {
  data <- bayes_test_data()

  fit_formula <- bayes_formula_fit(data = data, seed = 101)
  fit_default <- bayes_default_fit(data = data, seed = 101)

  expect_true(fit_formula$success)
  expect_true(fit_default$success)
  expect_equal(fit_formula$atom, fit_default$atom)
  expect_equal(coef(fit_formula), coef(fit_default), tolerance = 1e-8)
})

test_that("Bayesian default interface infers the atom when omitted", {
  data <- bayes_test_data()

  fit <- suppressWarnings(trunc_comp_bayes(
    data$Y,
    data$A,
    data$R,
    mixture_components = 3,
    chains = 2,
    iter_warmup = 60,
    iter_sampling = 60,
    seed = 202,
    refresh = 0,
    control = list(adapt_delta = 0.9, max_treedepth = 10),
    cores = 1
  ))

  expect_true(fit$success)
  expect_equal(fit$atom, 0)
})

test_that("Bayesian interface rejects unsupported inputs and invalid arguments", {
  data <- bayes_test_data()
  positive_data <- bayes_positive_test_data()
  nonbinary_data <- data
  nonbinary_data$R[[nrow(nonbinary_data)]] <- 2

  expect_error(
    trunc_comp_bayes(Y ~ R, atom = 0, data = nonbinary_data),
    "binary"
  )

  expect_error(
    trunc_comp_bayes(Y ~ R, atom = 0, data = transform(data, Y = abs(Y) + 1)),
    "Everything seems to have been observed"
  )

  missing_data <- data
  missing_data$Y[[1]] <- NA_real_
  expect_error(
    trunc_comp_bayes(Y ~ R, atom = 0, data = missing_data),
    "missing values"
  )

  expect_error(
    trunc_comp_bayes(Y ~ R, atom = 0, data = cbind(data, L = 1), adjust = ~ L),
    "Covariate adjustment is not implemented"
  )

  expect_error(
    trunc_comp_bayes(Y ~ R, atom = 0, data = data, mixture_components = 1),
    "mixture_components"
  )

  expect_error(
    trunc_comp_bayes(Y ~ R, atom = 0, data = data, iter_sampling = 0),
    "iter_sampling"
  )

  expect_error(
    trunc_comp_bayes(Y ~ R, atom = 0, data = data, continuous_support = "bad")
  )

  expect_error(
    trunc_comp_bayes(
      Y ~ R,
      atom = 0,
      data = transform(positive_data, Y = ifelse(A == 1L, -abs(Y), Y)),
      continuous_support = "positive_real"
    ),
    "strictly positive"
  )

  expect_error(
    trunc_comp_bayes(
      Y ~ R,
      atom = 0,
      data = data,
      prior = list(mean_meanlog = 0)
    ),
    "Unsupported prior fields"
  )

  expect_error(
    trunc_comp_bayes(
      Y ~ R,
      atom = 0,
      data = positive_data,
      continuous_support = "positive_real",
      prior = list(mu_mean = 0)
    ),
    "Unsupported prior fields"
  )
})

test_that("positive_real support is stored and real_line remains the default", {
  real_line_fit <- bayes_formula_fit(seed = 1201)
  positive_fit <- bayes_positive_formula_fit(seed = 1202)

  expect_true(real_line_fit$success)
  expect_true(positive_fit$success)
  expect_equal(real_line_fit$settings$continuous_support, "real_line")
  expect_equal(positive_fit$settings$continuous_support, "positive_real")
})
