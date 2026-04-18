test_that("Bayesian print, summary, coef, and confint methods work on successful fits", {
  fit <- bayes_formula_fit()

  expect_match(paste(capture.output(print(fit)), collapse = "\n"), "trunc_comp_bayes_fit")

  capture.output(summary_out <- summary(fit))
  expect_s3_class(summary_out, "trunc_comp_bayes_summary")
  expect_true(all(c("contrasts", "arms", "diagnostics", "settings") %in% names(summary_out)))

  coefficients <- coef(fit)
  expect_equal(names(coefficients), c("mu_delta", "delta_atom", "alpha_delta", "delta"))
  expect_true(all(is.finite(coefficients)))

  capture.output(intervals <- confint(fit))
  expect_equal(
    rownames(intervals),
    c("delta_atom", "mu_delta", "alpha_delta", "delta")
  )

  expect_error(confint(fit, parameter = "joint"), "arg")
  expect_error(confint(fit, parameter = "delta", method = "profile"), "Bayesian credible intervals")
})

test_that("Bayesian methods behave on failed fits", {
  bad_data <- data.frame(
    Y = c(0, 1.1, 0, 1.4),
    R = c(0, 0, 1, 1)
  )

  fit <- trunc_comp_bayes(
    Y ~ R,
    atom = 0,
    data = bad_data,
    mixture_components = 3,
    chains = 2,
    iter_warmup = 20,
    iter_sampling = 20,
    seed = 77,
    refresh = 0,
    control = list(adapt_delta = 0.9, max_treedepth = 10),
    cores = 1
  )

  expect_false(fit$success)
  expect_match(paste(capture.output(print(fit)), collapse = "\n"), "Error:")
  expect_match(paste(capture.output(summary(fit)), collapse = "\n"), "failed")
  expect_equal(coef(fit), stats::setNames(rep(NA_real_, 4), c("mu_delta", "delta_atom", "alpha_delta", "delta")))
  expect_error(confint(fit), "Estimation failed")
})
