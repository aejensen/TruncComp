test_that("Successful Bayesian fits return the expected object shape", {
  fit <- bayes_formula_fit()

  expect_s3_class(fit, "trunc_comp_bayes_fit")
  expect_true(fit$success)
  expect_true(inherits(fit$fit, "stanfit"))

  expect_true(all(c(
    "fit", "draws", "summary_table", "arm_table", "diagnostics",
    "settings", "conf.level", "success", "error", "data", "atom", "call"
  ) %in% names(fit)))

  expect_true(all(TruncComp2:::bayes_parameter_names("all") %in% names(fit$draws)))
  expect_gt(nrow(fit$draws), 0)

  contrast_estimates <- fit$summary_table[c("delta_atom", "mu_delta", "alpha_delta", "delta"), "estimate"]
  arm_estimates <- fit$arm_table[c("rho_0", "rho_1", "mu_0_c", "mu_1_c"), "estimate"]

  expect_true(all(is.finite(contrast_estimates)))
  expect_true(all(is.finite(arm_estimates)))
})

test_that("Positive-support Bayesian fits return the expected object shape", {
  fit <- bayes_positive_formula_fit(seed = 1301)

  expect_s3_class(fit, "trunc_comp_bayes_fit")
  expect_true(fit$success)
  expect_equal(fit$settings$continuous_support, "positive_real")
  expect_true(inherits(fit$fit, "stanfit"))
  expect_true(all(TruncComp2:::bayes_parameter_names("all") %in% names(fit$draws)))

  contrast_estimates <- fit$summary_table[c("delta_atom", "mu_delta", "alpha_delta", "delta"), "estimate"]
  arm_estimates <- fit$arm_table[c("rho_0", "rho_1", "mu_0_c", "mu_1_c"), "estimate"]

  expect_true(all(is.finite(contrast_estimates)))
  expect_true(all(is.finite(arm_estimates)))
})
