test_that("Successful Bayesian fits return the expected object shape", {
  fit <- bayes_formula_fit()

  expect_s3_class(fit, "trunc_comp_bayes_fit")
  expect_true(fit$success)
  expect_true(inherits(fit$fit, "stanfit"))

  expect_true(all(c(
    "fit", "draws", "summary_table", "arm_table", "diagnostics",
    "ppc_table", "ppc_settings",
    "settings", "conf.level", "success", "error", "data", "atom", "call"
  ) %in% names(fit)))

  expect_true(all(c(
    "auto_select_mixture_components",
    "mixture_components_initial",
    "mixture_components_final",
    "mixture_components_max",
    "mixture_component_path",
    "mixture_selection_history"
  ) %in% names(fit$settings)))
  expect_false(fit$settings$auto_select_mixture_components)
  expect_equal(fit$settings$mixture_components_initial, 3L)
  expect_equal(fit$settings$mixture_components_final, 3L)
  expect_equal(fit$settings$mixture_component_path, 3L)
  expect_s3_class(fit$settings$mixture_selection_history, "data.frame")

  expect_true(all(TruncComp2:::bayes_parameter_names("all") %in% names(fit$draws)))
  expect_gt(nrow(fit$draws), 0)

  contrast_estimates <- fit$summary_table[c("delta_atom", "mu_delta", "alpha_delta", "delta"), "estimate"]
  arm_estimates <- fit$arm_table[c("rho_0", "rho_1", "mu_0_c", "mu_1_c"), "estimate"]

  expect_true(all(is.finite(contrast_estimates)))
  expect_true(all(is.finite(arm_estimates)))
  expect_equal(rownames(fit$ppc_table), c("atom", "continuous"))
  expect_true(all(is.finite(fit$ppc_table$p_value)))
  expect_true(all(c("core_ok", "truncation_ok", "truncation", "parameter_table") %in% names(fit$diagnostics)))
  expect_s3_class(fit$diagnostics$truncation$parameter_table, "data.frame")
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
  expect_equal(rownames(fit$ppc_table), c("atom", "continuous"))
  expect_true(all(is.finite(fit$ppc_table$p_value)))
  expect_true(all(c("core_ok", "truncation_ok", "truncation", "parameter_table") %in% names(fit$diagnostics)))
})

test_that("Bounded Bayesian fits return expected settings and draw names", {
  continuous_fit <- bayes_bounded_continuous_formula_fit(seed = 2301)
  score_fit <- bayes_bounded_score_formula_fit(seed = 2302)

  expect_s3_class(continuous_fit, "trunc_comp_bayes_fit")
  expect_s3_class(score_fit, "trunc_comp_bayes_fit")
  expect_true(continuous_fit$success)
  expect_true(score_fit$success)

  expect_equal(continuous_fit$settings$continuous_support, "bounded_continuous")
  expect_equal(score_fit$settings$continuous_support, "bounded_score")
  expect_equal(continuous_fit$settings$bounded_kernel, "beta")
  expect_equal(score_fit$settings$bounded_kernel, "beta")
  expect_equal(continuous_fit$settings$score_min, 0)
  expect_equal(continuous_fit$settings$score_max, 100)
  expect_equal(score_fit$settings$score_step, 1)
  expect_equal(score_fit$settings$heaping, "shared")
  expect_true(all(c("score_values", "bin_lower", "bin_upper", "bin_valid") %in% names(score_fit$settings)))

  expect_true(all(TruncComp2:::bayes_parameter_names("all") %in% names(continuous_fit$draws)))
  expect_true(all(TruncComp2:::bayes_parameter_names("all") %in% names(score_fit$draws)))
  expect_true(all(is.finite(continuous_fit$summary_table$estimate)))
  expect_true(all(is.finite(score_fit$summary_table$estimate)))
  expect_equal(rownames(continuous_fit$ppc_table), c("atom", "continuous"))
  expect_equal(rownames(score_fit$ppc_table), c("atom", "continuous"))
})
