test_that("Bayesian smoke fit returns finite summaries and diagnostics", {
  fit <- bayes_formula_fit(seed = 404)

  expect_true(fit$success)
  expect_gt(nrow(fit$draws), 0)
  expect_true(is.list(fit$diagnostics))
  expect_true(all(c(
    "divergences", "max_rhat", "min_bulk_ess", "min_tail_ess",
    "core_ok", "truncation_ok", "diagnostic_ok", "truncation"
  ) %in% names(fit$diagnostics)))
  expect_true(all(is.finite(fit$summary_table$estimate)))
  expect_true(all(is.finite(fit$arm_table$estimate)))
})

test_that("Positive-support Bayesian smoke fit returns finite summaries and diagnostics", {
  fit <- bayes_positive_formula_fit(seed = 1404)

  expect_true(fit$success)
  expect_equal(fit$settings$continuous_support, "positive_real")
  expect_gt(nrow(fit$draws), 0)
  expect_true(all(c("core_ok", "truncation_ok", "diagnostic_ok", "truncation") %in% names(fit$diagnostics)))
  expect_true(all(is.finite(fit$summary_table$estimate)))
  expect_true(all(is.finite(fit$arm_table$estimate)))
})
