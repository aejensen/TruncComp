bayes_ppc_fit <- bayes_formula_fit(seed = 606)
bayes_positive_ppc_fit <- bayes_positive_formula_fit(seed = 1606)

test_that("posterior_predictive_check returns ggplots for the requested PPC type", {
  both_plots <- posterior_predictive_check(bayes_ppc_fit, seed = 1)
  atom_plot <- posterior_predictive_check(bayes_ppc_fit, type = "atom", seed = 1)
  continuous_plot <- posterior_predictive_check(bayes_ppc_fit, type = "continuous", seed = 1)

  expect_type(both_plots, "list")
  expect_named(both_plots, c("atom", "continuous"))
  expect_s3_class(both_plots$atom, "ggplot")
  expect_s3_class(both_plots$continuous, "ggplot")
  expect_s3_class(atom_plot, "ggplot")
  expect_s3_class(continuous_plot, "ggplot")
  expect_equal(both_plots$continuous$labels$y, "Density")
  expect_equal(both_plots$continuous$labels$x, "Outcome")
})

test_that("Bayesian PPC helper returns replicated draws with the expected shapes", {
  ppc_data <- TruncComp2:::bayes_ppc_data(
    bayes_ppc_fit,
    ndraws = 7,
    seed = 11
  )

  expect_equal(dim(ppc_data$yrep_atom), c(7, nrow(bayes_ppc_fit$data)))
  expect_equal(dim(ppc_data$yrep_cont), c(7, sum(bayes_ppc_fit$data$A == 1)))
  expect_true(all(ppc_data$yrep_atom %in% c(0, 1)))
  expect_true(all(is.finite(ppc_data$yrep_cont)))
  expect_true(all(ppc_data$group_atom %in% c("Control", "Treatment")))
  expect_true(all(ppc_data$group_cont %in% c("Control", "Treatment")))
})

test_that("Bayesian PPC helper is reproducible with a fixed seed and clips ndraws", {
  data_first <- TruncComp2:::bayes_ppc_data(
    bayes_ppc_fit,
    ndraws = 8,
    seed = 99
  )
  data_second <- TruncComp2:::bayes_ppc_data(
    bayes_ppc_fit,
    ndraws = 8,
    seed = 99
  )

  expect_equal(data_first$draw_indices, data_second$draw_indices)
  expect_equal(data_first$yrep_atom, data_second$yrep_atom)
  expect_equal(data_first$yrep_cont, data_second$yrep_cont)

  clipped <- TruncComp2:::bayes_ppc_data(
    bayes_ppc_fit,
    ndraws = nrow(bayes_ppc_fit$draws) + 100,
    seed = 123
  )

  expect_equal(nrow(clipped$yrep_atom), nrow(bayes_ppc_fit$draws))
  expect_equal(nrow(clipped$yrep_cont), nrow(bayes_ppc_fit$draws))
})

test_that("posterior_predictive_check validates failed fits and invalid inputs", {
  failed_fit <- bayes_ppc_fit
  failed_fit$success <- FALSE
  expect_error(posterior_predictive_check(failed_fit), "Estimation failed")

  missing_fit <- bayes_ppc_fit
  missing_fit$fit <- NULL
  expect_error(posterior_predictive_check(missing_fit), "Raw Stan mixture parameters")

  expect_error(
    posterior_predictive_check(bayes_ppc_fit, type = "bad"),
    "should be one of"
  )
  expect_error(
    posterior_predictive_check(bayes_ppc_fit, ndraws = 0),
    "ndraws must be a single integer >= 1"
  )
  expect_error(
    posterior_predictive_check(bayes_ppc_fit, seed = -1),
    "seed must be NULL or a single non-negative integer"
  )
})

test_that("Bayesian PPC helper supports positive-support fits", {
  both_plots <- posterior_predictive_check(bayes_positive_ppc_fit, seed = 2)
  ppc_data <- TruncComp2:::bayes_ppc_data(
    bayes_positive_ppc_fit,
    ndraws = 6,
    seed = 21
  )

  expect_named(both_plots, c("atom", "continuous"))
  expect_s3_class(both_plots$continuous, "ggplot")
  expect_equal(both_plots$continuous$labels$y, "Density")
  expect_equal(both_plots$continuous$labels$x, "log(Outcome)")
  expect_equal(dim(ppc_data$yrep_cont), c(6, sum(bayes_positive_ppc_fit$data$A == 1)))
  expect_true(all(is.finite(ppc_data$yrep_cont)))
  expect_true(all(ppc_data$yrep_cont > 0))
})

test_that("Positive-support PPC plot inputs are formed on the log scale", {
  ppc_data <- TruncComp2:::bayes_ppc_data(
    bayes_positive_ppc_fit,
    ndraws = 6,
    seed = 25
  )
  plot_inputs <- TruncComp2:::bayes_ppc_continuous_plot_inputs(
    ppc_data,
    continuous_support = "positive_real"
  )

  expect_equal(plot_inputs$x_label, "log(Outcome)")
  expect_true(all(is.finite(plot_inputs$y)))
  expect_true(all(is.finite(plot_inputs$yrep)))
  expect_equal(plot_inputs$y, log(ppc_data$y_cont))
  expect_equal(plot_inputs$yrep, log(ppc_data$yrep_cont))
})
