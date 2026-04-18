bayes_ppc_fit <- bayes_formula_fit(seed = 606)
bayes_positive_ppc_fit <- bayes_positive_formula_fit(seed = 1606)

test_that("posterior_predictive_check returns ggplots for the requested PPC type", {
  both_plots <- posterior_predictive_check(bayes_ppc_fit, seed = 1)
  atom_plot <- posterior_predictive_check(bayes_ppc_fit, type = "atom", seed = 1)
  continuous_plot <- posterior_predictive_check(bayes_ppc_fit, type = "continuous", seed = 1)
  ppc_table <- posterior_predictive_pvalues(bayes_ppc_fit, ndraws = 50, seed = 1)

  expect_type(both_plots, "list")
  expect_named(both_plots, c("atom", "continuous"))
  expect_s3_class(both_plots$atom, "ggplot")
  expect_s3_class(both_plots$continuous, "ggplot")
  expect_s3_class(atom_plot, "ggplot")
  expect_s3_class(continuous_plot, "ggplot")
  expect_equal(both_plots$continuous$labels$y, "Density")
  expect_equal(both_plots$continuous$labels$x, "Outcome")
  expect_match(
    both_plots$atom$labels$subtitle,
    sprintf("%.3f", ppc_table["atom", "p_value"])
  )
  expect_match(
    both_plots$continuous$labels$subtitle,
    sprintf("%.3f", ppc_table["continuous", "p_value"])
  )
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

test_that("posterior_predictive_pvalues returns structured and reproducible summaries", {
  ppc_first <- posterior_predictive_pvalues(
    bayes_ppc_fit,
    ndraws = 8,
    seed = 99
  )
  ppc_second <- posterior_predictive_pvalues(
    bayes_ppc_fit,
    ndraws = 8,
    seed = 99
  )
  clipped <- posterior_predictive_pvalues(
    bayes_ppc_fit,
    ndraws = nrow(bayes_ppc_fit$draws) + 100,
    seed = 123
  )

  expect_s3_class(ppc_first, "data.frame")
  expect_equal(rownames(ppc_first), c("atom", "continuous"))
  expect_true(all(c(
    "p_value", "statistic", "scale", "ndraws",
    "mean_observed_discrepancy", "mean_replicated_discrepancy"
  ) %in% names(ppc_first)))
  expect_true(all(ppc_first$p_value >= 0 & ppc_first$p_value <= 1))
  expect_equal(ppc_first, ppc_second)
  expect_equal(clipped$ndraws, c(nrow(bayes_ppc_fit$draws), nrow(bayes_ppc_fit$draws)))
})

test_that("posterior_predictive_pvalues validates failed fits and invalid inputs", {
  failed_fit <- bayes_ppc_fit
  failed_fit$success <- FALSE
  expect_error(posterior_predictive_pvalues(failed_fit), "Estimation failed")

  missing_fit <- bayes_ppc_fit
  missing_fit$fit <- NULL
  expect_error(posterior_predictive_pvalues(missing_fit), "Raw Stan mixture parameters")

  expect_error(
    posterior_predictive_pvalues(bayes_ppc_fit, ndraws = 0),
    "ndraws must be a single integer >= 1"
  )
  expect_error(
    posterior_predictive_pvalues(bayes_ppc_fit, seed = -1),
    "seed must be NULL or a single non-negative integer"
  )
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

test_that("Bayesian PPC discrepancy summaries match the implemented formulas", {
  ppc_summary <- TruncComp2:::bayes_ppc_summary(
    bayes_ppc_fit,
    ndraws = 6,
    seed = 11
  )
  ppc_data <- TruncComp2:::bayes_ppc_data(
    bayes_ppc_fit,
    ndraws = 6,
    seed = 11
  )
  parameters <- TruncComp2:::bayes_ppc_component_parameters(bayes_ppc_fit)

  atom_obs <- atom_rep <- cont_obs <- cont_rep <- numeric(length(ppc_data$draw_indices))
  for(s in seq_along(ppc_data$draw_indices)) {
    draw <- ppc_data$draw_indices[s]
    atom_obs[s] <- TruncComp2:::bayes_ppc_atom_discrepancy(
      y_atom = ppc_data$y_atom,
      arm_full = ppc_data$arm_full,
      rho_draw = parameters$rho[draw, ]
    )
    atom_rep[s] <- TruncComp2:::bayes_ppc_atom_discrepancy(
      y_atom = ppc_data$yrep_atom[s, ],
      arm_full = ppc_data$arm_full,
      rho_draw = parameters$rho[draw, ]
    )
    cont_obs[s] <- TruncComp2:::bayes_ppc_continuous_discrepancy(
      y_cont = ppc_data$y_cont,
      arm_obs = ppc_data$arm_obs,
      weights_draw = parameters$weights[draw, , ],
      means_draw = parameters$means[draw, , ],
      continuous_support = parameters$support,
      sds_draw = parameters$sds[draw, , ]
    )
    cont_rep[s] <- TruncComp2:::bayes_ppc_continuous_discrepancy(
      y_cont = ppc_data$yrep_cont[s, ],
      arm_obs = ppc_data$arm_obs,
      weights_draw = parameters$weights[draw, , ],
      means_draw = parameters$means[draw, , ],
      continuous_support = parameters$support,
      sds_draw = parameters$sds[draw, , ]
    )
  }

  expect_equal(ppc_summary$details$atom_observed_discrepancy, atom_obs)
  expect_equal(ppc_summary$details$atom_replicated_discrepancy, atom_rep)
  expect_equal(ppc_summary$details$continuous_observed_discrepancy, cont_obs)
  expect_equal(ppc_summary$details$continuous_replicated_discrepancy, cont_rep)
  expect_equal(ppc_summary$table["atom", "p_value"], mean(atom_rep >= atom_obs))
  expect_equal(ppc_summary$table["continuous", "p_value"], mean(cont_rep >= cont_obs))
  expect_equal(ppc_summary$table["continuous", "scale"], "Outcome")
})

test_that("Positive-support PPC p-values use the log-scale continuous discrepancy", {
  ppc_summary <- TruncComp2:::bayes_ppc_summary(
    bayes_positive_ppc_fit,
    ndraws = 6,
    seed = 25
  )
  ppc_data <- TruncComp2:::bayes_ppc_data(
    bayes_positive_ppc_fit,
    ndraws = 6,
    seed = 25
  )
  parameters <- TruncComp2:::bayes_ppc_component_parameters(bayes_positive_ppc_fit)

  manual_discrepancy <- TruncComp2:::bayes_ppc_continuous_discrepancy(
    y_cont = ppc_data$y_cont,
    arm_obs = ppc_data$arm_obs,
    weights_draw = parameters$weights[ppc_data$draw_indices[1], , ],
    means_draw = parameters$means[ppc_data$draw_indices[1], , ],
    continuous_support = parameters$support,
    shapes_draw = parameters$shapes[ppc_data$draw_indices[1], , ]
  )

  expect_true(is.finite(manual_discrepancy))
  expect_equal(ppc_summary$table["continuous", "scale"], "log(Outcome)")
})
