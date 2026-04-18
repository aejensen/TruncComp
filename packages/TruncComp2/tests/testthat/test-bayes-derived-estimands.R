test_that("Derived Bayesian estimands agree with the stored posterior draws", {
  fit <- bayes_formula_fit(seed = 303)
  draws <- fit$draws

  alpha_from_pi <- (draws$pi_1 / (1 - draws$pi_1)) / (draws$pi_0 / (1 - draws$pi_0))
  delta_atom_from_rho <- draws$rho_1 - draws$rho_0
  mu_delta_from_means <- draws$mu_1_c - draws$mu_0_c
  delta_from_parts <- (fit$atom * draws$rho_1 + draws$pi_1 * draws$mu_1_c) -
    (fit$atom * draws$rho_0 + draws$pi_0 * draws$mu_0_c)

  expect_equal(alpha_from_pi, draws$alpha_delta, tolerance = 1e-10)
  expect_equal(delta_atom_from_rho, draws$delta_atom, tolerance = 1e-10)
  expect_equal(mu_delta_from_means, draws$mu_delta, tolerance = 1e-10)
  expect_equal(delta_from_parts, draws$delta, tolerance = 1e-10)

  capture.output(intervals <- confint(fit))
  estimates <- coef(fit)

  expect_lte(intervals["mu_delta", 1], estimates["mu_delta"])
  expect_gte(intervals["mu_delta", 2], estimates["mu_delta"])
  expect_lte(intervals["delta_atom", 1], estimates["delta_atom"])
  expect_gte(intervals["delta_atom", 2], estimates["delta_atom"])
  expect_lte(intervals["alpha_delta", 1], estimates["alpha_delta"])
  expect_gte(intervals["alpha_delta", 2], estimates["alpha_delta"])
  expect_lte(intervals["delta", 1], estimates["delta"])
  expect_gte(intervals["delta", 2], estimates["delta"])
})

test_that("Positive-support arm means agree with the Gamma-mixture draw-wise formula", {
  fit <- bayes_positive_formula_fit(seed = 1302)
  extracted <- posterior::as_draws_df(rstan::extract(
    fit$fit,
    pars = c("w", "mean_comp"),
    permuted = FALSE,
    inc_warmup = FALSE
  ))
  component_index <- seq_len(fit$settings$mixture_components)
  w_0 <- do.call(cbind, lapply(component_index, function(h) extracted[[sprintf("w[1,%d]", h)]]))
  w_1 <- do.call(cbind, lapply(component_index, function(h) extracted[[sprintf("w[2,%d]", h)]]))
  mean_0 <- do.call(cbind, lapply(component_index, function(h) extracted[[sprintf("mean_comp[1,%d]", h)]]))
  mean_1 <- do.call(cbind, lapply(component_index, function(h) extracted[[sprintf("mean_comp[2,%d]", h)]]))

  mu_0_from_components <- fit$settings$y_scale * rowSums(w_0 * mean_0)
  mu_1_from_components <- fit$settings$y_scale * rowSums(w_1 * mean_1)
  delta_from_parts <- (fit$atom * fit$draws$rho_1 + fit$draws$pi_1 * fit$draws$mu_1_c) -
    (fit$atom * fit$draws$rho_0 + fit$draws$pi_0 * fit$draws$mu_0_c)

  expect_equal(mu_0_from_components, fit$draws$mu_0_c, tolerance = 1e-10)
  expect_equal(mu_1_from_components, fit$draws$mu_1_c, tolerance = 1e-10)
  expect_equal(delta_from_parts, fit$draws$delta, tolerance = 1e-10)
})
