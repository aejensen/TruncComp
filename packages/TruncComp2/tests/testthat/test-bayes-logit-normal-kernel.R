test_that("bounded_kernel dispatch and prior validation keep beta defaults", {
  continuous_data <- bayes_bounded_continuous_test_data()
  score_data <- bayes_bounded_score_test_data()
  real_data <- bayes_test_data()
  positive_data <- bayes_positive_test_data()

  expect_equal(
    TruncComp2:::bayes_model_name("bounded_continuous"),
    "trunc_comp_bayes_bounded_continuous"
  )
  expect_equal(
    TruncComp2:::bayes_model_name("bounded_score"),
    "trunc_comp_bayes_bounded_score"
  )
  expect_equal(
    TruncComp2:::bayes_model_name("bounded_continuous", bounded_kernel = "logit_normal"),
    "trunc_comp_bayes_bounded_continuous_logit_normal"
  )
  expect_equal(
    TruncComp2:::bayes_model_name("bounded_score", bounded_kernel = "logit_normal"),
    "trunc_comp_bayes_bounded_score_logit_normal"
  )

  expect_error(
    trunc_comp_bayes(
      Y ~ R,
      atom = 0,
      data = real_data,
      bounded_kernel = "beta"
    ),
    "bounded_kernel is only supported"
  )
  expect_error(
    trunc_comp_bayes(
      Y ~ R,
      atom = 0,
      data = positive_data,
      continuous_support = "positive_real",
      bounded_kernel = "logit_normal"
    ),
    "bounded_kernel is only supported"
  )
  expect_error(
    trunc_comp_bayes(
      Y ~ R,
      atom = -1,
      data = continuous_data,
      continuous_support = "bounded_continuous",
      score_min = 0,
      score_max = 100,
      bounded_kernel = "logit_normal",
      prior = list(m_alpha = 2)
    ),
    "Unsupported prior fields"
  )
  expect_error(
    trunc_comp_bayes(
      Y ~ R,
      atom = -1,
      data = continuous_data,
      continuous_support = "bounded_continuous",
      score_min = 0,
      score_max = 100,
      prior = list(mu_logit_mean = 0)
    ),
    "Unsupported prior fields"
  )

  logit_prior <- TruncComp2:::normalize_bayes_prior(
    NULL,
    continuous_support = "bounded_score",
    bounded_kernel = "logit_normal"
  )
  expect_true(all(c(
    "mu_logit_mean", "mu_logit_sd",
    "sigma_logit_meanlog", "sigma_logit_sdlog",
    "eta_prior"
  ) %in% names(logit_prior)))
  expect_false(any(c("m_alpha", "m_beta", "phi_meanlog", "phi_sdlog") %in% names(logit_prior)))

  beta_prior <- TruncComp2:::normalize_bayes_prior(
    NULL,
    continuous_support = "bounded_score"
  )
  expect_true(all(c("m_alpha", "m_beta", "phi_meanlog", "phi_sdlog", "eta_prior") %in% names(beta_prior)))
  expect_false(any(c("mu_logit_mean", "sigma_logit_meanlog") %in% names(beta_prior)))

  expect_true(all(score_data$Y[score_data$A == 1L] >= 0))
})

bayes_bounded_continuous_logit_fit <- bayes_bounded_continuous_formula_fit(
  seed = 2701,
  bounded_kernel = "logit_normal"
)
bayes_bounded_score_logit_fit <- bayes_bounded_score_formula_fit(
  seed = 2702,
  bounded_kernel = "logit_normal"
)

test_that("bounded logit-normal fits return successful objects and settings", {
  expect_s3_class(bayes_bounded_continuous_logit_fit, "trunc_comp_bayes_fit")
  expect_s3_class(bayes_bounded_score_logit_fit, "trunc_comp_bayes_fit")
  expect_true(bayes_bounded_continuous_logit_fit$success)
  expect_true(bayes_bounded_score_logit_fit$success)

  expect_equal(bayes_bounded_continuous_logit_fit$settings$continuous_support, "bounded_continuous")
  expect_equal(bayes_bounded_score_logit_fit$settings$continuous_support, "bounded_score")
  expect_equal(bayes_bounded_continuous_logit_fit$settings$bounded_kernel, "logit_normal")
  expect_equal(bayes_bounded_score_logit_fit$settings$bounded_kernel, "logit_normal")
  expect_equal(
    bayes_bounded_continuous_logit_fit$settings$model_name,
    "trunc_comp_bayes_bounded_continuous_logit_normal"
  )
  expect_equal(
    bayes_bounded_score_logit_fit$settings$model_name,
    "trunc_comp_bayes_bounded_score_logit_normal"
  )

  expect_true(all(TruncComp2:::bayes_parameter_names("all") %in% names(bayes_bounded_continuous_logit_fit$draws)))
  expect_true(all(TruncComp2:::bayes_parameter_names("all") %in% names(bayes_bounded_score_logit_fit$draws)))
  expect_true(all(is.finite(bayes_bounded_continuous_logit_fit$summary_table$estimate)))
  expect_true(all(is.finite(bayes_bounded_score_logit_fit$summary_table$estimate)))
  expect_true(any(bayes_bounded_score_logit_fit$data$Y[bayes_bounded_score_logit_fit$data$A == 1L] %in% c(0, 100)))
})

test_that("bounded-continuous logit-normal means agree with generated quantities", {
  fit <- bayes_bounded_continuous_logit_fit
  extracted <- posterior::as_draws_df(rstan::extract(
    fit$fit,
    pars = c("w", "mu_logit_comp", "sigma_logit_comp"),
    permuted = FALSE,
    inc_warmup = FALSE
  ))
  component_index <- seq_len(fit$settings$mixture_components)
  w_0 <- do.call(cbind, lapply(component_index, function(h) extracted[[sprintf("w[1,%d]", h)]]))
  w_1 <- do.call(cbind, lapply(component_index, function(h) extracted[[sprintf("w[2,%d]", h)]]))
  mu_0 <- do.call(cbind, lapply(component_index, function(h) extracted[[sprintf("mu_logit_comp[1,%d]", h)]]))
  mu_1 <- do.call(cbind, lapply(component_index, function(h) extracted[[sprintf("mu_logit_comp[2,%d]", h)]]))
  sigma_0 <- do.call(cbind, lapply(component_index, function(h) extracted[[sprintf("sigma_logit_comp[1,%d]", h)]]))
  sigma_1 <- do.call(cbind, lapply(component_index, function(h) extracted[[sprintf("sigma_logit_comp[2,%d]", h)]]))

  unit_mean_0 <- matrix(
    mapply(TruncComp2:::bayes_logitnormal_unit_mean, as.vector(mu_0), as.vector(sigma_0)),
    nrow = nrow(mu_0),
    ncol = ncol(mu_0)
  )
  unit_mean_1 <- matrix(
    mapply(TruncComp2:::bayes_logitnormal_unit_mean, as.vector(mu_1), as.vector(sigma_1)),
    nrow = nrow(mu_1),
    ncol = ncol(mu_1)
  )

  score_min <- fit$settings$score_min
  score_range <- fit$settings$score_range
  mu_0_from_components <- score_min + score_range * rowSums(w_0 * unit_mean_0)
  mu_1_from_components <- score_min + score_range * rowSums(w_1 * unit_mean_1)

  expect_equal(mu_0_from_components, fit$draws$mu_0_c, tolerance = 1e-10)
  expect_equal(mu_1_from_components, fit$draws$mu_1_c, tolerance = 1e-10)
})

test_that("bounded-score logit-normal means agree with generated score PMFs", {
  fit <- bayes_bounded_score_logit_fit
  extracted <- posterior::as_draws_df(rstan::extract(
    fit$fit,
    pars = "score_pmf",
    permuted = FALSE,
    inc_warmup = FALSE
  ))

  score_values <- fit$settings$score_values
  pmf_0 <- do.call(cbind, lapply(seq_along(score_values), function(j) {
    extracted[[sprintf("score_pmf[1,%d]", j)]]
  }))
  pmf_1 <- do.call(cbind, lapply(seq_along(score_values), function(j) {
    extracted[[sprintf("score_pmf[2,%d]", j)]]
  }))

  expect_equal(rowSums(pmf_0), rep(1, nrow(pmf_0)), tolerance = 1e-8)
  expect_equal(rowSums(pmf_1), rep(1, nrow(pmf_1)), tolerance = 1e-8)
  expect_equal(as.numeric(pmf_0 %*% score_values), fit$draws$mu_0_c, tolerance = 1e-10)
  expect_equal(as.numeric(pmf_1 %*% score_values), fit$draws$mu_1_c, tolerance = 1e-10)
})

test_that("PPC and posterior density helpers support bounded logit-normal kernels", {
  continuous_ppc <- TruncComp2:::bayes_ppc_data(
    bayes_bounded_continuous_logit_fit,
    ndraws = 6,
    seed = 51
  )
  score_ppc <- TruncComp2:::bayes_ppc_data(
    bayes_bounded_score_logit_fit,
    ndraws = 6,
    seed = 52
  )
  continuous_table <- posterior_predictive_pvalues(
    bayes_bounded_continuous_logit_fit,
    ndraws = 6,
    seed = 51
  )
  score_table <- posterior_predictive_pvalues(
    bayes_bounded_score_logit_fit,
    ndraws = 6,
    seed = 52
  )
  continuous_plot <- posterior_density_plot(bayes_bounded_continuous_logit_fit)
  score_plot <- posterior_density_plot(bayes_bounded_score_logit_fit)

  expect_true(all(continuous_ppc$yrep_cont >= 0 & continuous_ppc$yrep_cont <= 100))
  expect_true(all(score_ppc$yrep_cont %in% bayes_bounded_score_logit_fit$settings$score_values))
  expect_true(all(continuous_table$p_value >= 0 & continuous_table$p_value <= 1))
  expect_true(all(score_table$p_value >= 0 & score_table$p_value <= 1))
  expect_s3_class(continuous_plot, "ggplot")
  expect_s3_class(score_plot, "ggplot")

  continuous_plot_data <- TruncComp2:::bayes_density_plot_data(
    bayes_bounded_continuous_logit_fit,
    n = 30
  )
  score_plot_data <- TruncComp2:::bayes_density_plot_data(
    bayes_bounded_score_logit_fit,
    n = 30
  )
  expect_equal(continuous_plot_data$plot_type, "density")
  expect_equal(score_plot_data$plot_type, "mass")
  expect_true(all(is.finite(continuous_plot_data$density_data$density_mean)))
  expect_true(all(score_plot_data$density_data$density_mean >= 0))
})
