test_that("Bayesian mixture-component ladder doubles up to the requested maximum", {
  expect_equal(
    TruncComp2:::bayes_mixture_component_ladder(10, 40),
    c(10L, 20L, 40L)
  )
  expect_equal(
    TruncComp2:::bayes_mixture_component_ladder(6, 20),
    c(6L, 12L)
  )
  expect_equal(
    TruncComp2:::bayes_mixture_component_ladder(8, 8),
    8L
  )
})

test_that("Bayesian truncation thresholds use min(0.02, 1 / n_obs)", {
  data <- data.frame(
    Y = c(rep(1, 100), rep(0, 4), rep(1, 3), rep(0, 2)),
    A = c(rep(1L, 100), rep(0L, 4), rep(1L, 3), rep(0L, 2)),
    R = c(rep(0L, 104), rep(1L, 5))
  )

  thresholds <- TruncComp2:::bayes_truncation_thresholds(data)

  expect_equal(thresholds$n_obs, c(n_obs_0 = 100L, n_obs_1 = 3L))
  expect_equal(thresholds$thresholds, c(threshold_0 = 0.01, threshold_1 = 0.02))
})

test_that("Omitted-tail draws are computed from alpha and the final retained stick weight", {
  raw_draws <- data.frame(
    "alpha[1]" = c(1, 3),
    "alpha[2]" = c(2, 4),
    "w[1,10]" = c(0.1, 0.2),
    "w[2,10]" = c(0.3, 0.4),
    .chain = c(1, 1),
    .iteration = c(1, 2),
    .draw = c(1, 2),
    check.names = FALSE
  )

  truncation_draws <- TruncComp2:::bayes_compute_omitted_tail_draws(
    draws = raw_draws,
    mixture_components = 10
  )

  expect_equal(truncation_draws$omitted_tail_mass_0, c(0.05, 0.15))
  expect_equal(truncation_draws$omitted_tail_mass_1, c(0.2, 0.32))
  expect_equal(truncation_draws$.draw, raw_draws$.draw)
})

test_that("Truncation acceptance requires small omitted tails and clean diagnostics", {
  expect_true(TruncComp2:::bayes_truncation_ok(
    q95 = c(0.01, 0.015),
    thresholds = c(0.02, 0.02),
    convergence_ok = TRUE,
    divergences = 0L
  ))

  expect_false(TruncComp2:::bayes_truncation_ok(
    q95 = c(0.03, 0.015),
    thresholds = c(0.02, 0.02),
    convergence_ok = TRUE,
    divergences = 0L
  ))

  expect_false(TruncComp2:::bayes_truncation_ok(
    q95 = c(0.01, 0.015),
    thresholds = c(0.02, 0.02),
    convergence_ok = FALSE,
    divergences = 0L
  ))

  expect_false(TruncComp2:::bayes_truncation_ok(
    q95 = c(0.01, 0.015),
    thresholds = c(0.02, 0.02),
    convergence_ok = TRUE,
    divergences = 1L
  ))
})

test_that("Auto-selection stores the attempted truncation path and final retained H", {
  fake_fit <- function(mixture_components, settings, diagnostic_ok) {
    truncation_ok <- diagnostic_ok
    diagnostics <- list(
      divergences = 0L,
      max_rhat = 1,
      min_bulk_ess = 500,
      min_tail_ess = 500,
      parameter_table = data.frame(
        variable = c("delta_atom", "mu_delta", "alpha_delta", "delta"),
        rhat = rep(1, 4),
        ess_bulk = rep(500, 4),
        ess_tail = rep(500, 4),
        row.names = c("delta_atom", "mu_delta", "alpha_delta", "delta")
      ),
      core_parameter_table = data.frame(
        variable = c("delta_atom", "mu_delta", "alpha_delta", "delta"),
        rhat = rep(1, 4),
        ess_bulk = rep(500, 4),
        ess_tail = rep(500, 4),
        row.names = c("delta_atom", "mu_delta", "alpha_delta", "delta")
      ),
      core_ok = TRUE,
      truncation = list(
        n_obs = c(n_obs_0 = 8L, n_obs_1 = 9L),
        thresholds = c(threshold_0 = 0.02, threshold_1 = 0.02),
        q95 = c(
          q95_tail_0 = if(diagnostic_ok) 0.01 else 0.03,
          q95_tail_1 = if(diagnostic_ok) 0.01 else 0.03
        ),
        parameter_table = data.frame(
          variable = c("omitted_tail_mass_0", "omitted_tail_mass_1"),
          rhat = c(1, 1),
          ess_bulk = c(500, 500),
          ess_tail = c(500, 500),
          row.names = c("omitted_tail_mass_0", "omitted_tail_mass_1")
        ),
        max_rhat = 1,
        min_bulk_ess = 500,
        min_tail_ess = 500,
        convergence_ok = TRUE,
        truncation_ok = truncation_ok
      ),
      truncation_ok = truncation_ok,
      diagnostic_ok = diagnostic_ok
    )

    TruncComp2:::new_trunc_comp_bayes_fit(
      draws = data.frame(delta = 0),
      summary_table = data.frame(
        estimate = 0,
        conf.low = -1,
        conf.high = 1,
        posterior_prob = 0.5,
        row.names = "delta"
      ),
      arm_table = data.frame(
        estimate = 0,
        conf.low = -1,
        conf.high = 1,
        row.names = "rho_0"
      ),
      diagnostics = diagnostics,
      settings = utils::modifyList(settings, list(mixture_components = mixture_components)),
      conf.level = 0.95,
      success = TRUE,
      data = data.frame(Y = c(0, 1, 0, 2), A = c(0L, 1L, 0L, 1L), R = c(0L, 0L, 1L, 1L)),
      atom = 0
    )
  }

  data <- bayes_test_data()

  testthat::local_mocked_bindings(
    bayes_package_stanmodel = function(model_name) structure(list(name = model_name), class = "stanmodel"),
    validate_bayes_support_data = function(data, continuous_support) invisible(TRUE),
    fit_trunc_comp_bayes_once = function(data, atom, conf.level, continuous_support,
                                         mixture_components, chains, iter_warmup,
                                         iter_sampling, seed, refresh, control,
                                         prior, call, extra_args, model_object, settings) {
      fake_fit(
        mixture_components = mixture_components,
        settings = settings,
        diagnostic_ok = mixture_components >= 20
      )
    },
    bayes_finalize_trunc_comp_bayes_fit = function(fit_object) fit_object,
    .package = "TruncComp2"
  )

  fit_auto <- TruncComp2:::fit_trunc_comp_bayes(
    data = data,
    atom = 0,
    mixture_components = 10,
    auto_select_mixture_components = TRUE,
    mixture_components_max = 40,
    chains = 2,
    iter_warmup = 10,
    iter_sampling = 10
  )

  expect_true(fit_auto$success)
  expect_true(fit_auto$settings$auto_select_mixture_components)
  expect_equal(fit_auto$settings$mixture_components_initial, 10L)
  expect_equal(fit_auto$settings$mixture_components_final, 20L)
  expect_equal(fit_auto$settings$mixture_component_path, c(10L, 20L))
  expect_equal(
    fit_auto$settings$mixture_selection_history$accepted,
    c(FALSE, TRUE)
  )

  fit_fixed <- TruncComp2:::fit_trunc_comp_bayes(
    data = data,
    atom = 0,
    mixture_components = 10,
    auto_select_mixture_components = FALSE,
    chains = 2,
    iter_warmup = 10,
    iter_sampling = 10
  )

  expect_true(fit_fixed$success)
  expect_false(fit_fixed$settings$auto_select_mixture_components)
  expect_equal(fit_fixed$settings$mixture_components_initial, 10L)
  expect_equal(fit_fixed$settings$mixture_components_final, 10L)
  expect_equal(fit_fixed$settings$mixture_components_max, 10L)
  expect_equal(fit_fixed$settings$mixture_component_path, 10L)
  expect_equal(nrow(fit_fixed$settings$mixture_selection_history), 1L)
})
