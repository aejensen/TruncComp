fast_delta_example <- function() {
  data.frame(
    Y = c(
      0, 0, 1.2, 1.5, 1.8, 2.0, 2.2, 0, 1.7, 2.1,
      0, 2.0, 2.3, 2.5, 2.8, 3.0, 3.2, 0, 2.6, 2.9
    ),
    R = c(rep(0, 10), rep(1, 10))
  )
}

test_that("unadjusted fits expose Delta and compute Delta intervals on demand", {
  example_data <- fast_delta_example()

  fit_lrt <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = "lrt")
  fit_splrt <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = "splrt")
  empirical_delta <- mean(example_data$Y[example_data$R == 1]) -
    mean(example_data$Y[example_data$R == 0])
  welch_intervals <- list()

  for(fit in list(fit_lrt, fit_splrt)) {
    expect_true(fit$success)
    expect_equal(fit$delta, empirical_delta, tolerance = 1e-10)
    expect_false(any(c("DeltaCI", "DeltaMarginalCI", "DeltaProjectedCI", "DeltaProfileCI") %in% names(fit)))

    capture.output(delta_welch <- confint(fit, parameter = "delta", method = "welch"))
    delta_welch <- unname(delta_welch["delta (welch)", ])
    capture.output(delta_projected <- confint(fit, parameter = "delta", method = "projected", plot = FALSE))
    delta_projected <- unname(delta_projected["delta (projected)", ])
    capture.output(delta_profile <- confint(fit, parameter = "delta", method = "profile", plot = FALSE))
    delta_profile <- unname(delta_profile["delta (profile)", ])

    expect_true(all(is.finite(delta_welch)))
    expect_true(all(is.finite(delta_projected)))
    expect_true(all(is.finite(delta_profile)))
    expect_lte(delta_welch[1], fit$delta)
    expect_gte(delta_welch[2], fit$delta)
    expect_lte(delta_projected[1], fit$delta)
    expect_gte(delta_projected[2], fit$delta)
    expect_lte(delta_profile[1], fit$delta)
    expect_gte(delta_profile[2], fit$delta)

    welch_intervals[[length(welch_intervals) + 1]] <- delta_welch
  }

  expect_equal(welch_intervals[[1]], welch_intervals[[2]], tolerance = 1e-10)
})

test_that("grid-based DeltaProjectedCI agrees with high-resolution surface projection", {
  example_data <- fast_delta_example()

  for(method in c("lrt", "splrt")) {
    fit <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = method)
    surface_data <- TruncComp2:::delta_surface_for_inference(
      fit,
      conf.level = fit$conf.level,
      resolution = 41
    )
    grid_reference <- TruncComp2:::delta_interval_from_surface(
      surface_data,
      stats::qchisq(fit$conf.level, 2)
    )

    capture.output(projected <- confint(
      fit,
      parameter = "delta",
      method = "projected",
      algorithm = "grid",
      resolution = 41,
      plot = FALSE
    ))
    projected <- unname(projected["delta (projected)", ])
    expect_equal(projected, grid_reference, tolerance = 1e-10)
    expect_lte(projected[1], fit$delta)
    expect_gte(projected[2], fit$delta)
  }
})

test_that("grid-based DeltaProfileCI agrees with high-resolution surface profiling", {
  example_data <- fast_delta_example()

  for(method in c("lrt", "splrt")) {
    fit <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = method)
    surface_data <- TruncComp2:::delta_surface_for_inference(
      fit,
      conf.level = fit$conf.level,
      resolution = 41
    )
    grid_reference <- TruncComp2:::delta_interval_from_surface(
      surface_data,
      stats::qchisq(fit$conf.level, 1)
    )

    capture.output(profiled <- confint(
      fit,
      parameter = "delta",
      method = "profile",
      algorithm = "grid",
      resolution = 41,
      plot = FALSE
    ))
    profiled <- unname(profiled["delta (profile)", ])
    expect_equal(profiled, grid_reference, tolerance = 1e-10)
    expect_lte(profiled[1], fit$delta)
    expect_gte(profiled[2], fit$delta)
  }
})

test_that("optimizer-backed Delta profiling remains available for fast parametric smoke checks", {
  for(method in c("lrt", "splrt")) {
    fit <- trunc_comp(Y ~ R, atom = 0, data = fast_delta_example(), method = method)
    profile_fun <- TruncComp2:::delta_profile_factory(fit)
    expect_equal(profile_fun(fit$delta)$statistic, 0, tolerance = 1e-8)
  }

  fit <- trunc_comp(Y ~ R, atom = 0, data = fast_delta_example(), method = "lrt")
  capture.output(profiled <- confint(
    fit,
    parameter = "delta",
    method = "profile",
    algorithm = "optimize",
    plot = FALSE
  ))
  profiled <- unname(profiled["delta (profile)", ])
  expect_lte(profiled[1], fit$delta)
  expect_gte(profiled[2], fit$delta)
})

test_that("default grid-based Delta profile confint matches the direct helper", {
  example_data <- fast_delta_example()

  for(method in c("lrt", "splrt")) {
    fit <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = method)
    reference <- TruncComp2:::delta_profile_interval.grid(
      fit,
      conf.level = fit$conf.level,
      resolution = 35
    )
    capture.output(profiled <- confint(
      fit,
      parameter = "delta",
      method = "profile",
      algorithm = "grid",
      plot = FALSE
    ))
    profiled <- unname(profiled["delta (profile)", ])
    expect_equal(profiled, reference, tolerance = 1e-12)
  }
})

test_that("default interface stores and infers the atom value", {
  set.seed(20260415)
  simulated <- simulate_truncated_data(
    8,
    f0 = function(n) stats::rnorm(n, 2, 0.4),
    f1 = function(n) stats::rnorm(n, 3, 0.4),
    pi0 = 0.45,
    pi1 = 0.70,
    atom = -5
  )

  fit_formula <- trunc_comp(Y ~ R, atom = -5, data = simulated[, c("Y", "R")], method = "lrt")
  fit_default_explicit <- trunc_comp(simulated$Y, simulated$A, simulated$R, method = "lrt", atom = -5)
  fit_default_inferred <- trunc_comp(simulated$Y, simulated$A, simulated$R, method = "lrt")

  for(fit in list(fit_formula, fit_default_explicit, fit_default_inferred)) {
    expect_true(fit$success)
    expect_equal(fit$atom, -5)
    capture.output(delta_profile <- confint(fit, parameter = "delta", method = "profile", plot = FALSE))
    expect_true(all(is.finite(unname(delta_profile["delta (profile)", ]))))
  }

  expect_equal(fit_default_explicit$delta, fit_formula$delta, tolerance = 1e-10)
  expect_equal(fit_default_inferred$delta, fit_formula$delta, tolerance = 1e-10)

  ambiguous_y <- simulated$Y
  first_atom <- which(simulated$A == 0)[1]
  ambiguous_y[first_atom] <- -4
  expect_error(
    trunc_comp(ambiguous_y, simulated$A, simulated$R, method = "lrt"),
    "atom must be supplied"
  )
})
