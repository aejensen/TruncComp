fast_delta_example <- function() {
  data.frame(
    Y = c(
      0, 0, 1.2, 1.5, 1.8, 2.0, 2.2, 0, 1.7, 2.1,
      0, 2.0, 2.3, 2.5, 2.8, 3.0, 3.2, 0, 2.6, 2.9
    ),
    R = c(rep(0, 10), rep(1, 10))
  )
}

test_that("unadjusted fits expose Delta and the three Delta intervals", {
  example_data <- fast_delta_example()

  fit_lrt <- truncComp(Y ~ R, atom = 0, data = example_data, method = "LRT")
  fit_splrt <- truncComp(Y ~ R, atom = 0, data = example_data, method = "SPLRT")
  empirical_delta <- mean(example_data$Y[example_data$R == 1]) -
    mean(example_data$Y[example_data$R == 0])

  for(fit in list(fit_lrt, fit_splrt)) {
    expect_true(fit$success)
    expect_equal(fit$Delta, empirical_delta, tolerance = 1e-10)
    expect_equal(fit$DeltaCI, fit$DeltaProfileCI, tolerance = 1e-12)
    expect_true(all(is.finite(fit$DeltaMarginalCI)))
    expect_true(all(is.na(fit$DeltaProjectedCI)))
    expect_true(all(is.finite(fit$DeltaProfileCI)))
    expect_lte(fit$DeltaMarginalCI[1], fit$Delta)
    expect_gte(fit$DeltaMarginalCI[2], fit$Delta)
    expect_lte(fit$DeltaProfileCI[1], fit$Delta)
    expect_gte(fit$DeltaProfileCI[2], fit$Delta)
  }

  expect_equal(fit_lrt$DeltaMarginalCI, fit_splrt$DeltaMarginalCI, tolerance = 1e-10)
})

test_that("grid-based DeltaProjectedCI agrees with high-resolution surface projection", {
  example_data <- fast_delta_example()

  for(method in c("LRT", "SPLRT")) {
    fit <- truncComp(Y ~ R, atom = 0, data = example_data, method = method)
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
      type = "delta_projected",
      algorithm = "grid",
      resolution = 41,
      plot = FALSE
    ))
    projected <- unname(projected["Delta (projected)", ])
    expect_equal(projected, grid_reference, tolerance = 1e-10)
    expect_lte(projected[1], fit$Delta)
    expect_gte(projected[2], fit$Delta)
  }
})

test_that("grid-based DeltaProfileCI agrees with high-resolution surface profiling", {
  example_data <- fast_delta_example()

  for(method in c("LRT", "SPLRT")) {
    fit <- truncComp(Y ~ R, atom = 0, data = example_data, method = method)
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
      type = "delta_profile",
      algorithm = "grid",
      resolution = 41,
      plot = FALSE
    ))
    profiled <- unname(profiled["Delta (profile likelihood)", ])
    expect_equal(profiled, grid_reference, tolerance = 1e-10)
    expect_lte(profiled[1], fit$Delta)
    expect_gte(profiled[2], fit$Delta)
  }
})

test_that("optimizer-backed Delta profiling remains available for fast parametric smoke checks", {
  for(method in c("LRT", "SPLRT")) {
    fit <- truncComp(Y ~ R, atom = 0, data = fast_delta_example(), method = method)
    profile_fun <- TruncComp2:::delta_profile_factory(fit)
    expect_equal(profile_fun(fit$Delta)$statistic, 0, tolerance = 1e-8)
  }

  fit <- truncComp(Y ~ R, atom = 0, data = fast_delta_example(), method = "LRT")
  capture.output(profiled <- confint(
    fit,
    type = "delta_profile",
    algorithm = "optimize",
    plot = FALSE
  ))
  profiled <- unname(profiled["Delta (profile likelihood)", ])
  expect_lte(profiled[1], fit$Delta)
  expect_gte(profiled[2], fit$Delta)
})

test_that("stored DeltaProfileCI matches the default grid-based confint result", {
  example_data <- fast_delta_example()

  for(method in c("LRT", "SPLRT")) {
    fit <- truncComp(Y ~ R, atom = 0, data = example_data, method = method)
    capture.output(profiled <- confint(
      fit,
      type = "delta_profile",
      algorithm = "grid",
      plot = FALSE
    ))
    profiled <- unname(profiled["Delta (profile likelihood)", ])
    expect_equal(profiled, fit$DeltaProfileCI, tolerance = 1e-12)
  }
})

test_that("default interface stores and infers the atom value", {
  set.seed(20260415)
  simulated <- simulateTruncatedData(
    8,
    f0 = function(n) stats::rnorm(n, 2, 0.4),
    f1 = function(n) stats::rnorm(n, 3, 0.4),
    pi0 = 0.45,
    pi1 = 0.70,
    atom = -5
  )

  fit_formula <- truncComp(Y ~ R, atom = -5, data = simulated[, c("Y", "R")], method = "LRT")
  fit_default_explicit <- truncComp.default(simulated$Y, simulated$A, simulated$R, method = "LRT", atom = -5)
  fit_default_inferred <- truncComp.default(simulated$Y, simulated$A, simulated$R, method = "LRT")

  for(fit in list(fit_formula, fit_default_explicit, fit_default_inferred)) {
    expect_true(fit$success)
    expect_equal(fit$atom, -5)
    expect_true(all(is.finite(fit$DeltaProfileCI)))
  }

  expect_equal(fit_default_explicit$Delta, fit_formula$Delta, tolerance = 1e-10)
  expect_equal(fit_default_inferred$Delta, fit_formula$Delta, tolerance = 1e-10)

  ambiguous_y <- simulated$Y
  first_atom <- which(simulated$A == 0)[1]
  ambiguous_y[first_atom] <- -4
  expect_error(
    truncComp.default(ambiguous_y, simulated$A, simulated$R, method = "LRT"),
    "atom must be supplied"
  )
})
