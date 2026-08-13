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

test_that("full coded likelihood criterion agrees with the component profile surface", {
  example_data <- fast_delta_example()

  for(method in c("lrt", "splrt")) {
    fit <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = method)
    reference <- TruncComp:::delta_full_reference(fit)
    expect_false(is.null(reference))

    center <- TruncComp:::delta_full_candidate(reference, reference$center)
    expect_true(center$feasible)
    expect_equal(center$statistic, 0, tolerance = 1e-9)
    expect_equal(center$Delta, fit$delta, tolerance = 1e-10)

    mu_delta <- fit$mu_delta + 0.2
    log_or_delta <- log(fit$alpha_delta) - 0.15
    if(identical(method, "lrt")) {
      profiled <- TruncComp:::parametricJointCandidate(
        TruncComp:::prepareParametricJointReference(fit$data),
        mu_delta,
        log_or_delta
      )
    } else {
      profiled <- TruncComp:::splrtJointCandidate(
        TruncComp:::prepareSPLRTJointReference(fit$data),
        mu_delta,
        log_or_delta
      )
    }
    full <- TruncComp:::delta_full_candidate(
      reference,
      c(profiled$mu0, profiled$mu1, stats::qlogis(profiled$p0), stats::qlogis(profiled$p1))
    )
    expect_equal(full$statistic, profiled$statistic, tolerance = 1e-7)
  }
})

test_that("finite sublevel bounds and endpoint solutions attain their cutoffs", {
  example_data <- fast_delta_example()

  for(method in c("lrt", "splrt")) {
    fit <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = method)
    reference <- TruncComp:::delta_full_reference(fit)
    threshold <- stats::qchisq(fit$conf.level, 1)

    for(arm in 0:1) {
      bounds <- TruncComp:::delta_logit_sublevel_bounds(reference, arm, threshold)
      expect_true(all(is.finite(bounds)))
      expect_equal(
        vapply(bounds, function(eta) {
          TruncComp:::delta_bernoulli_arm_deviance(reference, arm, eta)
        }, numeric(1)),
        rep(threshold, 2),
        tolerance = 1e-7
      )
    }

    if(identical(method, "lrt")) {
      for(direction in c("lower", "upper")) {
        solution <- TruncComp:::delta_parametric_sublevel_endpoint(
          reference,
          threshold,
          direction = direction
        )
        expect_false(is.null(solution))
        expect_true(solution$candidate$feasible)
        expect_equal(solution$candidate$statistic, threshold, tolerance = 1e-6)
        expect_equal(
          solution$candidate$Delta,
          TruncComp:::delta_from_components(
            fit$atom,
            solution$candidate$mu0,
            solution$candidate$mu1,
            solution$candidate$p0,
            solution$candidate$p1
          ),
          tolerance = 1e-10
        )
      }
    } else {
      interval <- TruncComp:::delta_splrt_sublevel_interval(reference, threshold)
      coded0 <- c(rep(fit$atom, reference$n0 - reference$k0), reference$y0)
      coded1 <- c(rep(fit$atom, reference$n1 - reference$k1), reference$y1)
      expect_equal(
        vapply(interval, function(delta) {
          TruncComp:::el_mean_diff_statistic(coded1, coded0, delta)
        }, numeric(1)),
        rep(threshold, 2),
        tolerance = 1e-5
      )
    }
  }
})

test_that("full projected range contains accepted four parameter candidates", {
  example_data <- fast_delta_example()

  for(method in c("lrt", "splrt")) {
    fit <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = method)
    reference <- TruncComp:::delta_full_reference(fit)
    threshold <- stats::qchisq(fit$conf.level, 4)
    candidate <- TruncComp:::delta_full_candidate(
      reference,
      reference$center + c(0.05, -0.05, 0.05, -0.05)
    )
    expect_true(candidate$feasible)
    expect_lte(candidate$statistic, threshold)

    projected <- TruncComp:::delta_projected_interval.optimize(fit)
    expect_lte(projected[1], candidate$Delta + 1e-6)
    expect_gte(projected[2], candidate$Delta - 1e-6)
    expect_lte(projected[1], fit$delta)
    expect_gte(projected[2], fit$delta)
  }
})

test_that("profile ranges are nested in projected ranges", {
  for(method in c("lrt", "splrt")) {
    fit <- trunc_comp(Y ~ R, atom = 0, data = fast_delta_example(), method = method)
    profiled <- TruncComp:::delta_profile_interval(fit)
    projected <- TruncComp:::delta_projected_interval(fit)
    expect_true(all(is.finite(c(profiled, projected))))
    expect_gte(profiled[1], projected[1] - 1e-6)
    expect_lte(profiled[2], projected[2] + 1e-6)
    expect_lte(profiled[1], fit$delta)
    expect_gte(profiled[2], fit$delta)
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
  placeholder_y <- simulated$Y
  placeholder_y[simulated$A == 0] <- seq_len(sum(simulated$A == 0)) + 100
  fit_default_placeholder <- trunc_comp(
    placeholder_y,
    simulated$A,
    simulated$R,
    method = "lrt",
    atom = -5
  )

  for(fit in list(fit_formula, fit_default_explicit, fit_default_inferred,
                  fit_default_placeholder)) {
    expect_true(fit$success)
    expect_equal(fit$atom, -5)
    capture.output(delta_profile <- confint(fit, parameter = "delta", method = "profile", plot = FALSE))
    expect_true(all(is.finite(unname(delta_profile["delta (profile)", ]))))
  }

  expect_equal(fit_default_explicit$delta, fit_formula$delta, tolerance = 1e-10)
  expect_equal(fit_default_inferred$delta, fit_formula$delta, tolerance = 1e-10)
  expect_equal(fit_default_placeholder$delta, fit_formula$delta, tolerance = 1e-10)
  expect_equal(
    TruncComp:::delta_welch_interval(
      placeholder_y,
      simulated$A,
      simulated$R,
      -5,
      0.95
    ),
    TruncComp:::delta_welch_interval(
      simulated$Y,
      simulated$A,
      simulated$R,
      -5,
      0.95
    ),
    tolerance = 1e-12
  )

  ambiguous_y <- simulated$Y
  first_atom <- which(simulated$A == 0)[1]
  ambiguous_y[first_atom] <- -4
  expect_error(
    trunc_comp(ambiguous_y, simulated$A, simulated$R, method = "lrt"),
    "atom must be supplied"
  )
})
