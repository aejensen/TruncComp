test_that("el_one_sample_lambda handles degenerate and infeasible inputs", {
  zero_fit <- TruncComp2:::el_one_sample_lambda(c(0, 0, 0))
  expect_true(zero_fit$feasible)
  expect_equal(zero_fit$lambda, 0)
  expect_equal(zero_fit$statistic, 0)
  expect_equal(zero_fit$gradient, 0)

  expect_false(TruncComp2:::el_one_sample_lambda(c(1, 2, 3))$feasible)
  expect_false(TruncComp2:::el_one_sample_lambda(c(-3, -2, 0))$feasible)

  residuals <- c(-1.5, -0.25, 0.5, 2)
  fit <- TruncComp2:::el_one_sample_lambda(residuals)
  root_value <- sum(residuals / (1 + fit$lambda * residuals))
  expect_true(fit$feasible)
  expect_lt(abs(root_value), 1e-8)
  expect_gte(fit$statistic, 0)
})

test_that("el_mean_diff_statistic respects feasibility and point-mass cases", {
  x <- c(1, 2, 3)
  y <- c(2, 3, 4, 5)
  estimate <- mean(x) - mean(y)
  delta_range <- TruncComp2:::el_mean_diff_delta_range(x, y)

  expect_equal(TruncComp2:::el_mean_diff_statistic(x, y, estimate), 0, tolerance = 1e-12)
  expect_gte(TruncComp2:::el_mean_diff_statistic(x, y, estimate + 0.25), 0)
  expect_true(is.infinite(TruncComp2:::el_mean_diff_statistic(x, y, delta_range[1] - 0.1)))
  expect_true(is.infinite(TruncComp2:::el_mean_diff_statistic(x, y, delta_range[2] + 0.1)))

  same_fit <- TruncComp2:::el_mean_diff_fit(c(1, 1, 1), c(1, 1, 1))
  expect_equal(as.numeric(same_fit$estimate), 0)
  expect_equal(as.numeric(same_fit$statistic), 0)
  expect_equal(as.numeric(same_fit$conf.int), c(0, 0))

  separated_fit <- TruncComp2:::el_mean_diff_fit(c(1, 1, 1), c(2, 2, 2))
  expect_equal(as.numeric(separated_fit$estimate), -1)
  expect_true(is.infinite(as.numeric(separated_fit$statistic)))
  expect_equal(as.numeric(separated_fit$conf.int), c(-1, -1))
})

test_that("custom engine matches frozen upstream EL fixtures", {
  fixture <- readRDS(fixture_path("el_means_reference.rds"))

  expect_equal(fixture$generated_on, as.Date("2026-04-15"))
  expect_equal(fixture$upstream_package, "EL")

  for(case in fixture$cases) {
    fit <- TruncComp2:::el_mean_diff_fit(case$x, case$y)
    ref <- case$reference

    expect_true(abs(as.numeric(fit$estimate) - ref$estimate) < 1e-4,
                info = case$name)
    expect_true(all(abs(as.numeric(fit$conf.int) - ref$conf.int) < 5e-4),
                info = case$name)

    if(is.infinite(ref$statistic)) {
      expect_true(is.infinite(as.numeric(fit$statistic)), info = case$name)
    } else {
      expect_true(abs(as.numeric(fit$statistic) - ref$statistic) < 5e-4,
                  info = case$name)
    }

    expect_true(abs(as.numeric(fit$p.value) - ref$p.value) < 5e-4,
                info = case$name)
  }
})

test_that("empirical-likelihood engine satisfies symmetry and distributional properties", {
  set.seed(20260415)

  for(i in 1:20) {
    x <- stats::rnorm(sample(3:8, 1), mean = stats::runif(1, -0.5, 0.5))
    y <- stats::rnorm(sample(3:8, 1), mean = stats::runif(1, -0.5, 0.5))

    fit_xy <- TruncComp2:::el_mean_diff_fit(x, y)
    fit_yx <- TruncComp2:::el_mean_diff_fit(y, x)

    expect_equal(as.numeric(fit_xy$estimate), -as.numeric(fit_yx$estimate), tolerance = 1e-12)

    delta_range <- TruncComp2:::el_mean_diff_delta_range(x, y)
    margin <- 0.1 * (delta_range[2] - delta_range[1])
    delta <- stats::runif(1, delta_range[1] + margin, delta_range[2] - margin)

    stat_xy <- TruncComp2:::el_mean_diff_statistic(x, y, delta)
    stat_yx <- TruncComp2:::el_mean_diff_statistic(y, x, -delta)
    expect_equal(stat_xy, stat_yx, tolerance = 1e-8)

    ci <- as.numeric(fit_xy$conf.int)
    expect_lte(ci[1], as.numeric(fit_xy$estimate))
    expect_gte(ci[2], as.numeric(fit_xy$estimate))
    expect_gte(as.numeric(fit_xy$p.value), 0)
    expect_lte(as.numeric(fit_xy$p.value), 1)

    if(is.finite(as.numeric(fit_xy$statistic))) {
      expect_equal(as.numeric(fit_xy$statistic),
                   stats::qchisq(1 - as.numeric(fit_xy$p.value), 1),
                   tolerance = 1e-8)
    }
  }
})

test_that("empirical-likelihood helpers reject invalid inputs", {
  expect_error(TruncComp2:::el_mean_diff_fit("a", 1:3), "numeric")
  expect_error(TruncComp2:::el_mean_diff_fit(c(1, NA), 1:3), "finite")
  expect_error(TruncComp2:::el_mean_diff_fit(1, 1:3), "at least two")
  expect_error(TruncComp2:::el_mean_diff_fit(1:3, 1:3, conf.level = 1), "strictly between 0 and 1")
})

test_that("regression empirical likelihood is zero at the unconstrained OLS coefficient", {
  observed <- data.frame(
    Y = c(2.1, 1.9, 2.8, 3.0, 2.4, 3.3, 3.1, 4.0),
    R = c(0, 0, 0, 0, 1, 1, 1, 1),
    L = factor(c("a", "b", "a", "b", "a", "b", "a", "b"))
  )

  design <- TruncComp2:::el_regression_design(Y ~ R + L, observed)
  expect_true(design$success)

  profile <- TruncComp2:::el_regression_profile_fit(design, design$estimate)
  expect_true(profile$success)
  expect_equal(profile$statistic, 0, tolerance = 1e-12)
  expect_equal(profile$nuisance,
               as.numeric(design$coefficients[seq_along(design$coefficients) != design$term_index]),
               tolerance = 1e-12)
  expect_true(all(profile$dual$weights > 0))
  expect_equal(sum(profile$dual$weights), 1, tolerance = 1e-10)
  expect_true(max(abs(profile$dual$weighted_moments)) < 1e-10)
})

test_that("regression empirical likelihood is finite and non-negative on regular cases", {
  observed <- data.frame(
    Y = c(2.1, 1.9, 2.8, 3.0, 2.4, 3.3, 3.1, 4.0),
    R = c(0, 0, 0, 0, 1, 1, 1, 1),
    L = factor(c("a", "b", "a", "b", "a", "b", "a", "b"))
  )

  design <- TruncComp2:::el_regression_design(Y ~ R + L, observed)
  statistic <- TruncComp2:::el_regression_statistic(design, design$estimate + 0.25)

  expect_true(is.finite(statistic))
  expect_gte(statistic, 0)
})

test_that("regression empirical likelihood reduces to the two-sample mean-difference EL without covariates", {
  x <- c(2.2, 2.7, 2.9, 3.1, 2.5)
  y <- c(1.8, 2.0, 2.4, 2.6, 2.9)
  observed <- data.frame(Y = c(y, x), R = c(rep(0, length(y)), rep(1, length(x))))

  mean_fit <- TruncComp2:::el_mean_diff_fit(x, y)
  regression_fit <- TruncComp2:::el_regression_fit(observed, Y ~ R, mu = 0)

  expect_true(regression_fit$success)
  expect_equal(regression_fit$estimate, as.numeric(mean_fit$estimate), tolerance = 1e-10)
  expect_equal(regression_fit$statistic, as.numeric(mean_fit$statistic), tolerance = 1e-8)
  expect_equal(as.numeric(regression_fit$conf.int), as.numeric(mean_fit$conf.int), tolerance = 1e-7)
})
