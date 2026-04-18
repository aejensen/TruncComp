test_that("SPLRT reproduces the TruncComp2 example outputs without loading EL", {
  expect_false("EL" %in% loadedNamespaces())

  example_data <- trunc_comp_example
  model <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = "splrt")

  expect_false("EL" %in% loadedNamespaces())
  expect_true(model$success)
  expect_equal(model$mu_delta, 1.85643, tolerance = 1e-4)
  expect_equal(model$mu_delta_ci, c(1.163886, 2.480132), tolerance = 1e-4)
  expect_equal(as.numeric(model$alpha_delta), 0.5238095, tolerance = 1e-7)
  expect_equal(model$alpha_delta_ci, c(0.1660407, 1.59682), tolerance = 1e-6)
  expect_equal(model$atom, 0)
  expect_false(any(c("DeltaCI", "DeltaMarginalCI", "DeltaProjectedCI", "DeltaProfileCI") %in% names(model)))
  expect_equal(model$statistic, 31.09545, tolerance = 1e-4)
  expect_lt(abs(model$p.value - 1.768924e-07), 1e-10)
})

test_that("marginal and simultaneous confidence helpers still work", {
  example_data <- trunc_comp_example
  model <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = "splrt")

  capture.output(ci <- confint(model))
  expect_equal(rownames(ci),
               c("mu_delta",
                 "alpha_delta"))
  expect_equal(unname(ci["mu_delta", ]), model$mu_delta_ci)
  expect_equal(unname(ci["alpha_delta", ]), model$alpha_delta_ci)

  capture.output(delta_welch <- confint(model, parameter = "delta", method = "welch"))
  expect_equal(rownames(delta_welch), "delta (welch)")
  expect_lte(unname(delta_welch[1, 1]), model$delta)
  expect_gte(unname(delta_welch[1, 2]), model$delta)

  capture.output(delta_projected <- confint(model, parameter = "delta", method = "projected", plot = FALSE))
  expect_equal(rownames(delta_projected), "delta (projected)")
  expect_lte(unname(delta_projected[1, 1]), model$delta)
  expect_gte(unname(delta_projected[1, 2]), model$delta)

  capture.output(delta_profile <- confint(model, parameter = "delta", method = "profile", plot = FALSE))
  expect_equal(rownames(delta_profile), "delta (profile)")
  expect_lte(unname(delta_profile[1, 1]), model$delta)
  expect_gte(unname(delta_profile[1, 2]), model$delta)

  expect_error(confint(model, conf.level = 0.9),
               "stored at the fitted confidence level")

  joint <- joint_contrast_surface(model, plot = FALSE)
  expect_equal(dim(joint$surface), c(35, 35))
  expect_true(all(is.finite(joint$surface)))
})

test_that("joint contrast surface matches the null test and vanishes at the fitted point", {
  example_data <- trunc_comp_example
  model <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = "splrt")

  expect_equal(TruncComp2:::jointContrastLRT(model$data, 0, 0), model$statistic, tolerance = 1e-6)
  expect_equal(TruncComp2:::jointContrastLRT(model$data, model$mu_delta, log(model$alpha_delta)),
               0, tolerance = 1e-6)
})

test_that("joint_contrast_surface returns invisibly when plotting", {
  example_data <- trunc_comp_example
  model <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = "splrt")

  plotted <- withVisible(joint_contrast_surface(model, plot = TRUE, resolution = 5))
  unplotted <- withVisible(joint_contrast_surface(model, plot = FALSE, resolution = 5))

  expect_false(plotted$visible)
  expect_true(unplotted$visible)
  expect_equal(dim(plotted$value$surface), c(5, 5))
  expect_equal(dim(unplotted$value$surface), c(5, 5))
})

test_that("existing TruncComp2 data validation still stops before EL helper use", {
  bad_data <- data.frame(
    Y = c(0, 0, 1, 0, 2, 3),
    R = c(0, 0, 0, 1, 1, 1)
  )

  expect_warning(
    fit <- trunc_comp(Y ~ R, atom = 0, data = bad_data, method = "splrt"),
    "data error"
  )
  expect_s3_class(fit, "trunc_comp_fit")
  expect_false(fit$success)
  expect_equal(fit$method, "splrt")

  expect_match(paste(capture.output(summary(fit)), collapse = "\n"),
               "The estimation procedure failed")
  expect_match(paste(capture.output(print(fit)), collapse = "\n"),
               "Error: Estimation failed due to data error")
  expect_error(confint(fit), "Estimation failed")
})

test_that("joint_contrast_surface rejects unsupported fits", {
  example_data <- trunc_comp_example

  splrt_model <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = "splrt")
  expect_silent(joint_contrast_surface(splrt_model, plot = FALSE, offset = 1, resolution = 5))

  expect_warning(
    failed_model <- trunc_comp(Y ~ R, atom = 0,
                               data = data.frame(Y = c(0, 0, 1, 0, 2, 3),
                                                 R = c(0, 0, 0, 1, 1, 1)),
                               method = "splrt"),
    "data error"
  )
  expect_error(joint_contrast_surface(failed_model, plot = FALSE, offset = 1, resolution = 5),
               "Estimation failed")
  expect_error(confint(failed_model, parameter = "delta", method = "projected"), "Estimation failed")
  expect_error(confint(failed_model, parameter = "delta", method = "profile"), "Estimation failed")
})

test_that("joint_contrast_surface falls back to finite default grids", {
  example_data <- trunc_comp_example
  model <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = "splrt")

  alpha_fallback <- model
  alpha_fallback$alpha_delta_ci <- c(0, Inf)
  joint_alpha <- joint_contrast_surface(alpha_fallback, plot = FALSE, offset = 2, resolution = 5)
  expect_equal(joint_alpha$log_or_delta,
               seq(log(as.numeric(model$alpha_delta)) - 2,
                   log(as.numeric(model$alpha_delta)) + 2,
                   length.out = 5))

  alpha_symmetric <- model
  alpha_symmetric$alpha_delta_ci <- c(0, Inf)
  alpha_symmetric$alpha_delta <- Inf
  joint_alpha_sym <- joint_contrast_surface(alpha_symmetric, plot = FALSE, offset = 1.5, resolution = 5)
  expect_equal(joint_alpha_sym$log_or_delta, seq(-1.5, 1.5, length.out = 5))

  mu_fallback <- model
  mu_fallback$mu_delta_ci <- c(NA_real_, Inf)
  joint_mu <- joint_contrast_surface(mu_fallback, plot = FALSE, offset = 1.25, resolution = 5)
  expect_equal(joint_mu$mu_delta,
               seq(model$mu_delta - 1.25, model$mu_delta + 1.25, length.out = 5))
})

test_that("joint_contrast_surface derives default offsets from fitted SPLRT data", {
  example_data <- trunc_comp_example
  model <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = "splrt")

  offsets <- TruncComp2:::jointContrastDefaultOffsets(model)
  joint <- joint_contrast_surface(model, plot = FALSE, resolution = 5)

  expect_equal(joint$mu_delta,
               seq(model$mu_delta_ci[1] - offsets[1],
                   model$mu_delta_ci[2] + offsets[1],
                   length.out = 5))
  expect_equal(joint$log_or_delta,
               seq(log(model$alpha_delta_ci[1]) - offsets[2],
                   log(model$alpha_delta_ci[2]) + offsets[2],
                   length.out = 5))

  vector_joint <- joint_contrast_surface(model, plot = FALSE, offset = c(0.5, 0.25), resolution = 5)
  expect_equal(vector_joint$mu_delta,
               seq(model$mu_delta_ci[1] - 0.5,
                   model$mu_delta_ci[2] + 0.5,
                   length.out = 5))
  expect_equal(vector_joint$log_or_delta,
               seq(log(model$alpha_delta_ci[1]) - 0.25,
                   log(model$alpha_delta_ci[2]) + 0.25,
                   length.out = 5))
})

test_that("joint contrast plotting helper returns a ggplot object", {
  example_data <- trunc_comp_example
  model <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = "splrt")
  joint <- TruncComp2:::jointContrastSurfaceData(model, plot = FALSE, resolution = 5)

  plot_obj <- TruncComp2:::jointContrastPlot(
    joint$mu_delta,
    joint$log_or_delta,
    joint$surface,
    model,
    model$conf.level
  )

  expect_s3_class(plot_obj, "ggplot")
})

test_that("simulate_truncated_data validates inputs and returns the canonical shape", {
  f0 <- function(n) rep(1, n)
  f1 <- function(n) rep(2, n)

  simulated <- simulate_truncated_data(5, f0 = f0, f1 = f1, pi0 = 0.4, pi1 = 0.7)
  expect_equal(names(simulated), c("R", "A", "Y"))
  expect_equal(nrow(simulated), 10)
  expect_true(all(simulated$A %in% c(0, 1)))
  expect_true(all(simulated$Y[simulated$A == 0] == 0))
  expect_true(all(simulated$Y[simulated$A == 1] != 0))

  atom_simulated <- simulate_truncated_data(5, f0 = f0, f1 = f1, pi0 = 0.4, pi1 = 0.7, atom = -3)
  expect_true(all(atom_simulated$Y[atom_simulated$A == 0] == -3))
  expect_true(all(atom_simulated$Y[atom_simulated$A == 1] != -3))

  scalar_simulated <- simulate_truncated_data(4,
                                            f0 = function(n) 1,
                                            f1 = function(n) 2,
                                            pi0 = 1,
                                            pi1 = 1)
  expect_equal(scalar_simulated$Y, c(rep(1, 4), rep(2, 4)))

  set.seed(123)
  random_scalar_simulated <- simulate_truncated_data(
    5,
    f0 = function(n) stats::rnorm(1),
    f1 = function(n) stats::rnorm(1),
    pi0 = 1,
    pi1 = 1
  )
  set.seed(123)
  expected_scalar_draws <- c(replicate(5, stats::rnorm(1)),
                             replicate(5, stats::rnorm(1)))
  expect_equal(random_scalar_simulated$Y, expected_scalar_draws)

  expect_error(
    simulate_truncated_data(5,
                          f0 = function(n) rep(0, n),
                          f1 = f1,
                          pi0 = 1,
                          pi1 = 1),
    "must not return the atom value"
  )

  expect_error(simulate_truncated_data(0, f0 = f0, f1 = f1, pi0 = 0.4, pi1 = 0.7),
               "positive integer")
  expect_error(simulate_truncated_data(5, f0 = 1, f1 = f1, pi0 = 0.4, pi1 = 0.7),
               "must be functions")
  expect_error(simulate_truncated_data(5, f0 = f0, f1 = f1, pi0 = 1.4, pi1 = 0.7),
               "between 0 and 1")
  expect_error(simulate_truncated_data(5, f0 = f0, f1 = f1, pi0 = c(0.4, 0.5), pi1 = 0.7),
               "between 0 and 1")
  expect_error(simulate_truncated_data(5,
                                     f0 = function(n) c(1, 2),
                                     f1 = f1,
                                     pi0 = 0.4,
                                     pi1 = 0.7),
               "length n")

  expect_error(TruncComp2:::simTruncData(5, mu0 = c(0, 1), mu1 = 2, pi0 = 0.4, pi1 = 0.7),
               "mu0 must be a single finite numeric value")
  expect_error(TruncComp2:::simTruncData(5, mu0 = 0, mu1 = 2, pi0 = 0.4, pi1 = 0.7, sigma = -1),
               "sigma must be a single positive finite numeric value")
  expect_error(TruncComp2:::simTruncData(5, mu0 = 0, mu1 = 2, pi0 = 0.4, pi1 = 0.7, dist = "t-sq", df = 1),
               "df must be a single finite numeric value greater than 1")
})

test_that("cached logistic profile matches the uncached helper", {
  example_data <- trunc_comp_example
  model <- trunc_comp(Y ~ R, atom = 0, data = example_data, method = "splrt")
  logit_reference <- TruncComp2:::logit.prepare(model$data)

  expect_equal(TruncComp2:::logit.LRT.prepared(logit_reference, 0),
               TruncComp2:::logit.LRT(model$data, 0),
               tolerance = 1e-10)
  expect_equal(TruncComp2:::logit.LRT.prepared(logit_reference, log(as.numeric(model$alpha_delta))),
               TruncComp2:::logit.LRT(model$data, log(as.numeric(model$alpha_delta))),
               tolerance = 1e-10)
})

test_that("cached logistic profile adapts to extreme observation rates", {
  n <- 400
  a <- integer(2 * n)
  a[c(1, 2, n + 1, n + 2)] <- 1
  y <- numeric(2 * n)
  y[c(1, 2, n + 1, n + 2)] <- c(1.0, 1.2, 1.0, 1.2)
  r <- c(rep(0, n), rep(1, n))

  fit <- trunc_comp(y, a, r, method = "splrt")
  expect_true(fit$success)

  logit_reference <- TruncComp2:::logit.prepare(fit$data)
  wide_profile <- stats::optimize(
    function(b0) TruncComp2:::logit.likelihood.prepared(logit_reference, c(b0, 0)),
    interval = c(-20, 20),
    maximum = FALSE
  )

  expect_equal(TruncComp2:::logit.likelihood.profile.prepared(logit_reference, 0),
               wide_profile$objective,
               tolerance = 1e-10)
  expect_equal(TruncComp2:::logit.LRT.prepared(logit_reference, 0), 0, tolerance = 1e-10)
  expect_equal(TruncComp2:::jointContrastLRT(fit$data, fit$mu_delta, log(as.numeric(fit$alpha_delta))),
               0,
               tolerance = 1e-10)
})

test_that("trunc_comp_example exposes the packaged example data", {
  example_data <- trunc_comp_example

  expect_s3_class(example_data, "data.frame")
  expect_equal(names(example_data), c("R", "Y"))
  expect_equal(nrow(example_data), 50)
})
