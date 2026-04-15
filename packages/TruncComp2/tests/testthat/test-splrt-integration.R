test_that("SPLRT reproduces the TruncComp2 example outputs without loading EL", {
  expect_false("EL" %in% loadedNamespaces())

  example_data <- loadTruncComp2Example()
  model <- truncComp(Y ~ R, atom = 0, data = example_data, method = "SPLRT")

  expect_false("EL" %in% loadedNamespaces())
  expect_true(model$success)
  expect_equal(model$muDelta, 1.85643, tolerance = 1e-4)
  expect_equal(model$muDeltaCI, c(1.163886, 2.480132), tolerance = 1e-4)
  expect_equal(as.numeric(model$alphaDelta), 0.5238095, tolerance = 1e-7)
  expect_equal(model$alphaDeltaCI, c(0.1660407, 1.59682), tolerance = 1e-6)
  expect_equal(model$W, 31.09545, tolerance = 1e-4)
  expect_lt(abs(model$p - 1.768924e-07), 1e-10)
})

test_that("marginal and simultaneous confidence helpers still work", {
  example_data <- loadTruncComp2Example()
  model <- truncComp(Y ~ R, atom = 0, data = example_data, method = "SPLRT")

  capture.output(ci <- confint(model, type = "marginal"))
  expect_equal(rownames(ci),
               c("Difference in means among the observed:",
                 "Odds ratio of being observed:",
                 "log Odds ratio of being observed:"))
  expect_equal(unname(ci["Difference in means among the observed:", ]), model$muDeltaCI)
  expect_equal(unname(ci["Odds ratio of being observed:", ]), model$alphaDeltaCI)
  expect_equal(unname(ci["log Odds ratio of being observed:", ]), log(model$alphaDeltaCI))

  expect_error(confint(model, type = "marginal", conf.level = 0.9),
               "stored at the fitted confidence level")

  joint <- jointContrastCI(model, plot = FALSE)
  expect_equal(dim(joint$surface), c(35, 35))
  expect_true(all(is.finite(joint$surface)))
})

test_that("joint contrast surface matches the null test and vanishes at the fitted point", {
  example_data <- loadTruncComp2Example()
  model <- truncComp(Y ~ R, atom = 0, data = example_data, method = "SPLRT")

  expect_equal(TruncComp2:::jointContrastLRT(model$data, 0, 0), model$W, tolerance = 1e-6)
  expect_equal(TruncComp2:::jointContrastLRT(model$data, model$muDelta, log(model$alphaDelta)),
               0, tolerance = 1e-6)
})

test_that("existing TruncComp2 data validation still stops before EL helper use", {
  bad_data <- data.frame(
    Y = c(0, 0, 1, 0, 2, 3),
    R = c(0, 0, 0, 1, 1, 1)
  )

  expect_warning(
    fit <- truncComp(Y ~ R, atom = 0, data = bad_data, method = "SPLRT"),
    "data error"
  )
  expect_s3_class(fit, "TruncComp2")
  expect_false(fit$success)
  expect_equal(fit$method, "Semi-empirical Likelihood Ratio Test")

  expect_match(paste(capture.output(summary(fit)), collapse = "\n"),
               "The estimation procedure failed")
  expect_match(paste(capture.output(print(fit)), collapse = "\n"),
               "The estimation procedure failed")
  expect_error(confint(fit), "Estimation failed")
})

test_that("jointContrastCI rejects unsupported fits", {
  example_data <- loadTruncComp2Example()

  splrt_model <- truncComp(Y ~ R, atom = 0, data = example_data, method = "SPLRT")
  expect_silent(jointContrastCI(splrt_model, plot = FALSE, offset = 1, resolution = 5))

  param_model <- splrt_model
  param_model$method <- "Parametric Likelihood Ratio Test"
  expect_error(jointContrastCI(param_model, plot = FALSE, offset = 1, resolution = 5),
               "semi-parametric")

  expect_warning(
    failed_model <- truncComp(Y ~ R, atom = 0,
                              data = data.frame(Y = c(0, 0, 1, 0, 2, 3),
                                                R = c(0, 0, 0, 1, 1, 1)),
                              method = "SPLRT"),
    "data error"
  )
  expect_error(jointContrastCI(failed_model, plot = FALSE, offset = 1, resolution = 5),
               "Estimation failed")
})

test_that("jointContrastCI falls back to finite default grids", {
  example_data <- loadTruncComp2Example()
  model <- truncComp(Y ~ R, atom = 0, data = example_data, method = "SPLRT")

  alpha_fallback <- model
  alpha_fallback$alphaDeltaCI <- c(0, Inf)
  joint_alpha <- jointContrastCI(alpha_fallback, plot = FALSE, offset = 2, resolution = 5)
  expect_equal(joint_alpha$logORdelta,
               seq(log(as.numeric(model$alphaDelta)) - 2,
                   log(as.numeric(model$alphaDelta)) + 2,
                   length.out = 5))

  alpha_symmetric <- model
  alpha_symmetric$alphaDeltaCI <- c(0, Inf)
  alpha_symmetric$alphaDelta <- Inf
  joint_alpha_sym <- jointContrastCI(alpha_symmetric, plot = FALSE, offset = 1.5, resolution = 5)
  expect_equal(joint_alpha_sym$logORdelta, seq(-1.5, 1.5, length.out = 5))

  mu_fallback <- model
  mu_fallback$muDeltaCI <- c(NA_real_, Inf)
  joint_mu <- jointContrastCI(mu_fallback, plot = FALSE, offset = 1.25, resolution = 5)
  expect_equal(joint_mu$muDelta,
               seq(model$muDelta - 1.25, model$muDelta + 1.25, length.out = 5))
})

test_that("simulateTruncatedData validates inputs and returns the canonical shape", {
  f0 <- function(n) rep(1, n)
  f1 <- function(n) rep(2, n)

  simulated <- simulateTruncatedData(5, f0 = f0, f1 = f1, pi0 = 0.4, pi1 = 0.7)
  expect_equal(names(simulated), c("R", "A", "Y"))
  expect_equal(nrow(simulated), 10)
  expect_true(all(simulated$A %in% c(0, 1)))
  expect_true(all(simulated$Y[simulated$A == 0] == 0))

  scalar_simulated <- simulateTruncatedData(4,
                                            f0 = function(n) 1,
                                            f1 = function(n) 2,
                                            pi0 = 1,
                                            pi1 = 1)
  expect_equal(scalar_simulated$Y, c(rep(1, 4), rep(2, 4)))

  expect_error(simulateTruncatedData(0, f0 = f0, f1 = f1, pi0 = 0.4, pi1 = 0.7),
               "positive integer")
  expect_error(simulateTruncatedData(5, f0 = 1, f1 = f1, pi0 = 0.4, pi1 = 0.7),
               "must be functions")
  expect_error(simulateTruncatedData(5, f0 = f0, f1 = f1, pi0 = 1.4, pi1 = 0.7),
               "between 0 and 1")
  expect_error(simulateTruncatedData(5,
                                     f0 = function(n) c(1, 2),
                                     f1 = f1,
                                     pi0 = 0.4,
                                     pi1 = 0.7),
               "length n")
})

test_that("cached logistic profile matches the uncached helper", {
  example_data <- loadTruncComp2Example()
  model <- truncComp(Y ~ R, atom = 0, data = example_data, method = "SPLRT")
  logit_reference <- TruncComp2:::logit.prepare(model$data)

  expect_equal(TruncComp2:::logit.LRT.prepared(logit_reference, 0),
               TruncComp2:::logit.LRT(model$data, 0),
               tolerance = 1e-10)
  expect_equal(TruncComp2:::logit.LRT.prepared(logit_reference, log(as.numeric(model$alphaDelta))),
               TruncComp2:::logit.LRT(model$data, log(as.numeric(model$alphaDelta))),
               tolerance = 1e-10)
})

test_that("cached logistic profile adapts to extreme observation rates", {
  n <- 400
  a <- integer(2 * n)
  a[c(1, 2, n + 1, n + 2)] <- 1
  y <- numeric(2 * n)
  y[c(1, 2, n + 1, n + 2)] <- c(1.0, 1.2, 1.0, 1.2)
  r <- c(rep(0, n), rep(1, n))

  fit <- truncComp.default(y, a, r, method = "SPLRT")
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
  expect_equal(TruncComp2:::jointContrastLRT(fit$data, fit$muDelta, log(as.numeric(fit$alphaDelta))),
               0,
               tolerance = 1e-10)
})

test_that("loadTruncComp2Example returns the packaged example data", {
  example_data <- loadTruncComp2Example()

  expect_s3_class(example_data, "data.frame")
  expect_equal(names(example_data), c("R", "Y"))
  expect_equal(nrow(example_data), 50)
})
