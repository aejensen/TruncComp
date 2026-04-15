test_that("SPLRT reproduces the TruncComp example outputs without loading EL", {
  expect_false("EL" %in% loadedNamespaces())

  data("TruncCompExample", package = "TruncComp")
  model <- truncComp(Y ~ R, atom = 0, data = TruncCompExample, method = "SPLRT")

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
  data("TruncCompExample", package = "TruncComp")
  model <- truncComp(Y ~ R, atom = 0, data = TruncCompExample, method = "SPLRT")

  capture.output(ci <- confint(model, type = "marginal"))
  expect_equal(unname(ci["Difference in means among the observed:", ]), model$muDeltaCI)
  expect_equal(unname(ci["Odds ratio of being observed:", ]), model$alphaDeltaCI)

  joint <- jointContrastCI(model, plot = FALSE, resolution = 5, offset = 1)
  expect_equal(dim(joint$surface), c(5, 5))
  expect_true(all(is.finite(joint$surface)))
})

test_that("joint contrast surface matches the null test and vanishes at the fitted point", {
  data("TruncCompExample", package = "TruncComp")
  model <- truncComp(Y ~ R, atom = 0, data = TruncCompExample, method = "SPLRT")

  expect_equal(TruncComp:::jointContrastLRT(model$data, 0, 0), model$W, tolerance = 1e-6)
  expect_equal(TruncComp:::jointContrastLRT(model$data, model$muDelta, log(model$alphaDelta)),
               0, tolerance = 1e-6)
})

test_that("existing TruncComp data validation still stops before EL helper use", {
  bad_data <- data.frame(
    Y = c(0, 0, 1, 0, 2, 3),
    R = c(0, 0, 0, 1, 1, 1)
  )

  expect_warning(
    fit <- truncComp(Y ~ R, atom = 0, data = bad_data, method = "SPLRT"),
    "data error"
  )
  expect_s3_class(fit, "TruncComp")
  expect_false(fit$success)
  expect_equal(fit$method, "Semi-empirical Likelihood Ratio Test")

  expect_match(paste(capture.output(summary(fit)), collapse = "\n"),
               "The estimation procedure failed")
  expect_match(paste(capture.output(print(fit)), collapse = "\n"),
               "The estimation procedure failed")
  expect_error(confint(fit), "Estimation failed")
})

test_that("jointContrastCI rejects unsupported fits", {
  data("TruncCompExample", package = "TruncComp")

  splrt_model <- truncComp(Y ~ R, atom = 0, data = TruncCompExample, method = "SPLRT")
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
  data("TruncCompExample", package = "TruncComp")
  model <- truncComp(Y ~ R, atom = 0, data = TruncCompExample, method = "SPLRT")

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

test_that("cached logistic profile matches the uncached helper", {
  data("TruncCompExample", package = "TruncComp")
  model <- truncComp(Y ~ R, atom = 0, data = TruncCompExample, method = "SPLRT")
  logit_reference <- TruncComp:::logit.prepare(model$data)

  expect_equal(TruncComp:::logit.LRT.prepared(logit_reference, 0),
               TruncComp:::logit.LRT(model$data, 0),
               tolerance = 1e-10)
  expect_equal(TruncComp:::logit.LRT.prepared(logit_reference, log(as.numeric(model$alphaDelta))),
               TruncComp:::logit.LRT(model$data, log(as.numeric(model$alphaDelta))),
               tolerance = 1e-10)
})
