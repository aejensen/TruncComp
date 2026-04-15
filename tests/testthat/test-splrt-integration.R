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
  expect_false(fit$success)
})
