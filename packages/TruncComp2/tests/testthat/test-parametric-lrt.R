interior_lrt_case <- function(seed, n, mu0, mu1, sigma, pi0, pi1) {
  for(offset in 0:1000) {
    set.seed(seed + offset)
    a0 <- stats::rbinom(n, 1, pi0)
    a1 <- stats::rbinom(n, 1, pi1)

    if(sum(a0) < 2 || sum(a1) < 2 || sum(a0) >= n || sum(a1) >= n) {
      next
    }

    y0 <- numeric(n)
    y1 <- numeric(n)
    y0[a0 == 1] <- stats::rnorm(sum(a0), mu0, sigma)
    y1[a1 == 1] <- stats::rnorm(sum(a1), mu1, sigma)

    if(stats::var(y0[a0 == 1]) <= 0 || stats::var(y1[a1 == 1]) <= 0) {
      next
    }

    return(data.frame(Y = c(y0, y1), A = c(a0, a1), R = c(rep(0, n), rep(1, n))))
  }

  stop("Unable to generate an interior parametric LRT case.")
}

parametric_null_nll <- function(par, data) {
  alpha <- par[1]
  mu <- par[2]
  sigma <- exp(par[3])
  pi <- stats::plogis(alpha)

  y_alive <- data$Y[data$A == 1]
  loglik_a <- sum(data$A * log(pi) + (1 - data$A) * log1p(-pi))
  loglik_y <- sum(stats::dnorm(y_alive, mean = mu, sd = sigma, log = TRUE))

  -(loglik_a + loglik_y)
}

parametric_alt_nll <- function(par, data) {
  alpha0 <- par[1]
  alpha1 <- par[2]
  mu0 <- par[3]
  mu1 <- par[4]
  sigma <- exp(par[5])

  pi <- ifelse(data$R == 0, stats::plogis(alpha0), stats::plogis(alpha1))
  mu <- ifelse(data$R == 0, mu0, mu1)

  loglik_a <- sum(data$A * log(pi) + (1 - data$A) * log1p(-pi))
  loglik_y <- sum(stats::dnorm(data$Y[data$A == 1],
                               mean = mu[data$A == 1],
                               sd = sigma,
                               log = TRUE))

  -(loglik_a + loglik_y)
}

optim_parametric_lrt <- function(data) {
  observed <- data[data$A == 1, , drop = FALSE]
  start_null <- c(stats::qlogis(mean(data$A)),
                  mean(observed$Y),
                  log(stats::sd(observed$Y)))
  start_alt <- c(stats::qlogis(mean(data$A[data$R == 0])),
                 stats::qlogis(mean(data$A[data$R == 1])),
                 mean(observed$Y[observed$R == 0]),
                 mean(observed$Y[observed$R == 1]),
                 log(stats::sd(observed$Y)))

  opt_null <- stats::optim(start_null,
                           parametric_null_nll,
                           data = data,
                           method = "BFGS",
                           control = list(reltol = 1e-12, maxit = 1000))
  opt_alt <- stats::optim(start_alt,
                          parametric_alt_nll,
                          data = data,
                          method = "BFGS",
                          control = list(reltol = 1e-12, maxit = 1000))

  2 * (-opt_alt$value + opt_null$value)
}

runtime_model_reference <- function(data, conf.level = 0.95) {
  observed <- data[data$A == 1, c("Y", "R"), drop = FALSE]

  glm_null <- stats::glm(A ~ 1, family = stats::binomial(), data = data)
  glm_alt <- stats::glm(A ~ R, family = stats::binomial(), data = data)
  lm_null <- stats::lm(Y ~ 1, data = observed)
  lm_alt <- stats::lm(Y ~ R, data = observed)

  z <- stats::qnorm((1 + conf.level) / 2)
  log_or <- stats::coef(glm_alt)[["R"]]
  mu_delta <- stats::coef(lm_alt)[["R"]]

  list(
    W_A = 2 * (as.numeric(stats::logLik(glm_alt)) - as.numeric(stats::logLik(glm_null))),
    W_Y = 2 * (as.numeric(stats::logLik(lm_alt, REML = FALSE)) -
      as.numeric(stats::logLik(lm_null, REML = FALSE))),
    alphaDelta = exp(as.numeric(log_or)),
    muDelta = as.numeric(mu_delta),
    alphaDeltaCI = exp(as.numeric(log_or) + c(-1, 1) * z * sqrt(stats::vcov(glm_alt)["R", "R"])),
    muDeltaCI = as.numeric(mu_delta) + c(-1, 1) * z * sqrt(stats::vcov(lm_alt)["R", "R"])
  )
}

test_that("model-backed parametric engine matches direct glm/lm references on regular data", {
  cases <- list(
    interior_lrt_case(20260411, 12, 2.5, 2.5, 1.0, 0.55, 0.55),
    interior_lrt_case(20260412, 10, 2.2, 3.1, 0.8, 0.45, 0.70),
    interior_lrt_case(20260413, 14, 1.8, 2.8, 1.2, 0.35, 0.65)
  )

  for(case in cases) {
    fit <- TruncComp2:::parametric_lrt_fit(case)
    ref <- runtime_model_reference(case)

    expect_equal(fit$W_A, ref$W_A, tolerance = 1e-10)
    expect_equal(fit$W_Y, ref$W_Y, tolerance = 1e-10)
    expect_equal(fit$W, ref$W_A + ref$W_Y, tolerance = 1e-10)
    expect_equal(fit$alphaDelta, ref$alphaDelta, tolerance = 1e-10)
    expect_equal(fit$muDelta, ref$muDelta, tolerance = 1e-10)
    expect_equal(fit$alphaDeltaCI, ref$alphaDeltaCI, tolerance = 1e-10)
    expect_equal(fit$muDeltaCI, ref$muDeltaCI, tolerance = 1e-10)
    expect_equal(fit$p, stats::pchisq(fit$W, df = 2, lower.tail = FALSE), tolerance = 1e-12)
    expect_gte(fit$W, 0)
  }
})

test_that("model-backed parametric engine matches direct numerical optimization", {
  cases <- list(
    interior_lrt_case(20260421, 12, 2.5, 2.5, 1.0, 0.55, 0.55),
    interior_lrt_case(20260422, 10, 2.2, 3.1, 0.8, 0.45, 0.70),
    interior_lrt_case(20260423, 14, 1.8, 2.8, 1.2, 0.35, 0.65)
  )

  for(case in cases) {
    fit <- TruncComp2:::parametric_lrt_fit(case)
    numeric_lrt <- optim_parametric_lrt(case)

    expect_equal(fit$W, numeric_lrt, tolerance = 1e-6)
    expect_equal(fit$W, fit$W_A + fit$W_Y, tolerance = 1e-12)
    expect_equal(fit$p, stats::pchisq(fit$W, df = 2, lower.tail = FALSE), tolerance = 1e-12)
  }
})

test_that("public LRT path returns the expected object and keeps init compatible", {
  case <- interior_lrt_case(20260424, 12, 2.4, 3.0, 0.9, 0.50, 0.70)
  reference <- runtime_model_reference(case)
  custom_init <- list(alpha = 0, alphaDelta = 0, mu = 0, muDelta = 0, sigma = 1)

  fit <- truncComp(Y ~ R, atom = 0, data = case[, c("Y", "R")], method = "LRT", init = custom_init)

  expect_s3_class(fit, "TruncComp2")
  expect_true(fit$success)
  expect_equal(fit$method, "Parametric Likelihood Ratio Test")
  expect_equal(fit$DeltaCI, c(NA_real_, NA_real_))
  expect_equal(fit$init, custom_init)
  expect_equal(fit$muDelta, reference$muDelta, tolerance = 1e-10)
  expect_equal(fit$alphaDelta, reference$alphaDelta, tolerance = 1e-10)
  expect_equal(fit$W, reference$W_A + reference$W_Y, tolerance = 1e-10)
  expect_equal(fit$muDeltaCI, reference$muDeltaCI, tolerance = 1e-10)
  expect_equal(fit$alphaDeltaCI, reference$alphaDeltaCI, tolerance = 1e-10)

  expect_match(paste(capture.output(summary(fit)), collapse = "\n"), "Joint test statistic")
  expect_match(paste(capture.output(print(fit)), collapse = "\n"), "Parametric Likelihood Ratio Test")

  capture.output(ci <- confint(fit, type = "marginal"))
  expect_equal(unname(ci["Difference in means among the observed:", ]), fit$muDeltaCI)
  expect_equal(unname(ci["Odds ratio of being observed:", ]), fit$alphaDeltaCI)
})

test_that("LRT boundary cases keep the statistic even when Wald intervals are undefined", {
  same_constants <- data.frame(
    Y = c(1, 1, 0, 1, 1, 0),
    A = c(1, 1, 0, 1, 1, 0),
    R = c(0, 0, 0, 1, 1, 1)
  )
  same_fit <- truncComp.default(same_constants$Y, same_constants$A, same_constants$R, method = "LRT")
  expect_true(same_fit$success)
  expect_equal(same_fit$muDelta, 0)
  expect_equal(same_fit$alphaDelta, 1)
  expect_equal(same_fit$W, 0)
  expect_true(all(is.na(same_fit$muDeltaCI)))

  separated_constants <- data.frame(
    Y = c(1, 1, 0, 2, 2, 0),
    A = c(1, 1, 0, 1, 1, 0),
    R = c(0, 0, 0, 1, 1, 1)
  )
  separated_fit <- truncComp.default(separated_constants$Y,
                                     separated_constants$A,
                                     separated_constants$R,
                                     method = "LRT")
  expect_true(separated_fit$success)
  expect_equal(separated_fit$muDelta, 1)
  expect_true(is.infinite(separated_fit$W))
  expect_equal(separated_fit$p, 0)
  expect_true(all(is.na(separated_fit$muDeltaCI)))

  all_observed_arm <- data.frame(
    Y = c(1.0, 1.2, 1.4, 1.6, 0, 2.0, 2.3, 2.5),
    A = c(1, 1, 1, 1, 0, 1, 1, 1),
    R = c(0, 0, 0, 0, 1, 1, 1, 1)
  )
  boundary_fit <- truncComp.default(all_observed_arm$Y,
                                    all_observed_arm$A,
                                    all_observed_arm$R,
                                    method = "LRT")
  expect_true(boundary_fit$success)
  expect_equal(boundary_fit$alphaDelta, 0)
  expect_true(all(is.na(boundary_fit$alphaDeltaCI)))
  expect_true(is.finite(boundary_fit$W) || is.infinite(boundary_fit$W))
})

test_that("LRT still returns a failed TruncComp2 object on invalid data", {
  expect_warning(
    fit <- truncComp.default(c(0, 1, 0, 2),
                             c(0, 1, 0, 1),
                             c(0, 0, 1, 1),
                             method = "LRT"),
    "data error"
  )

  expect_s3_class(fit, "TruncComp2")
  expect_false(fit$success)
  expect_error(confint(fit), "Estimation failed")
})
