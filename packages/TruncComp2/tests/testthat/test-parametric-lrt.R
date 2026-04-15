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
  analytic <- TruncComp2:::parametric_lrt_fit(data)
  start_null <- c(stats::qlogis(analytic$bernoulli$pi),
                  analytic$normal$mu,
                  log(sqrt(analytic$normal$sigma2_hat)))
  start_alt <- c(stats::qlogis(analytic$bernoulli$pi0),
                 stats::qlogis(analytic$bernoulli$pi1),
                 analytic$normal$mu0,
                 analytic$normal$mu1,
                 log(sqrt(analytic$normal$sigma2_hat)))

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

test_that("parametric helpers handle safe logs and grouped Bernoulli likelihoods", {
  expect_equal(TruncComp2:::parametric_safe_xlogy(c(0, 2), c(0, 0.5)),
               c(0, 2 * log(0.5)))
  expect_true(is.infinite(TruncComp2:::parametric_safe_xlogy(1, 0)))

  a <- c(rep(1, 3), rep(0, 7), rep(1, 6), rep(0, 6))
  r <- c(rep(0, 10), rep(1, 12))
  fit <- TruncComp2:::parametric_bernoulli_summary(a, r, conf.level = 0.95)

  pi0 <- 3 / 10
  pi1 <- 6 / 12
  pi <- 9 / 22
  expected <- 2 * (
    (3 * log(pi0) + 7 * log(1 - pi0) + 6 * log(pi1) + 6 * log(1 - pi1)) -
      (9 * log(pi) + 13 * log(1 - pi))
  )

  expect_equal(fit$W, expected, tolerance = 1e-12)
  expect_equal(fit$alphaDelta, (pi1 / (1 - pi1)) / (pi0 / (1 - pi0)), tolerance = 1e-12)
})

test_that("parametric normal helper matches the closed-form LR and boundary rules", {
  y0 <- c(1, 2, 4)
  y1 <- c(2, 3, 5)
  fit <- TruncComp2:::parametric_normal_summary(y0, y1, conf.level = 0.95)

  sse1 <- sum((y0 - mean(y0))^2) + sum((y1 - mean(y1))^2)
  sse0 <- sum((c(y0, y1) - mean(c(y0, y1)))^2)
  expected <- length(c(y0, y1)) * log(sse0 / sse1)

  expect_equal(fit$W, expected, tolerance = 1e-12)
  expect_equal(fit$muDelta, mean(y1) - mean(y0), tolerance = 1e-12)

  swapped <- TruncComp2:::parametric_normal_summary(y1, y0, conf.level = 0.95)
  expect_equal(swapped$W, fit$W, tolerance = 1e-12)
  expect_equal(swapped$muDelta, -fit$muDelta, tolerance = 1e-12)

  same_constants <- TruncComp2:::parametric_normal_summary(c(1, 1), c(1, 1), conf.level = 0.95)
  expect_equal(same_constants$W, 0)
  expect_true(all(is.na(same_constants$muDeltaCI)))

  separated_constants <- TruncComp2:::parametric_normal_summary(c(1, 1), c(2, 2), conf.level = 0.95)
  expect_true(is.infinite(separated_constants$W))
  expect_true(all(is.na(separated_constants$muDeltaCI)))
})

test_that("analytic parametric engine matches direct numerical optimization", {
  cases <- list(
    interior_lrt_case(20260411, 12, 2.5, 2.5, 1.0, 0.55, 0.55),
    interior_lrt_case(20260412, 10, 2.2, 3.1, 0.8, 0.45, 0.70),
    interior_lrt_case(20260413, 14, 1.8, 2.8, 1.2, 0.35, 0.65)
  )

  for(case in cases) {
    analytic <- TruncComp2:::parametric_lrt_fit(case)
    numeric_lrt <- optim_parametric_lrt(case)

    expect_equal(analytic$W, numeric_lrt, tolerance = 1e-6)
    expect_equal(analytic$W, analytic$W_A + analytic$W_Y, tolerance = 1e-12)
    expect_equal(analytic$p, stats::pchisq(analytic$W, df = 2, lower.tail = FALSE), tolerance = 1e-12)
    expect_gte(analytic$W, 0)
  }
})

test_that("analytic parametric engine matches frozen optimizer fixtures on regular cases", {
  fixture <- readRDS(fixture_path("lrt_reference.rds"))

  expect_equal(fixture$generated_on, as.Date("2026-04-15"))
  expect_equal(fixture$reference_package, "bbmle")

  for(case in fixture$cases) {
    fit <- truncComp(Y ~ R, atom = 0, data = case$data, method = "LRT")
    ref <- case$reference

    expect_true(fit$success, info = case$name)
    expect_equal(fit$muDelta, ref$muDelta, tolerance = 1e-5, info = case$name)
    expect_equal(as.numeric(fit$alphaDelta), ref$alphaDelta, tolerance = 2e-4, info = case$name)
    expect_equal(fit$W, ref$W, tolerance = 2e-6, info = case$name)
    expect_equal(fit$p, ref$p, tolerance = 1e-6, info = case$name)
    expect_equal(fit$muDeltaCI, ref$muDeltaCI, tolerance = 1e-5, info = case$name)
    expect_equal(fit$alphaDeltaCI, ref$alphaDeltaCI, tolerance = 5e-4, info = case$name)
  }
})

test_that("public LRT path returns the expected object and keeps init compatible", {
  fixture <- readRDS(fixture_path("lrt_reference.rds"))
  case <- fixture$cases[[1]]
  custom_init <- list(alpha = 0, alphaDelta = 0, mu = 0, muDelta = 0, sigma = 1)

  fit <- truncComp(Y ~ R, atom = 0, data = case$data, method = "LRT", init = custom_init)

  expect_s3_class(fit, "TruncComp2")
  expect_true(fit$success)
  expect_equal(fit$method, "Parametric Likelihood Ratio Test")
  expect_equal(fit$DeltaCI, c(NA_real_, NA_real_))
  expect_equal(fit$init, custom_init)

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
  expect_true(is.infinite(boundary_fit$alphaDelta) || boundary_fit$alphaDelta == 0)
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
