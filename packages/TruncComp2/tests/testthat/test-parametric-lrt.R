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

adjusted_lrt_case <- function(seed, n, treatment_effect = 0.7, confounded = FALSE) {
  for(offset in 0:1000) {
    set.seed(seed + offset)
    total_n <- 2 * n
    l1 <- stats::rnorm(total_n)
    l2 <- factor(sample(c("a", "b"), total_n, replace = TRUE))

    if(confounded) {
      r <- stats::rbinom(total_n, 1, stats::plogis(0.9 * l1 + 0.8 * (l2 == "b")))
      if(sum(r == 0) < n / 2 || sum(r == 1) < n / 2) {
        next
      }
    } else {
      r <- c(rep(0, n), rep(1, n))
    }

    lp_a <- -0.4 + 0.5 * r + 0.7 * l1 - 0.5 * (l2 == "b")
    a <- stats::rbinom(total_n, 1, stats::plogis(lp_a))

    if(sum(a[r == 0]) < 3 || sum(a[r == 1]) < 3 || all(a == 1)) {
      next
    }

    mu <- 1.6 + treatment_effect * r + 0.6 * l1 + 0.5 * (l2 == "b")
    y <- numeric(total_n)
    y[a == 1] <- stats::rnorm(sum(a), mu[a == 1], 0.7)

    observed <- data.frame(Y = y[a == 1], R = r[a == 1], L1 = l1[a == 1], L2 = l2[a == 1])
    gaussian_fit <- stats::lm(Y ~ R + L1 + L2, data = observed)
    bernoulli_fit <- suppressWarnings(stats::glm(A ~ R + L1 + L2,
                                                 family = stats::binomial(),
                                                 data = data.frame(A = a, R = r, L1 = l1, L2 = l2)))

    if(anyNA(stats::coef(gaussian_fit)) ||
       anyNA(stats::coef(bernoulli_fit)) ||
       summary(gaussian_fit)$sigma <= 0) {
      next
    }

    fitted_probs <- stats::fitted(bernoulli_fit)
    if(any(!is.finite(fitted_probs)) || any(fitted_probs <= 1e-6 | fitted_probs >= 1 - 1e-6)) {
      next
    }

    return(data.frame(Y = y, A = a, R = r, L1 = l1, L2 = l2))
  }

  stop("Unable to generate an interior adjusted parametric LRT case.")
}

model_data_for_formula <- function(data, adjust = NULL) {
  variables <- c("Y", "R")
  if(!is.null(adjust)) {
    variables <- unique(c(variables, all.vars(adjust)))
  }

  data[, variables, drop = FALSE]
}

runtime_model_reference <- function(data, conf.level = 0.95, adjust = NULL) {
  formulas <- TruncComp2:::parametric_model_formulas(adjust)
  observed_variables <- c("Y", "R")
  if(!is.null(adjust)) {
    observed_variables <- unique(c(observed_variables, all.vars(adjust)))
  }
  observed <- data[data$A == 1, observed_variables, drop = FALSE]

  glm_null <- suppressWarnings(stats::glm(formulas$bernoulli_null,
                                          family = stats::binomial(),
                                          data = data))
  glm_alt <- suppressWarnings(stats::glm(formulas$bernoulli_alt,
                                         family = stats::binomial(),
                                         data = data))
  lm_null <- stats::lm(formulas$normal_null, data = observed)
  lm_alt <- stats::lm(formulas$normal_alt, data = observed)

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

runtime_parametric_joint_reference <- function(data, muDelta, logORdelta) {
  observed <- droplevels(data[data$A == 1, c("Y", "R"), drop = FALSE])

  glm_alt <- suppressWarnings(stats::glm(A ~ R, family = stats::binomial(), data = data))
  glm_constrained <- suppressWarnings(stats::glm(
    A ~ 1,
    family = stats::binomial(),
    data = data,
    offset = logORdelta * data$R
  ))
  lm_alt <- suppressWarnings(stats::lm(Y ~ R, data = observed))
  lm_constrained <- suppressWarnings(stats::lm(
    Y ~ 1,
    data = observed,
    offset = muDelta * observed$R
  ))

  W_A <- 2 * (as.numeric(stats::logLik(glm_alt)) - as.numeric(stats::logLik(glm_constrained)))
  W_Y <- 2 * (as.numeric(stats::logLik(lm_alt, REML = FALSE)) -
    as.numeric(stats::logLik(lm_constrained, REML = FALSE)))

  list(
    W_A = W_A,
    W_Y = W_Y,
    W = W_A + W_Y
  )
}

parametric_null_nll <- function(par, data, adjust = NULL) {
  formulas <- TruncComp2:::parametric_model_formulas(adjust)
  adjust_vars <- if(is.null(adjust)) character(0) else all.vars(adjust)
  bernoulli_x <- stats::model.matrix(formulas$bernoulli_null, data = data)
  observed <- data[data$A == 1, unique(c("Y", "R", adjust_vars)), drop = FALSE]
  gaussian_x <- stats::model.matrix(formulas$normal_null, data = observed)

  n_bernoulli <- ncol(bernoulli_x)
  n_gaussian <- ncol(gaussian_x)

  beta_bernoulli <- par[seq_len(n_bernoulli)]
  beta_gaussian <- par[n_bernoulli + seq_len(n_gaussian)]
  sigma <- exp(par[n_bernoulli + n_gaussian + 1])

  pi <- stats::plogis(drop(bernoulli_x %*% beta_bernoulli))
  mu <- drop(gaussian_x %*% beta_gaussian)

  loglik_a <- sum(data$A * log(pi) + (1 - data$A) * log1p(-pi))
  loglik_y <- sum(stats::dnorm(observed$Y, mean = mu, sd = sigma, log = TRUE))

  -(loglik_a + loglik_y)
}

parametric_alt_nll <- function(par, data, adjust = NULL) {
  formulas <- TruncComp2:::parametric_model_formulas(adjust)
  adjust_vars <- if(is.null(adjust)) character(0) else all.vars(adjust)
  bernoulli_x <- stats::model.matrix(formulas$bernoulli_alt, data = data)
  observed <- data[data$A == 1, unique(c("Y", "R", adjust_vars)), drop = FALSE]
  gaussian_x <- stats::model.matrix(formulas$normal_alt, data = observed)

  n_bernoulli <- ncol(bernoulli_x)
  n_gaussian <- ncol(gaussian_x)

  beta_bernoulli <- par[seq_len(n_bernoulli)]
  beta_gaussian <- par[n_bernoulli + seq_len(n_gaussian)]
  sigma <- exp(par[n_bernoulli + n_gaussian + 1])

  pi <- stats::plogis(drop(bernoulli_x %*% beta_bernoulli))
  mu <- drop(gaussian_x %*% beta_gaussian)

  loglik_a <- sum(data$A * log(pi) + (1 - data$A) * log1p(-pi))
  loglik_y <- sum(stats::dnorm(observed$Y, mean = mu, sd = sigma, log = TRUE))

  -(loglik_a + loglik_y)
}

optim_parametric_lrt <- function(data, adjust = NULL) {
  fits <- TruncComp2:::parametric_fit_models(data, adjust = adjust)
  sigma_start <- log(summary(fits$normal_alt)$sigma)

  start_null <- c(unname(stats::coef(fits$bernoulli_null)),
                  unname(stats::coef(fits$normal_null)),
                  sigma_start)
  start_alt <- c(unname(stats::coef(fits$bernoulli_alt)),
                 unname(stats::coef(fits$normal_alt)),
                 sigma_start)

  opt_null <- stats::optim(start_null,
                           parametric_null_nll,
                           data = data,
                           adjust = adjust,
                           method = "BFGS",
                           control = list(reltol = 1e-12, maxit = 5000))
  opt_alt <- stats::optim(start_alt,
                          parametric_alt_nll,
                          data = data,
                          adjust = adjust,
                          method = "BFGS",
                          control = list(reltol = 1e-12, maxit = 5000))

  2 * (-opt_alt$value + opt_null$value)
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

    expect_true(fit$success)
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

    expect_true(fit$success)
    expect_equal(fit$W, numeric_lrt, tolerance = 1e-6)
    expect_equal(fit$W, fit$W_A + fit$W_Y, tolerance = 1e-12)
    expect_equal(fit$p, stats::pchisq(fit$W, df = 2, lower.tail = FALSE), tolerance = 1e-12)
  }
})

test_that("adjusted parametric engine matches direct glm/lm references", {
  cases <- list(
    list(data = adjusted_lrt_case(20260501, 16), adjust = ~ L1),
    list(data = adjusted_lrt_case(20260502, 16), adjust = ~ L2),
    list(data = adjusted_lrt_case(20260503, 18), adjust = ~ L1 + L2)
  )

  for(case in cases) {
    fit <- TruncComp2:::parametric_lrt_fit(case$data, adjust = case$adjust)
    ref <- runtime_model_reference(case$data, adjust = case$adjust)

    expect_true(fit$success)
    expect_equal(fit$W_A, ref$W_A, tolerance = 1e-10)
    expect_equal(fit$W_Y, ref$W_Y, tolerance = 1e-10)
    expect_equal(fit$W, ref$W_A + ref$W_Y, tolerance = 1e-10)
    expect_equal(fit$alphaDelta, ref$alphaDelta, tolerance = 1e-10)
    expect_equal(fit$muDelta, ref$muDelta, tolerance = 1e-10)
    expect_equal(fit$alphaDeltaCI, ref$alphaDeltaCI, tolerance = 1e-10)
    expect_equal(fit$muDeltaCI, ref$muDeltaCI, tolerance = 1e-10)
    expect_true(is.na(fit$Delta))
    expect_equal(fit$DeltaCI, c(NA_real_, NA_real_))
  }
})

test_that("adjusted parametric engine matches direct numerical optimization", {
  cases <- list(
    list(data = adjusted_lrt_case(20260511, 16), adjust = ~ L1),
    list(data = adjusted_lrt_case(20260512, 18), adjust = ~ L1 + L2)
  )

  for(case in cases) {
    fit <- TruncComp2:::parametric_lrt_fit(case$data, adjust = case$adjust)
    numeric_lrt <- optim_parametric_lrt(case$data, adjust = case$adjust)

    expect_true(fit$success)
    expect_equal(fit$W, numeric_lrt, tolerance = 1e-6)
    expect_equal(fit$W, fit$W_A + fit$W_Y, tolerance = 1e-12)
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

test_that("parametric simultaneous confidence helpers work for unadjusted fits", {
  case <- interior_lrt_case(20260601, 16, 2.1, 2.8, 0.75, 0.45, 0.70)
  fit <- truncComp(Y ~ R, atom = 0, data = case[, c("Y", "R")], method = "LRT")

  expect_true(fit$success)

  capture.output(joint <- confint(fit, type = "simultaneous", plot = FALSE, offset = 1, resolution = 5))
  expect_equal(dim(joint$surface), c(5, 5))
  expect_equal(length(joint$muDelta), 5)
  expect_equal(length(joint$logORdelta), 5)

  direct_joint <- jointContrastCI(fit, plot = FALSE, offset = 1, resolution = 5)
  expect_equal(direct_joint$surface, joint$surface, tolerance = 1e-10)

  parametric_reference <- TruncComp2:::parametricJointReference(fit$data)
  fitted_log_or <- log(as.numeric(fit$alphaDelta))

  expect_equal(
    TruncComp2:::jointContrastLRT.parametric.cached(parametric_reference, 0, 0),
    fit$W,
    tolerance = 1e-10
  )
  expect_equal(
    TruncComp2:::jointContrastLRT.parametric.cached(parametric_reference, fit$muDelta, fitted_log_or),
    0,
    tolerance = 1e-10
  )
})

test_that("parametric simultaneous surface matches direct constrained model fits", {
  cases <- list(
    list(data = interior_lrt_case(20260611, 18, 2.3, 2.9, 0.8, 0.35, 0.65),
         mu = 0.2, log_or = -0.15),
    list(data = interior_lrt_case(20260612, 45, 2.0, 2.2, 0.7, 0.06, 0.12),
         mu = 0.1, log_or = 0.05),
    list(data = interior_lrt_case(20260613, 20, 2.5, 2.7, 0.08, 0.50, 0.65),
         mu = -0.03, log_or = 0.12)
  )

  for(case in cases) {
    fit <- truncComp.default(case$data$Y, case$data$A, case$data$R, method = "LRT")
    expect_true(fit$success)

    parametric_reference <- TruncComp2:::parametricJointReference(case$data)
    runtime_reference <- runtime_parametric_joint_reference(case$data, case$mu, case$log_or)

    expect_equal(
      TruncComp2:::jointContrastLRT.parametric.cached(parametric_reference, case$mu, case$log_or),
      runtime_reference$W,
      tolerance = 1e-10
    )
    expect_equal(
      TruncComp2:::jointContrastLRT.parametric(case$data, case$mu, case$log_or),
      runtime_reference$W,
      tolerance = 1e-10
    )
  }
})

test_that("adjusted formula interface stores adjustment metadata and conditional outputs", {
  case <- adjusted_lrt_case(20260521, 18)
  fit <- truncComp(Y ~ R,
                   atom = 0,
                   data = model_data_for_formula(case, adjust = ~ L1 + L2),
                   method = "LRT",
                   adjust = ~ L1 + L2)

  expect_s3_class(fit, "TruncComp2")
  expect_true(fit$success)
  expect_equal(fit$adjust, "L1 + L2")
  expect_true(all(c("L1", "L2") %in% names(fit$data)))
  expect_true(is.na(fit$Delta))
  expect_equal(fit$DeltaCI, c(NA_real_, NA_real_))

  summary_text <- paste(capture.output(summary(fit)), collapse = "\n")
  expect_match(summary_text, "Adjusted for: L1 \\+ L2")

  capture.output(ci <- confint(fit, type = "marginal"))
  expect_equal(unname(ci["Difference in means among the observed:", ]), fit$muDeltaCI)
  expect_equal(unname(ci["Odds ratio of being observed:", ]), fit$alphaDeltaCI)
  expect_error(confint(fit, type = "simultaneous"), "adjusted fits")
  expect_error(jointContrastCI(fit, plot = FALSE), "adjusted fits")
})

test_that("packaged adjusted example illustrates attenuation after covariate adjustment", {
  example_data <- loadTruncComp2AdjustedExample()

  expect_s3_class(example_data, "data.frame")
  expect_equal(names(example_data), c("R", "L", "Y"))
  expect_equal(nrow(example_data), 50)
  expect_equal(as.integer(table(example_data$R)), c(25, 25))
  expect_true(is.factor(example_data$L))
  expect_equal(levels(example_data$L), c("low", "mid", "high"))

  expect_gte(sum(example_data$Y[example_data$R == 0] != 0), 2)
  expect_gte(sum(example_data$Y[example_data$R == 1] != 0), 2)
  expect_false(all(example_data$Y != 0))
  expect_equal(unname(as.integer(table(example_data$L, example_data$R)[, "0"])),
               c(14, 6, 5))
  expect_equal(unname(as.integer(table(example_data$L, example_data$R)[, "1"])),
               c(3, 8, 14))

  fit_unadjusted <- truncComp(Y ~ R,
                              atom = 0,
                              data = example_data[, c("Y", "R")],
                              method = "LRT")
  fit_adjusted <- truncComp(Y ~ R,
                            atom = 0,
                            data = example_data,
                            method = "LRT",
                            adjust = ~ L)

  expect_true(fit_unadjusted$success)
  expect_true(fit_adjusted$success)
  expect_equal(fit_adjusted$adjust, "L")

  expect_lt(fit_unadjusted$p, 0.05)
  expect_gt(fit_adjusted$p, 0.05)
  expect_lt(fit_adjusted$W, fit_unadjusted$W)
  expect_lt(abs(fit_adjusted$muDelta), abs(fit_unadjusted$muDelta))
  expect_lt(abs(log(as.numeric(fit_adjusted$alphaDelta))),
            abs(log(as.numeric(fit_unadjusted$alphaDelta))))

  expect_equal(fit_unadjusted$p, 0.04787936, tolerance = 1e-4)
  expect_equal(fit_adjusted$p, 0.14451384, tolerance = 1e-4)
  expect_equal(fit_unadjusted$W, 6.078142, tolerance = 1e-4)
  expect_equal(fit_adjusted$W, 3.86876, tolerance = 1e-4)
  expect_equal(fit_unadjusted$muDelta, 0.6720578, tolerance = 1e-4)
  expect_equal(fit_adjusted$muDelta, 0.5511571, tolerance = 1e-4)
  expect_equal(as.numeric(fit_unadjusted$alphaDelta), 3.160494, tolerance = 1e-5)
  expect_equal(as.numeric(fit_adjusted$alphaDelta), 1.836725, tolerance = 1e-5)
})

test_that("adjust equals ~1 reproduces the unadjusted parametric fit", {
  case <- adjusted_lrt_case(20260522, 18)

  unadjusted <- truncComp(Y ~ R, atom = 0, data = case[, c("Y", "R")], method = "LRT")
  intercept_only <- truncComp(Y ~ R,
                              atom = 0,
                              data = case[, c("Y", "R")],
                              method = "LRT",
                              adjust = ~ 1)

  expect_true(unadjusted$success)
  expect_true(intercept_only$success)
  expect_equal(intercept_only$adjust, NULL)
  expect_equal(intercept_only$muDelta, unadjusted$muDelta, tolerance = 1e-12)
  expect_equal(intercept_only$alphaDelta, unadjusted$alphaDelta, tolerance = 1e-12)
  expect_equal(intercept_only$W, unadjusted$W, tolerance = 1e-12)
  expect_equal(intercept_only$muDeltaCI, unadjusted$muDeltaCI, tolerance = 1e-12)
  expect_equal(intercept_only$alphaDeltaCI, unadjusted$alphaDeltaCI, tolerance = 1e-12)
})

test_that("adjusted default interface accepts data.frame and matrix covariates", {
  case <- adjusted_lrt_case(20260523, 18)

  fit_df <- truncComp.default(case$Y,
                              case$A,
                              case$R,
                              method = "LRT",
                              adjust = case["L1"])
  fit_matrix <- truncComp.default(case$Y,
                                  case$A,
                                  case$R,
                                  method = "LRT",
                                  adjust = as.matrix(case["L1"]))
  fit_formula <- truncComp(Y ~ R,
                           atom = 0,
                           data = model_data_for_formula(case, adjust = ~ L1),
                           method = "LRT",
                           adjust = ~ L1)

  expect_true(fit_df$success)
  expect_true(fit_matrix$success)
  expect_true(fit_formula$success)
  expect_equal(fit_df$muDelta, fit_formula$muDelta, tolerance = 1e-10)
  expect_equal(fit_df$alphaDelta, fit_formula$alphaDelta, tolerance = 1e-10)
  expect_equal(fit_df$W, fit_formula$W, tolerance = 1e-10)
  expect_equal(fit_matrix$muDelta, fit_formula$muDelta, tolerance = 1e-10)
  expect_equal(fit_matrix$alphaDelta, fit_formula$alphaDelta, tolerance = 1e-10)
  expect_equal(fit_matrix$W, fit_formula$W, tolerance = 1e-10)
})

test_that("covariate adjustment changes the conditional fit under strong confounding", {
  set.seed(20260524)
  n <- 40
  l1 <- c(stats::rnorm(n, -1, 0.25), stats::rnorm(n, 1, 0.25))
  r <- c(rep(0, n), rep(1, n))
  a <- stats::rbinom(2 * n, 1, stats::plogis(-0.2 + 1.2 * l1))
  y <- numeric(2 * n)
  y[a == 1] <- stats::rnorm(sum(a == 1), mean = 2 + 1.4 * l1[a == 1], sd = 0.4)
  case <- data.frame(Y = y, A = a, R = r, L1 = l1)

  expect_gte(sum(a[r == 0]), 3)
  expect_gte(sum(a[r == 1]), 3)
  expect_false(all(a == 1))

  unadjusted <- truncComp(Y ~ R, atom = 0, data = case[, c("Y", "R")], method = "LRT")
  adjusted <- truncComp(Y ~ R,
                        atom = 0,
                        data = case[, c("Y", "R", "L1")],
                        method = "LRT",
                        adjust = ~ L1)

  expect_true(unadjusted$success)
  expect_true(adjusted$success)
  expect_lt(abs(adjusted$muDelta), abs(unadjusted$muDelta))
  expect_lt(adjusted$W, unadjusted$W)
})

test_that("adjusted interface validates formulas and default adjustment shapes", {
  case <- adjusted_lrt_case(20260525, 16)
  formula_data <- model_data_for_formula(case, adjust = ~ L1 + L2)

  expect_error(
    truncComp(Y ~ R, atom = 0, data = formula_data, method = "LRT", adjust = Y ~ L1),
    "one-sided formula"
  )
  expect_error(
    truncComp(Y ~ R, atom = 0, data = formula_data, method = "LRT", adjust = ~ R + L1),
    "must not include the outcome or treatment variable"
  )
  expect_error(
    truncComp(Y ~ R, atom = 0, data = formula_data, method = "LRT", adjust = ~ L1 * L2),
    "additive"
  )
  expect_error(
    truncComp(Y ~ R, atom = 0, data = formula_data, method = "SPLRT", adjust = ~ L1:L2),
    "additive"
  )

  missing_outcome <- formula_data
  missing_outcome$Y[1] <- NA_real_
  expect_error(
    truncComp(Y ~ R, atom = 0, data = missing_outcome, method = "LRT", adjust = ~ L1 + L2),
    "missing values"
  )

  missing_covariate <- formula_data
  missing_covariate$L1[1] <- NA_real_
  expect_error(
    truncComp(Y ~ R, atom = 0, data = missing_covariate, method = "LRT", adjust = ~ L1 + L2),
    "missing values"
  )

  expect_error(
    truncComp.default(case$Y, case$A, case$R, method = "LRT", adjust = case["L1"][1:5, , drop = FALSE]),
    "same number of rows"
  )
})

test_that("observed-only factor levels are dropped before adjusted observed-outcome fits", {
  data <- data.frame(
    Y = c(1.1, 1.3, 0, 0, 2.1, 2.4, 0, 0, 0, 0),
    A = c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0),
    R = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 1),
    L = factor(c("a", "b", "c", "c", "a", "b", "c", "c", "c", "c"),
               levels = c("a", "b", "c"))
  )

  fits <- TruncComp2:::parametric_fit_models(data, adjust = ~ L)

  expect_equal(levels(fits$normal_data$L), c("a", "b"))

  normal_matrix <- stats::model.matrix(fits$formulas$normal_alt, fits$normal_data)
  expect_equal(qr(normal_matrix)$rank, ncol(normal_matrix))

  el_design <- TruncComp2:::el_regression_design(fits$formulas$normal_alt, fits$normal_data)
  expect_true(el_design$success)
})

test_that("adjusted fits fail cleanly when treatment is aliased or the logistic model separates", {
  case <- adjusted_lrt_case(20260526, 18)

  aliased_fit <- truncComp.default(case$Y,
                                   case$A,
                                   case$R,
                                   method = "LRT",
                                   adjust = data.frame(Rcopy = case$R))
  expect_s3_class(aliased_fit, "TruncComp2")
  expect_false(aliased_fit$success)
  expect_match(aliased_fit$error, "not estimable")

  separation_fit <- truncComp.default(case$Y,
                                      case$A,
                                      case$R,
                                      method = "LRT",
                                      adjust = data.frame(Lsep = case$A))
  expect_s3_class(separation_fit, "TruncComp2")
  expect_false(separation_fit$success)
  expect_match(separation_fit$error, "not estimable")
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
