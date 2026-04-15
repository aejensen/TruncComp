adjusted_splrt_case <- function(seed, n, adjust = ~ L1 + L2,
                                treatment_effect = 0.5,
                                logit_treatment = 0.35,
                                confounded = FALSE) {
  adjust_vars <- all.vars(adjust)
  min_observed <- max(4, length(adjust_vars) + 2)

  for(offset in 0:2000) {
    set.seed(seed + offset)
    total_n <- 2 * n
    l1 <- stats::rnorm(total_n)
    l2 <- factor(sample(c("a", "b", "c"),
                        total_n,
                        replace = TRUE,
                        prob = c(0.45, 0.35, 0.20)))

    if(confounded) {
      r <- stats::rbinom(total_n, 1, stats::plogis(0.9 * l1 + 0.4 * (l2 == "b") + 0.9 * (l2 == "c")))
      if(sum(r == 0) < n / 2 || sum(r == 1) < n / 2) {
        next
      }
    } else {
      r <- c(rep(0, n), rep(1, n))
    }

    lp_a <- -0.7 + logit_treatment * r + 0.6 * l1 - 0.5 * (l2 == "b") + 0.8 * (l2 == "c")
    a <- stats::rbinom(total_n, 1, stats::plogis(lp_a))

    if(sum(a[r == 0]) < min_observed || sum(a[r == 1]) < min_observed || all(a == 1)) {
      next
    }

    mu <- 1.5 + treatment_effect * r + 0.5 * l1 + 0.4 * (l2 == "b") + 0.9 * (l2 == "c")
    y <- numeric(total_n)
    y[a == 1] <- stats::rnorm(sum(a), mu[a == 1], 0.75)

    data <- data.frame(Y = y, A = a, R = r, L1 = l1, L2 = l2)
    fit <- truncComp.default(data$Y,
                             data$A,
                             data$R,
                             method = "SPLRT",
                             adjust = data[, adjust_vars, drop = FALSE])

    if(isTRUE(fit$success)) {
      return(data)
    }
  }

  stop("Unable to generate a regular adjusted SPLRT case.")
}

adjusted_splrt_reference <- function(data, adjust, conf.level = 0.95) {
  formulas <- TruncComp2:::parametric_model_formulas(adjust)
  adjust_vars <- all.vars(adjust)
  observed <- data[data$A == 1, unique(c("Y", "R", adjust_vars)), drop = FALSE]

  glm_null <- suppressWarnings(stats::glm(formulas$bernoulli_null,
                                          family = stats::binomial(),
                                          data = data))
  glm_alt <- suppressWarnings(stats::glm(formulas$bernoulli_alt,
                                         family = stats::binomial(),
                                         data = data))
  el_fit <- TruncComp2:::el_regression_fit(observed,
                                           formulas$normal_alt,
                                           mu = 0,
                                           conf.level = conf.level)

  z <- stats::qnorm((1 + conf.level) / 2)
  log_or <- stats::coef(glm_alt)[["R"]]
  log_or_se <- sqrt(stats::vcov(glm_alt)["R", "R"])

  list(
    W_alpha = 2 * (as.numeric(stats::logLik(glm_alt)) - as.numeric(stats::logLik(glm_null))),
    alphaDelta = exp(as.numeric(log_or)),
    alphaDeltaCI = exp(as.numeric(log_or) + c(-1, 1) * z * log_or_se),
    mu = el_fit
  )
}

combined_el_objective <- function(par, design, delta) {
  denom <- 1 + as.numeric(par)
  if(any(!is.finite(denom)) || any(denom <= 0)) {
    return(1e12)
  }

  -2 * sum(log(denom))
}

direct_dual_statistic <- function(moments) {
  lambda_fit <- stats::optim(
    par = rep(0, ncol(moments)),
    fn = function(lambda) {
      denom <- 1 + as.numeric(moments %*% lambda)
      if(any(!is.finite(denom)) || any(denom <= 0)) {
        return(1e12)
      }
      -2 * sum(log(denom))
    },
    method = "Nelder-Mead",
    control = list(reltol = 1e-10, maxit = 20000)
  )

  -lambda_fit$value
}

direct_profile_statistic <- function(design, delta) {
  nuisance_count <- ncol(design$X) - 1

  if(nuisance_count == 0) {
    beta <- TruncComp2:::el_regression_beta(delta, numeric(0), design$term_index, ncol(design$X))
    return(direct_dual_statistic(TruncComp2:::el_regression_moments(design, beta)))
  }

  nuisance_start <- TruncComp2:::el_regression_restricted_ls(design, delta)
  fit <- stats::optim(
    nuisance_start,
    fn = function(nuisance) {
      beta <- TruncComp2:::el_regression_beta(delta, nuisance, design$term_index, ncol(design$X))
      moments <- TruncComp2:::el_regression_moments(design, beta)
      direct_dual_statistic(moments)
    },
    method = "Nelder-Mead",
    control = list(reltol = 1e-10, maxit = 20000)
  )

  fit$value
}

combined_el_objective <- function(par, design, delta) {
  nuisance_count <- ncol(design$X) - 1
  nuisance <- if(nuisance_count > 0) par[seq_len(nuisance_count)] else numeric(0)
  lambda <- par[nuisance_count + seq_len(ncol(design$X))]
  beta <- TruncComp2:::el_regression_beta(delta, nuisance, design$term_index, ncol(design$X))
  moments <- TruncComp2:::el_regression_moments(design, beta)
  denom <- 1 + as.numeric(moments %*% lambda)

  if(any(!is.finite(denom)) || any(denom <= 0)) {
    return(1e12)
  }

  2 * sum(log(denom))
}

test_that("adjust equals ~1 reproduces the unadjusted SPLRT fit", {
  example_data <- loadTruncComp2AdjustedExample()

  unadjusted <- truncComp(Y ~ R,
                          atom = 0,
                          data = example_data[, c("Y", "R")],
                          method = "SPLRT")
  intercept_only <- truncComp(Y ~ R,
                              atom = 0,
                              data = example_data,
                              method = "SPLRT",
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

test_that("adjusted SPLRT matches direct logistic and regression-EL references", {
  cases <- list(
    list(data = adjusted_splrt_case(202604151, 16, adjust = ~ L1), adjust = ~ L1),
    list(data = adjusted_splrt_case(202604152, 16, adjust = ~ L2), adjust = ~ L2),
    list(data = adjusted_splrt_case(202604153, 18, adjust = ~ L1 + L2), adjust = ~ L1 + L2)
  )

  for(case in cases) {
    fit <- truncComp(Y ~ R,
                     atom = 0,
                     data = case$data[, c("Y", "R", all.vars(case$adjust)), drop = FALSE],
                     method = "SPLRT",
                     adjust = case$adjust)
    ref <- adjusted_splrt_reference(case$data, case$adjust)

    expect_true(fit$success)
    expect_true(ref$mu$success)
    expect_equal(fit$muDelta, ref$mu$estimate, tolerance = 1e-8)
    expect_equal(fit$muDeltaCI, as.numeric(ref$mu$conf.int), tolerance = 1e-6)
    expect_equal(fit$alphaDelta, ref$alphaDelta, tolerance = 1e-10)
    expect_equal(fit$alphaDeltaCI, ref$alphaDeltaCI, tolerance = 1e-8)
    expect_equal(fit$W, ref$mu$statistic + ref$W_alpha, tolerance = 1e-8)
    expect_equal(fit$p, stats::pchisq(fit$W, df = 2, lower.tail = FALSE), tolerance = 1e-12)
    expect_true(is.na(fit$Delta))
    expect_equal(fit$DeltaCI, c(NA_real_, NA_real_))
  }
})

test_that("regression EL profile matches a direct combined optimization on a small case", {
  observed <- data.frame(
    Y = c(2.1, 1.9, 2.8, 3.0, 2.4, 3.3, 3.1, 4.0),
    R = c(0, 0, 0, 0, 1, 1, 1, 1),
    L = factor(c("a", "b", "a", "b", "a", "b", "a", "b"))
  )

  design <- TruncComp2:::el_regression_design(Y ~ R + L, observed)
  delta <- design$estimate + 0.2
  profile <- TruncComp2:::el_regression_profile_fit(design, delta)

  expect_true(profile$success)
  expect_equal(profile$statistic, direct_profile_statistic(design, delta), tolerance = 1e-4)
})

test_that("adjusted SPLRT handles the packaged adjusted example and attenuates after adjustment", {
  example_data <- loadTruncComp2AdjustedExample()

  fit_unadjusted <- truncComp(Y ~ R,
                              atom = 0,
                              data = example_data[, c("Y", "R")],
                              method = "SPLRT")
  fit_adjusted <- truncComp(Y ~ R,
                            atom = 0,
                            data = example_data,
                            method = "SPLRT",
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

  expect_equal(fit_unadjusted$p, 0.03698883, tolerance = 1e-4)
  expect_equal(fit_adjusted$p, 0.09027938, tolerance = 1e-4)
  expect_equal(fit_unadjusted$W, 6.594279, tolerance = 1e-4)
  expect_equal(fit_adjusted$W, 4.809692, tolerance = 1e-4)
  expect_equal(fit_unadjusted$muDelta, 0.6720578, tolerance = 1e-4)
  expect_equal(fit_adjusted$muDelta, 0.5511571, tolerance = 1e-4)
  expect_equal(as.numeric(fit_unadjusted$alphaDelta), 3.160494, tolerance = 1e-5)
  expect_equal(as.numeric(fit_adjusted$alphaDelta), 1.836725, tolerance = 1e-5)

  summary_text <- paste(capture.output(summary(fit_adjusted)), collapse = "\n")
  expect_match(summary_text, "Adjusted for: L")

  capture.output(ci <- confint(fit_adjusted, type = "marginal"))
  expect_equal(unname(ci["Difference in means among the observed:", ]), fit_adjusted$muDeltaCI)
  expect_equal(unname(ci["Odds ratio of being observed:", ]), fit_adjusted$alphaDeltaCI)
  expect_error(confint(fit_adjusted, type = "simultaneous"),
               "not implemented for adjusted fits")
  expect_error(jointContrastCI(fit_adjusted, plot = FALSE),
               "not implemented for adjusted fits")
})

test_that("adjusted SPLRT fails cleanly on non-regular adjusted fits", {
  case <- adjusted_splrt_case(202604154, 18, adjust = ~ L1 + L2)

  aliased_fit <- truncComp.default(case$Y,
                                   case$A,
                                   case$R,
                                   method = "SPLRT",
                                   adjust = data.frame(Rcopy = case$R))
  expect_s3_class(aliased_fit, "TruncComp2")
  expect_false(aliased_fit$success)
  expect_match(aliased_fit$error, "rank deficient|not estimable")

  separation_fit <- truncComp.default(case$Y,
                                      case$A,
                                      case$R,
                                      method = "SPLRT",
                                      adjust = data.frame(Lsep = case$A))
  expect_s3_class(separation_fit, "TruncComp2")
  expect_false(separation_fit$success)
  expect_match(separation_fit$error, "not estimable")

  sparse_case <- data.frame(
    Y = c(1.2, 1.6, 0, 0, 2.1, 2.4, 0, 0),
    A = c(1, 1, 0, 0, 1, 1, 0, 0),
    R = c(0, 0, 0, 0, 1, 1, 1, 1)
  )
  sparse_adjust <- data.frame(
    L1 = c(-1.0, -0.5, -0.2, -0.1, 0.1, 0.4, 0.6, 0.9),
    L2 = factor(c("a", "b", "a", "b", "a", "b", "a", "b")),
    L3 = c(1, 0, 1, 0, 1, 0, 1, 0)
  )
  sparse_fit <- truncComp.default(sparse_case$Y,
                                  sparse_case$A,
                                  sparse_case$R,
                                  method = "SPLRT",
                                  adjust = sparse_adjust)
  expect_s3_class(sparse_fit, "TruncComp2")
  expect_false(sparse_fit$success)
  expect_match(sparse_fit$error, "too few observed outcomes|not estimable")

  missing_covariate <- case[, c("Y", "R", "L1")]
  missing_covariate$L1[1] <- NA_real_
  expect_error(
    truncComp(Y ~ R, atom = 0, data = missing_covariate, method = "SPLRT", adjust = ~ L1),
    "missing values"
  )
})
