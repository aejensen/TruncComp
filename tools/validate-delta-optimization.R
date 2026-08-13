#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(pkgload))
pkgload::load_all(".", quiet = TRUE)

delta_validation_data <- function() {
  data.frame(
    Y = c(
      0, 0, 1.2, 1.5, 1.8, 2.0, 2.2, 0, 1.7, 2.1,
      0, 2.0, 2.3, 2.5, 2.8, 3.0, 3.2, 0, 2.6, 2.9
    ),
    R = c(rep(0, 10), rep(1, 10))
  )
}

assert_close <- function(actual, expected, tolerance, label) {
  difference <- max(abs(as.numeric(actual) - as.numeric(expected)))
  if(!is.finite(difference) || difference > tolerance) {
    stop(sprintf(
      "%s differed by %.6g, exceeding tolerance %.6g.",
      label,
      difference,
      tolerance
    ))
  }
}

coded_samples <- function(reference) {
  list(
    arm0 = c(rep(reference$atom, reference$n0 - reference$k0), reference$y0),
    arm1 = c(rep(reference$atom, reference$n1 - reference$k1), reference$y1)
  )
}

coded_el_arm_components <- function(z, nonatom, target_mean) {
  fit <- el_one_sample_lambda(z - target_mean)
  if(!isTRUE(fit$feasible)) {
    stop("Coded empirical-likelihood reconstruction was infeasible.")
  }
  weights <- 1 / (length(z) * (1 + fit$lambda * (z - target_mean)))
  probability <- sum(weights[nonatom])
  conditional_mean <- sum(weights[nonatom] * z[nonatom]) / probability
  list(
    probability = probability,
    conditional_mean = conditional_mean,
    weights = weights,
    statistic = fit$statistic
  )
}

splrt_factored_candidate <- function(reference, target_delta) {
  coded <- coded_samples(reference)
  theta_fit <- el_mean_diff_theta(target_delta, coded$arm1, coded$arm0)
  if(!isTRUE(theta_fit$feasible)) {
    stop("Two-sample coded empirical likelihood was infeasible.")
  }
  target1 <- theta_fit$theta
  target0 <- target1 - target_delta
  nonatom0 <- seq_along(coded$arm0) > reference$n0 - reference$k0
  nonatom1 <- seq_along(coded$arm1) > reference$n1 - reference$k1
  arm0 <- coded_el_arm_components(coded$arm0, nonatom0, target0)
  arm1 <- coded_el_arm_components(coded$arm1, nonatom1, target1)
  candidate <- delta_full_candidate(
    reference,
    c(
      arm0$conditional_mean,
      arm1$conditional_mean,
      stats::qlogis(arm0$probability),
      stats::qlogis(arm1$probability)
    )
  )
  list(
    candidate = candidate,
    coded_statistic = theta_fit$statistic,
    coded_delta = target_delta
  )
}

validate_center <- function(fit, reference) {
  center <- delta_full_candidate(reference, reference$center)
  if(!isTRUE(center$feasible)) {
    stop(sprintf("%s full criterion was infeasible at its center.", fit$method))
  }
  assert_close(center$statistic, 0, 1e-8, paste(fit$method, "center statistic"))
  assert_close(center$Delta, fit$delta, 1e-10, paste(fit$method, "center Delta"))
}

validate_parametric_endpoints <- function(fit, reference) {
  set.seed(20260718)
  for(df in c(1, 4)) {
    threshold <- stats::qchisq(fit$conf.level, df)
    lower <- delta_parametric_sublevel_endpoint(
      reference, threshold, direction = "lower"
    )
    upper <- delta_parametric_sublevel_endpoint(
      reference, threshold, direction = "upper"
    )
    if(is.null(lower) || is.null(upper)) {
      stop(sprintf("Parametric df=%d endpoint calculation failed.", df))
    }
    interval <- c(lower$candidate$Delta, upper$candidate$Delta)
    for(solution in list(lower, upper)) {
      if(!isTRUE(solution$candidate$feasible)) {
        stop(sprintf("Parametric df=%d endpoint was infeasible.", df))
      }
      assert_close(
        solution$candidate$statistic,
        threshold,
        2e-6,
        sprintf("parametric df=%d endpoint cutoff", df)
      )
      assert_close(
        solution$candidate$Delta,
        delta_from_components(
          reference$atom,
          solution$candidate$mu0,
          solution$candidate$mu1,
          solution$candidate$p0,
          solution$candidate$p1
        ),
        1e-10,
        sprintf("parametric df=%d endpoint functional", df)
      )
    }
    if(interval[1] > reference$delta_hat || interval[2] < reference$delta_hat) {
      stop(sprintf("Parametric df=%d interval did not contain its center.", df))
    }

    eta_bounds <- list(
      lower = c(
        delta_logit_sublevel_bounds(reference, 0L, threshold)[1],
        delta_logit_sublevel_bounds(reference, 1L, threshold)[1]
      ),
      upper = c(
        delta_logit_sublevel_bounds(reference, 0L, threshold)[2],
        delta_logit_sublevel_bounds(reference, 1L, threshold)[2]
      )
    )
    for(i in seq_len(200)) {
      phi <- stats::runif(1, 0, 2 * pi)
      radial_limit <- delta_parametric_radial_limit(
        reference, threshold, phi, eta_bounds
      )
      radius <- stats::runif(1, 0, radial_limit)
      eta <- c(reference$eta0_hat, reference$eta1_hat) +
        radius * c(cos(phi), sin(phi))
      for(direction in c("lower", "upper")) {
        candidate <- delta_parametric_eta_candidate(
          reference, eta, threshold, direction = direction
        )
        if(!isTRUE(candidate$feasible) ||
           candidate$Delta < interval[1] - 2e-5 ||
           candidate$Delta > interval[2] + 2e-5) {
          stop(sprintf(
            "Parametric df=%d radial endpoint missed an accepted candidate.",
            df
          ))
        }
      }
    }
    cat(sprintf(
      "lrt df=%d interval=[%.8f, %.8f] cutoff residuals=(%.3g, %.3g)\n",
      df,
      interval[1],
      interval[2],
      lower$candidate$statistic - threshold,
      upper$candidate$statistic - threshold
    ))
  }
}

validate_splrt_equivalence <- function(fit, reference) {
  coded <- coded_samples(reference)
  estimate <- mean(coded$arm1) - mean(coded$arm0)
  support <- el_mean_diff_delta_range(coded$arm1, coded$arm0)
  targets <- estimate + c(-0.15, 0.15) * diff(support)

  for(target in targets) {
    reconstructed <- splrt_factored_candidate(reference, target)
    if(!isTRUE(reconstructed$candidate$feasible)) {
      stop("Reconstructed factored SPLRT candidate was infeasible.")
    }
    assert_close(
      reconstructed$candidate$Delta,
      target,
      1e-8,
      "factored versus coded Delta"
    )
    assert_close(
      reconstructed$candidate$statistic,
      reconstructed$coded_statistic,
      2e-7,
      "factored versus coded empirical-likelihood statistic"
    )
  }

  for(df in c(1, 4)) {
    threshold <- stats::qchisq(fit$conf.level, df)
    interval <- delta_splrt_sublevel_interval(reference, threshold)
    endpoint_statistics <- vapply(interval, function(target) {
      el_mean_diff_statistic(coded$arm1, coded$arm0, target)
    }, numeric(1))
    assert_close(
      endpoint_statistics,
      rep(threshold, 2),
      2e-6,
      sprintf("SPLRT df=%d endpoint cutoff", df)
    )
    for(target in interval) {
      reconstructed <- splrt_factored_candidate(reference, target)
      assert_close(
        reconstructed$candidate$statistic,
        threshold,
        2e-6,
        sprintf("SPLRT df=%d factored endpoint cutoff", df)
      )
    }
    if(interval[1] > estimate || interval[2] < estimate) {
      stop(sprintf("SPLRT df=%d interval did not contain its center.", df))
    }
    cat(sprintf(
      "splrt df=%d interval=[%.8f, %.8f] cutoff residuals=(%.3g, %.3g)\n",
      df,
      interval[1],
      interval[2],
      endpoint_statistics[1] - threshold,
      endpoint_statistics[2] - threshold
    ))
  }
}

cat("Validating exact coded scale sublevel calculations\n")
data <- delta_validation_data()
fit_lrt <- trunc_comp(Y ~ R, atom = 0, data = data, method = "lrt")
fit_splrt <- trunc_comp(Y ~ R, atom = 0, data = data, method = "splrt")
for(fit in list(fit_lrt, fit_splrt)) {
  if(!isTRUE(fit$success)) {
    stop(sprintf("%s fit failed during coded interval validation.", fit$method))
  }
  reference <- delta_full_reference(fit)
  if(is.null(reference)) {
    stop(sprintf("%s coded likelihood reference was unavailable.", fit$method))
  }
  validate_center(fit, reference)
  if(identical(fit$method, "lrt")) {
    validate_parametric_endpoints(fit, reference)
  } else {
    validate_splrt_equivalence(fit, reference)
  }
}
cat("Exact coded scale validation complete\n")
