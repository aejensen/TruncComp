#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(pkgload))
pkgload::load_all(".", quiet = TRUE)

validation_integer <- function(name, default, minimum) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if(is.na(value) || value < minimum) default else value
}

simulate_two_part <- function(n_per_arm, mu0, mu1, p0, p1, atom,
                              distribution = c("normal", "skewed")) {
  distribution <- match.arg(distribution)
  treatment <- rep(0:1, each = n_per_arm)
  probability <- ifelse(treatment == 0, p0, p1)
  nonatom_indicator <- stats::rbinom(length(treatment), 1, probability)
  mean_value <- ifelse(treatment == 0, mu0, mu1)
  nonatom <- if(distribution == "normal") {
    stats::rnorm(length(treatment), mean_value, 1)
  } else {
    mean_value + stats::rexp(length(treatment)) - 1
  }
  data.frame(
    Y = ifelse(nonatom_indicator == 1, nonatom, atom),
    R = treatment
  )
}

parametric_profile_interval <- function(reference, threshold,
                                        angle_resolution,
                                        radial_resolution) {
  lower <- delta_parametric_sublevel_endpoint(
    reference,
    threshold,
    direction = "lower",
    angle_resolution = angle_resolution,
    radial_resolution = radial_resolution
  )
  upper <- delta_parametric_sublevel_endpoint(
    reference,
    threshold,
    direction = "upper",
    angle_resolution = angle_resolution,
    radial_resolution = radial_resolution
  )
  if(is.null(lower) || is.null(upper)) return(c(NA_real_, NA_real_))
  c(lower$candidate$Delta, upper$candidate$Delta)
}

validate_scenario <- function(method, scenario, repetitions, n_per_arm,
                              distribution, angle_resolution,
                              radial_resolution) {
  profile_covered <- rep(NA, repetitions)
  full_covered <- rep(NA, repetitions)
  full_statistics <- rep(NA_real_, repetitions)
  d_true <- delta_from_components(
    scenario$atom,
    scenario$mu0,
    scenario$mu1,
    scenario$p0,
    scenario$p1
  )
  true_parameter <- c(
    scenario$mu0,
    scenario$mu1,
    stats::qlogis(scenario$p0),
    stats::qlogis(scenario$p1)
  )
  profile_threshold <- stats::qchisq(0.95, 1)
  full_threshold <- stats::qchisq(0.95, 4)

  for(replication in seq_len(repetitions)) {
    data <- simulate_two_part(
      n_per_arm,
      scenario$mu0,
      scenario$mu1,
      scenario$p0,
      scenario$p1,
      scenario$atom,
      distribution = distribution
    )
    fit <- trunc_comp(Y ~ R, atom = scenario$atom, data = data, method = method)
    if(!isTRUE(fit$success)) next
    reference <- delta_full_reference(fit)
    if(is.null(reference)) next

    full_candidate <- delta_full_candidate(reference, true_parameter)
    if(!isTRUE(full_candidate$feasible)) next
    full_statistics[replication] <- full_candidate$statistic
    full_covered[replication] <- full_candidate$statistic <= full_threshold

    interval <- if(identical(method, "lrt")) {
      parametric_profile_interval(
        reference,
        profile_threshold,
        angle_resolution,
        radial_resolution
      )
    } else {
      delta_splrt_sublevel_interval(reference, profile_threshold)
    }
    if(all(is.finite(interval))) {
      profile_covered[replication] <-
        interval[1] <= d_true && d_true <= interval[2]
    }
  }

  keep <- !is.na(profile_covered) & !is.na(full_covered) &
    is.finite(full_statistics)
  if(sum(keep) < 0.9 * repetitions) {
    stop(sprintf("%s %s produced too few regular validation fits.", method, scenario$name))
  }
  profile_coverage <- mean(profile_covered[keep])
  full_coverage <- mean(full_covered[keep])
  coverage_floor <- max(
    0.75,
    0.95 - 4 * sqrt(0.95 * 0.05 / sum(keep))
  )
  if(profile_coverage < coverage_floor || full_coverage < coverage_floor) {
    stop(sprintf(
      "%s %s showed gross asymptotic undercoverage: profile %.3f, full %.3f.",
      method,
      scenario$name,
      profile_coverage,
      full_coverage
    ))
  }

  cat(sprintf(
    "%s %-20s n=%d reps=%d profile=%.3f full=%.3f full-q95=%.3f\n",
    method,
    scenario$name,
    2 * n_per_arm,
    sum(keep),
    profile_coverage,
    full_coverage,
    unname(stats::quantile(full_statistics[keep], 0.95))
  ))
}

scenarios <- list(
  list(name = "joint null", mu0 = 2, mu1 = 2, p0 = 0.5, p1 = 0.5, atom = 0),
  list(name = "mean effect", mu0 = 2, mu1 = 3, p0 = 0.5, p1 = 0.5, atom = 0),
  list(name = "atom effect", mu0 = 2, mu1 = 2, p0 = 0.3, p1 = 0.7, atom = 0),
  list(name = "coded cancellation", mu0 = 4, mu1 = 1, p0 = 0.2, p1 = 0.8, atom = 0)
)

set.seed(20260718)
parametric_repetitions <- validation_integer("TRUNCCOMP_DELTA_PARAM_REPS", 40L, 20L)
empirical_repetitions <- validation_integer("TRUNCCOMP_DELTA_EL_REPS", 30L, 20L)
angle_resolution <- validation_integer("TRUNCCOMP_DELTA_ANGLE_RESOLUTION", 360L, 90L)
radial_resolution <- validation_integer("TRUNCCOMP_DELTA_RADIAL_RESOLUTION", 33L, 17L)

cat("Checking coded scale likelihood calibration\n")
for(scenario in scenarios) {
  validate_scenario(
    "lrt", scenario, parametric_repetitions, 500, "normal",
    angle_resolution, radial_resolution
  )
}
for(scenario in scenarios) {
  validate_scenario(
    "splrt", scenario, empirical_repetitions, 500, "skewed",
    angle_resolution, radial_resolution
  )
}
cat("Coded scale asymptotic validation complete\n")
