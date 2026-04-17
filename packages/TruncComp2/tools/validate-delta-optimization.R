#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pkgload)
})

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

validate_delta_optimizer <- function(method, tolerance = 0.2) {
  data <- delta_validation_data()
  fit <- truncComp(Y ~ R, atom = 0, data = data, method = method)
  if(!isTRUE(fit$success)) {
    stop(sprintf("%s fit failed unexpectedly during delta optimizer validation.", method))
  }

  surface_data <- delta_surface_for_inference(fit, conf.level = fit$conf.level, resolution = 41)
  grid_projected <- delta_interval_from_surface(surface_data, stats::qchisq(fit$conf.level, 2))
  grid_profile <- delta_interval_from_surface(surface_data, stats::qchisq(fit$conf.level, 1))

  projected_time <- system.time(
    projected_opt <- confint(
      fit,
      parameter = "Delta",
      method = "projected",
      algorithm = "optimize",
      plot = FALSE
    )
  )[["elapsed"]]
  profile_time <- system.time(
    profile_opt <- confint(
      fit,
      parameter = "Delta",
      method = "profile",
      algorithm = "optimize",
      plot = FALSE
    )
  )[["elapsed"]]

  projected_opt <- unname(projected_opt["Delta (projected)", ])
  profile_opt <- unname(profile_opt["Delta (profile)", ])

  if(max(abs(projected_opt - grid_projected)) > tolerance) {
    stop(sprintf(
      "%s projected optimizer interval differs from grid reference by more than %.3f.",
      method,
      tolerance
    ))
  }
  if(max(abs(profile_opt - grid_profile)) > tolerance) {
    stop(sprintf(
      "%s profile optimizer interval differs from grid reference by more than %.3f.",
      method,
      tolerance
    ))
  }

  cat(sprintf(
    "%s: projected optimize %.2fs, profile optimize %.2fs, projected diff %.4f, profile diff %.4f\n",
    method,
    projected_time,
    profile_time,
    max(abs(projected_opt - grid_projected)),
    max(abs(profile_opt - grid_profile))
  ))
}

cat("Validating optimizer-backed Delta intervals against grid references\n")
for(method in c("LRT", "SPLRT")) {
  validate_delta_optimizer(method)
}
cat("Delta optimizer validation complete\n")
