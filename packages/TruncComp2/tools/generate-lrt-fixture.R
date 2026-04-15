script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if(length(script_arg) > 0) {
  normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath("tools/generate-lrt-fixture.R", winslash = "/", mustWork = FALSE)
}
package_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
output_path <- if(length(args) > 0) args[1] else {
  file.path(package_root, "tests", "testthat", "fixtures", "lrt_reference.rds")
}
output_path <- normalizePath(output_path, winslash = "/", mustWork = FALSE)

old_wd <- getwd()
setwd(package_root)
on.exit(setwd(old_wd), add = TRUE)

lib_dir <- tempfile("lrt-lib-")
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
install.packages("bbmle", repos = "https://cran.r-project.org", lib = lib_dir, quiet = TRUE)

old_lib_paths <- .libPaths(c(lib_dir, .libPaths()))
on.exit(.libPaths(old_lib_paths), add = TRUE)

old_getLRstartingValues <- function(d) {
  logOdds <- function(p) log(p / (1 - p))

  mu0Init <- base::mean(d[d$R == 0 & d$A == 1, "Y"])
  muDeltaInit <- base::mean(d[d$R == 1 & d$A == 1, "Y"]) - mu0Init
  sigmaInit <- stats::sd(d[d$A == 1, "Y"])
  alphaInit <- logOdds(base::mean(d[d$R == 0, "A"]))
  alphaDeltaInit <- logOdds(base::mean(d[d$R == 1, "A"])) - alphaInit

  list(mu = mu0Init, muDelta = muDeltaInit,
       alpha = alphaInit, alphaDelta = alphaDeltaInit,
       sigma = sigmaInit)
}

old_lrt_fit <- function(data, conf.level = 0.95) {
  a <- data$A
  y <- data$Y
  z <- data$R

  expit <- function(x) exp(x) / (exp(x) + 1)

  dQoL <- function(A, Y, pA, eta, sigma) {
    stats::dnorm(Y, eta, sigma)^A * pA^A * (1 - pA)^(1 - A)
  }

  logLikH0 <- function(alpha, mu, logSigma) {
    -sum(log(dQoL(a, y, expit(alpha), mu, exp(logSigma))))
  }

  logLikHA <- function(alpha, alphaDelta, mu, muDelta, logSigma) {
    piA <- expit(alpha + alphaDelta * z)
    eta <- mu + muDelta * z
    sigma <- exp(logSigma)
    -sum(log(dQoL(a, y, piA, eta, sigma)))
  }

  init <- old_getLRstartingValues(data)
  initH0 <- list(alpha = init$alpha, mu = init$mu, logSigma = log(init$sigma))
  mH0 <- bbmle::mle2(logLikH0, start = initH0)

  initHA <- list(alpha = init$alpha, alphaDelta = init$alphaDelta, mu = init$mu,
                 muDelta = init$muDelta, logSigma = log(init$sigma))
  mHA <- bbmle::mle2(logLikHA, start = initHA)

  coefHA <- bbmle::coef(mHA)
  confintHA <- bbmle::confint(mHA, level = conf.level, method = "quad", quietly = TRUE)
  lrt <- bbmle::anova(mH0, mHA)

  list(
    muDelta = as.numeric(coefHA["muDelta"]),
    muDeltaCI = as.numeric(confintHA["muDelta", ]),
    alphaDelta = as.numeric(exp(coefHA["alphaDelta"])),
    alphaDeltaCI = as.numeric(exp(confintHA["alphaDelta", ])),
    W = as.numeric(lrt[2, "Chisq"]),
    p = as.numeric(lrt[2, "Pr(>Chisq)"])
  )
}

simulate_case <- function(seed, n, mu0, mu1, sigma, pi0, pi1) {
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

    return(data.frame(
      Y = c(y0, y1),
      A = c(a0, a1),
      R = c(rep(0, n), rep(1, n))
    ))
  }

  stop("Unable to generate a valid interior LRT fixture case.")
}

load("data/TruncComp2Example.RData")

cases <- list(
  no_effect = simulate_case(20260401, 24, 2.5, 2.5, 1.0, 0.55, 0.55),
  mean_effect_only = simulate_case(20260402, 24, 2.5, 3.3, 1.0, 0.60, 0.60),
  observation_effect_only = simulate_case(20260403, 24, 2.5, 2.5, 1.0, 0.40, 0.75),
  both_effects = simulate_case(20260404, 24, 2.3, 3.2, 0.8, 0.45, 0.70),
  small_sample = simulate_case(20260405, 8, 1.8, 2.6, 1.1, 0.50, 0.65),
  low_variance = simulate_case(20260406, 20, 3.0, 3.25, 0.2, 0.60, 0.55),
  unequal_observation = simulate_case(20260407, 28, 2.0, 2.8, 0.9, 0.30, 0.80),
  trunccomp_example = transform(TruncComp2Example, A = as.numeric(Y != 0))
)

fixture <- list(
  generated_on = as.Date("2026-04-15"),
  reference_implementation = "optimizer-based parametric LRT via bbmle",
  reference_package = "bbmle",
  reference_version = as.character(utils::packageVersion("bbmle")),
  conf.level = 0.95,
  cases = lapply(names(cases), function(name) {
    standardized <- cases[[name]][, c("Y", "A", "R")]
    reference <- suppressWarnings(old_lrt_fit(standardized, conf.level = 0.95))

    list(
      name = name,
      data = standardized[, c("Y", "R")],
      reference = reference
    )
  })
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(fixture, output_path, version = 2)
