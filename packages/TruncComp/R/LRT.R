parametric_safe_xlogy <- function(x, y) {
  out <- numeric(length(x))
  nonzero <- x != 0

  if(any(nonzero)) {
    out[nonzero] <- x[nonzero] * log(y[nonzero])
  }

  out
}

parametric_odds_ratio <- function(pi1, pi0) {
  if((pi1 == 0 && pi0 == 0) || (pi1 == 1 && pi0 == 1)) {
    return(1)
  }

  if(pi1 == 1 || pi0 == 0) {
    return(Inf)
  }

  if(pi1 == 0 || pi0 == 1) {
    return(0)
  }

  (pi1 / (1 - pi1)) / (pi0 / (1 - pi0))
}

parametric_wald_interval <- function(estimate, se, conf.level) {
  if(!(length(se) == 1 && is.finite(se) && se >= 0)) {
    return(c(NA_real_, NA_real_))
  }

  z <- stats::qnorm((1 + conf.level) / 2)
  estimate + c(-1, 1) * z * se
}

parametric_loglik_bernoulli <- function(k, n, p) {
  sum(parametric_safe_xlogy(c(k, n - k), c(p, 1 - p)))
}

parametric_bernoulli_summary <- function(a, r, conf.level, tol = 1e-12) {
  n0 <- sum(r == 0)
  n1 <- sum(r == 1)
  k0 <- sum(a[r == 0])
  k1 <- sum(a[r == 1])

  pi0 <- k0 / n0
  pi1 <- k1 / n1
  pi <- (k0 + k1) / (n0 + n1)

  ell_h0 <- parametric_loglik_bernoulli(k0 + k1, n0 + n1, pi)
  ell_ha <- parametric_loglik_bernoulli(k0, n0, pi0) +
    parametric_loglik_bernoulli(k1, n1, pi1)

  W <- as.numeric(2 * (ell_ha - ell_h0))
  if(is.finite(W) && W < tol) {
    W <- 0
  }

  alphaDelta <- parametric_odds_ratio(pi1, pi0)

  regular_cells <- all(c(k0, n0 - k0, k1, n1 - k1) > 0)
  if(regular_cells) {
    log_or <- stats::qlogis(pi1) - stats::qlogis(pi0)
    se_log_or <- sqrt(1 / k1 + 1 / (n1 - k1) + 1 / k0 + 1 / (n0 - k0))
    alphaDeltaCI <- exp(parametric_wald_interval(log_or, se_log_or, conf.level))
  } else {
    alphaDeltaCI <- c(NA_real_, NA_real_)
  }

  list(
    n0 = n0,
    n1 = n1,
    k0 = k0,
    k1 = k1,
    pi0 = pi0,
    pi1 = pi1,
    pi = pi,
    ell_h0 = ell_h0,
    ell_ha = ell_ha,
    W = W,
    alphaDelta = alphaDelta,
    alphaDeltaCI = alphaDeltaCI
  )
}

parametric_normal_loglik_mle <- function(sample_size, sse, tol = 1e-12) {
  if(sample_size == 0) {
    return(0)
  }

  sse <- max(as.numeric(sse), 0)
  if(sse <= tol) {
    return(Inf)
  }

  -0.5 * sample_size * (log(2 * pi) + 1 + log(sse / sample_size))
}

parametric_normal_summary <- function(y0, y1, conf.level, tol = 1e-12) {
  m0 <- length(y0)
  m1 <- length(y1)
  m <- m0 + m1

  mu0 <- mean(y0)
  mu1 <- mean(y1)
  mu <- mean(c(y0, y1))

  sse1 <- sum((y0 - mu0)^2) + sum((y1 - mu1)^2)
  sse0 <- sum((c(y0, y1) - mu)^2)

  sse1 <- max(as.numeric(sse1), 0)
  sse0 <- max(as.numeric(sse0), 0)

  if(sse0 < sse1 && (sse1 - sse0) <= tol * max(1, sse1)) {
    sse0 <- sse1
  }

  ell_h0 <- parametric_normal_loglik_mle(m, sse0, tol = tol)
  ell_ha <- parametric_normal_loglik_mle(m, sse1, tol = tol)

  if(sse1 <= tol) {
    W <- if(sse0 <= tol) 0 else Inf
    muDeltaCI <- c(NA_real_, NA_real_)
    sigma2_hat <- 0
  } else {
    W <- as.numeric(m * log(sse0 / sse1))
    if(is.finite(W) && W < tol) {
      W <- 0
    }

    sigma2_hat <- sse1 / m
    se_mu_delta <- sqrt(sigma2_hat * (1 / m0 + 1 / m1))
    muDeltaCI <- parametric_wald_interval(mu1 - mu0, se_mu_delta, conf.level)
  }

  list(
    m0 = m0,
    m1 = m1,
    m = m,
    mu0 = mu0,
    mu1 = mu1,
    mu = mu,
    sse0 = sse0,
    sse1 = sse1,
    sigma2_hat = sigma2_hat,
    ell_h0 = ell_h0,
    ell_ha = ell_ha,
    W = W,
    muDelta = mu1 - mu0,
    muDeltaCI = muDeltaCI
  )
}

parametric_lrt_fit <- function(data, conf.level = 0.95, tol = 1e-12) {
  y0 <- data$Y[data$R == 0 & data$A == 1]
  y1 <- data$Y[data$R == 1 & data$A == 1]

  bernoulli <- parametric_bernoulli_summary(data$A, data$R, conf.level, tol = tol)
  normal <- parametric_normal_summary(y0, y1, conf.level, tol = tol)

  W <- as.numeric(bernoulli$W + normal$W)
  if(is.finite(W) && W < tol) {
    W <- 0
  }

  list(
    muDelta = normal$muDelta,
    muDeltaCI = normal$muDeltaCI,
    alphaDelta = bernoulli$alphaDelta,
    alphaDeltaCI = bernoulli$alphaDeltaCI,
    Delta = mean(data$Y[data$R == 1]) - mean(data$Y[data$R == 0]),
    DeltaCI = c(NA_real_, NA_real_),
    W_A = bernoulli$W,
    W_Y = normal$W,
    W = W,
    p = stats::pchisq(W, df = 2, lower.tail = FALSE),
    bernoulli = bernoulli,
    normal = normal
  )
}

LRT <- function(data, init = NULL, conf.level = 0.95) {
  fit <- parametric_lrt_fit(data, conf.level = conf.level)

  newTruncComp(muDelta = fit$muDelta,
               muDeltaCI = fit$muDeltaCI,
               alphaDelta = fit$alphaDelta,
               alphaDeltaCI = fit$alphaDeltaCI,
               Delta = fit$Delta,
               DeltaCI = fit$DeltaCI,
               W = fit$W,
               p = fit$p,
               method = "LRT",
               conf.level = conf.level,
               success = TRUE,
               init = init)
}
