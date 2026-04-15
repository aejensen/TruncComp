#logistic log-likelihood function
logit_log1pexp <- function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

logit_profile_score.prepared <- function(logitReference, beta0, delta) {
  eta <- as.numeric(logitReference$X[, 1] * beta0 + logitReference$X[, 2] * delta)
  sum(stats::plogis(eta) - logitReference$y)
}

logit_profile_interval.prepared <- function(logitReference, delta, start = NULL,
                                            initial_width = 2, max_abs = 100) {
  if(is.null(start) || !is.finite(start)) {
    if(!is.null(logitReference$beta) && length(logitReference$beta) >= 1 &&
       is.finite(logitReference$beta[[1]])) {
      start <- as.numeric(logitReference$beta[[1]])
    } else {
      mean_y <- mean(logitReference$y)
      start <- stats::qlogis(min(max(mean_y, 1e-10), 1 - 1e-10))
    }
  }

  width <- initial_width
  lower <- start - width
  upper <- start + width
  score_lower <- logit_profile_score.prepared(logitReference, lower, delta)
  score_upper <- logit_profile_score.prepared(logitReference, upper, delta)

  while((!is.finite(score_lower) || !is.finite(score_upper) ||
         score_lower > 0 || score_upper < 0) &&
        max(abs(c(lower, upper))) < max_abs) {
    width <- width * 2
    lower <- start - width
    upper <- start + width
    score_lower <- logit_profile_score.prepared(logitReference, lower, delta)
    score_upper <- logit_profile_score.prepared(logitReference, upper, delta)
  }

  c(lower, upper)
}

logit_profile_beta0.prepared <- function(logitReference, delta, interval = NULL) {
  if(is.null(interval)) {
    interval <- logit_profile_interval.prepared(logitReference, delta)
  }

  score_lower <- logit_profile_score.prepared(logitReference, interval[1], delta)
  score_upper <- logit_profile_score.prepared(logitReference, interval[2], delta)

  if(all(is.finite(c(score_lower, score_upper))) &&
     score_lower <= 0 && score_upper >= 0) {
    return(stats::uniroot(
      function(b0) logit_profile_score.prepared(logitReference, b0, delta),
      interval = interval,
      tol = 1e-10
    )$root)
  }

  profileOpt <- function(b0) logit.likelihood.prepared(logitReference, c(b0, delta))
  opt <- optimize(profileOpt, interval = interval, maximum = FALSE)
  opt$minimum
}

logit.prepare <- function(data) {
  X <- model.matrix(~R, data = data)
  y <- data$A
  fit <- stats::glm(A ~ R, family = stats::binomial(), data = data)
  beta <- stats::coef(fit)

  list(X = X,
       y = y,
       fit = fit,
       beta = beta,
       ll1 = logit.likelihood.prepared(list(X = X, y = y), beta))
}

logit.likelihood.prepared <- function(logitReference, beta) {
  X <- logitReference$X
  y <- logitReference$y

  eta <- as.numeric(X %*% beta)
  sum(logit_log1pexp(eta) - y * eta)
}

logit.likelihood <- function(data, beta) {
  logitReference <- list(X = model.matrix(~R, data = data), y = data$A)
  logit.likelihood.prepared(logitReference, beta)
}

logit.likelihood.profile.prepared <- function(logitReference, delta, interval = NULL) {
  beta0 <- logit_profile_beta0.prepared(logitReference, delta, interval = interval)
  logit.likelihood.prepared(logitReference, c(beta0, delta))
}

logit.likelihood.profile <- function(data, delta, interval = NULL) {
  logitReference <- list(X = model.matrix(~R, data = data), y = data$A)
  logit.likelihood.profile.prepared(logitReference, delta, interval = interval)
}

logit.LRT.prepared <- function(logitReference, delta, interval = NULL) {
  ll2 <- logit.likelihood.profile.prepared(logitReference, delta, interval = interval)

  2 * (ll2 - logitReference$ll1)
}

logit.LRT <- function(data, delta) {
  logitReference <- logit.prepare(data)
  logit.LRT.prepared(logitReference, delta)
}
