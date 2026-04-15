#logistic log-likelihood function
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

  exb <- exp(X%*%beta)
  prob1 <- exb/(1+exb)

  logexb <- log(prob1)

  y0 <- 1 - y
  logexb0 <- log(1 - prob1)

  yt <- t(y)
  y0t <- t(y0)

  logl <- -sum(yt%*%logexb + y0t%*%logexb0)

  return(logl)
}

logit.likelihood <- function(data, beta) {
  logitReference <- list(X = model.matrix(~R, data = data), y = data$A)
  logit.likelihood.prepared(logitReference, beta)
}

logit.likelihood.profile.prepared <- function(logitReference, delta, interval=c(-5,5)) {
  profileOpt <- function(b0) logit.likelihood.prepared(logitReference, c(b0, delta))

  opt <- optimize(profileOpt, interval=interval, maximum=FALSE)

  logit.likelihood.prepared(logitReference, c(opt$minimum, delta))
}

logit.likelihood.profile <- function(data, delta, interval=c(-5,5)) {
  logitReference <- list(X = model.matrix(~R, data = data), y = data$A)
  logit.likelihood.profile.prepared(logitReference, delta, interval = interval)
}

logit.LRT.prepared <- function(logitReference, delta, interval=c(-5,5)) {
  ll2 <- logit.likelihood.profile.prepared(logitReference, delta, interval = interval)

  2 * (ll2 - logitReference$ll1)
}

logit.LRT <- function(data, delta) {
  logitReference <- logit.prepare(data)
  logit.LRT.prepared(logitReference, delta)
}
