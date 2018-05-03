#logistic log-likelihood function
logit.likelihood <- function(data, beta) {
  X <- model.matrix(~R, data=data)
  y <- data$A

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

logit.likelihood.profile <- function(data, delta, interval=c(-5,5)) {
  profileOpt <- function(b0) logit.likelihood(data, c(b0, delta))

  opt <- optimize(profileOpt, interval=interval, maximum=FALSE)

  logit.likelihood(data, c(opt$minimum, delta))
}

logit.LRT <- function(data, delta) {
  m1 <- stats::glm(A ~ R, family=stats::binomial(), data=data)
  ll1 <- logit.likelihood(data, coef(m1))

  ll2 <- logit.likelihood.profile(data, delta)

  2*(ll2 - ll1)
}

