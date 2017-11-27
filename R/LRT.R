getLRstartingValues <- function(d) {
  logOdds <- function(p) log(p / (1-p))

  mu0Init <- base::mean(d[d$Z == 0 & d$A == 1, "Y"])
  muDeltaInit <- base::mean(d[d$Z == 1 & d$A == 1, "Y"])  - mu0Init
  sigmaInit <- stats::sd(d[d$A == 1, "Y"])
  alphaInit <- logOdds(base::mean(d[d$Z == 0, "A"]))
  alphaDeltaInit <- logOdds(base::mean(d[d$Z == 1, "A"]))- alphaInit

  list(mu = mu0Init, muDelta = muDeltaInit,
       alpha = alphaInit, alphaDelta = alphaDeltaInit,
       sigma = sigmaInit)
}


LRT <- function(data, init, conf.level = 0.95) {
  a <- data$A
  y <- data$Y
  z <- data$Z

  #A = 1(y_i is observed)
  expit <- function(x) exp(x)/(exp(x) + 1)

  dQoL <- function(A, Y, pA, eta, sigma) {
    stats::dnorm(Y, eta, sigma)^A * pA^A * (1 - pA)^(1 - A)
  }

  logLikH0 <- function(alpha, mu, logSigma) {
    #The log-likelihood under H_0: no treatment effect
    -sum(log(dQoL(a, y, expit(alpha), mu, exp(logSigma))))
  }

  logLikHA <- function(alpha, alphaDelta, mu, muDelta, logSigma) {
    #The log-likelihood under H_A: effects of treatment
    piA <- expit(alpha + alphaDelta * z)
    eta <- mu + muDelta * z
    sigma <- exp(logSigma)
    -sum(log(dQoL(a, y, piA, eta, sigma)))
  }

  #First, fit the model under H_0
  initH0 <- list(alpha = init$alpha, mu = init$mu, logSigma = log(init$sigma))
  mH0 <-  tryCatch(bbmle::mle2(logLikH0, start = initH0),
                  error = function(e) { NULL })
  if(is.null(mH0)) {
    return(NULL) #Return immediately on error
  }

  #Then fit the model under H_A
  initHA <- list(alpha = init$alpha, alphaDelta = init$alphaDelta, mu = init$mu,
                 muDelta = init$muDelta, logSigma = log(init$sigma))
  mHA <- tryCatch(bbmle::mle2(logLikHA, start = initHA),
                  error = function(e) { NULL })
  if(is.null(mHA)) {
    return(NULL)  #Return immediately on error
  }

  #We have arrived here without errors.
  #Get difference contrasts
  coefHA <- bbmle::coef(mHA)
  muDelta <- as.numeric(coefHA["muDelta"])
  alphaDelta <- as.numeric(exp(coefHA["alphaDelta"]))

  #Get confidence intervals for differences
  confintHA <- bbmle::confint(mHA, level = conf.level, method = "quad", quietly = TRUE)
  muDeltaCI <- as.numeric(confintHA["muDelta", ])
  alphaDeltaCI <- as.numeric(exp(confintHA["alphaDelta",]))

  #LR comparison between H_0 and H_A
  LRT <- bbmle::anova(mH0, mHA)
  W <- LRT[2,"Chisq"]
  p <- LRT[2, "Pr(>Chisq)"]

  out <- list(muDelta = muDelta,
              muDeltaCI = muDeltaCI,
              alphaDelta = alphaDelta,
              alphaDeltaCI = alphaDeltaCI,
              W = W,
              p = p,
              method = "Parametric Likelihood Ratio Test",
              conf.level = conf.level,
              success = TRUE,
              error = "",
              init = init)
  class(out) <- append(class(out), "TruncComp")
  out
}

