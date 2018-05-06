SPLRT <- function(data, conf.level = 0.95) {
  yAlive1 <- data[data$R == 0 & data$A == 1, "Y"]
  yAlive2 <- data[data$R == 1 & data$A == 1, "Y"]

  #Empirical likelihood ratio test
  ELRT <- EL::EL.means(yAlive2, yAlive1, conf.level = conf.level)

  muDelta <- as.numeric(ELRT$estimate)
  muDeltaCI <- as.numeric(ELRT$conf.int)
  muW <- as.numeric(ELRT$statistic)
  #muP <- as.numeric(ELRT$p.value)

  #Fit logistic models
  m0 <- stats::glm(A ~ 1, family=stats::binomial(), data=data)
  m1 <- stats::glm(A ~ R, family=stats::binomial(), data=data)

  binomConfint <- suppressMessages(stats::confint(m1, level = conf.level))
  alphaDelta <- exp(stats::coef(m1)["R"])
  alphaDeltaCI <- as.numeric(exp(binomConfint["R",]))

  binomTest <- stats::anova(m0, m1, test="LRT")
  alphaW <- as.numeric(binomTest$Deviance[2])
  #alphaP <- binomTest$"Pr(>Chi)"[2]
  #W.binom <- 2 * (logLik(m1) - logLik(m0))

  #Get Delta
  delta <- mean(data$Y[data$R == 1]) - mean(data$Y[data$R == 0])
  deltaCI <- c(NA, NA)

  #Joint likelihood ratio test
  W <- as.numeric(muW + alphaW) #Joint test statistic
  p <- 1 - stats::pchisq(W, 2)  #Joint p-value

  out <- list(muDelta = muDelta,
              muDeltaCI = muDeltaCI,
              alphaDelta = alphaDelta,
              alphaDeltaCI = alphaDeltaCI,
              Delta = delta,
              DeltaCI = deltaCI,
              W = W,
              p = p,
              method ="Semi-empirical Likelihood Ratio Test",
              conf.level = conf.level,
              success = TRUE,
              error = "",
              init = NULL)

  class(out) <- append(class(out), "TruncComp")
  out
}

