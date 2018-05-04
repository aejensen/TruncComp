library(TruncComp)
library(parallel)

#########################
rm(list=ls())

par(mfrow=c(3,3))
for(k in 1:9) {
  data <- TruncComp:::simTruncData(50, 0, 0.3, 0.85, 0.9)
  m <- truncComp.default(data$Y, data$A, data$R, method="SPLRT")
  a <- jointContrastCI(m)
}

#########################
simCoverage <- function(nsim, n, method, seed=12345) {
  set.seed(seed)
  coverageOut <- t(sapply(1:nsim, function(sim) {
    #cat(sim, "\n")
    data <- TruncComp:::simTruncData(n, 0, 0.3, 0.85, 0.9)
    m <- truncComp.default(data$Y, data$A, data$R, method=method)

    muDeltaTrue <- 0.3
    orDeltaTrue <- (0.9/(1-0.9)) / (0.85/(1-0.85))

    coverage.muDelta <- (m$muDeltaCI[1] < muDeltaTrue) & m$muDeltaCI[2] > muDeltaTrue
    coverage.OR <- (m$alphaDeltaCI[1] < orDeltaTrue) & (m$alphaDeltaCI[2] > orDeltaTrue)

    c(coverage.muDelta, coverage.OR)
  }))
  c(mean(coverageOut[,1] == TRUE, na.rm=TRUE), mean(coverageOut[,2] == TRUE, na.rm=TRUE))
}

simCoverage(nsim = 500, n = 30, method="SPLRT")
simCoverage(nsim = 500, n = 50, method="SPLRT")
simCoverage(nsim = 500, n = 100, method="SPLRT")

#########################
simJointCoverage <- function(nsim, n, ncores) {
  unlist(parallel::mclapply(1:nsim, function(sim) {
    data <- TruncComp:::simTruncData(n, 0, 0.3, 0.85, 0.9)
    muDeltaTrue <- 0.3
    logorDeltaTrue <- log((0.9/(1-0.9)) / (0.85/(1-0.85)))

    m <- truncComp.default(data$Y, data$A, data$R, method="SPLRT")
    a <- jointContrastCI(m, muDelta = seq(-1, 2, length.out = 20), logORdelta = seq(-5, 5, length.out = 20), plot=FALSE)
    cont <- contourLines(a$muDelta, a$logORdelta, a$surface, levels=stats::qchisq(0.95, 2))
    sp::point.in.polygon(muDeltaTrue, logorDeltaTrue, cont[[1]]$x, cont[[1]]$y)
  }, mc.cores=ncores))
}


#########################
deltaSeq <- seq(-1, 3, length.out=50)
plot(deltaSeq, sapply(deltaSeq, function(delta) TruncComp:::logit.LRT(d, delta)), type="l")
abline(h = stats::qchisq(1-0.05, 1))
abline(v=log(m$alphaDeltaCI[1]))
abline(v=log(m$alphaDeltaCI[2]))

TruncComp:::logit.LRT(d, 0)
anova(stats::glm(A ~ 1, family=stats::binomial(), data=d),
      stats::glm(A ~ Z, family=stats::binomial(), data=d), test="LRT")$Deviance[2]

