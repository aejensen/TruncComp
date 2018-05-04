library(TruncComp)
library(parallel)
library(sp)

#########################
rm(list=ls())

par(mfrow=c(3,3))
for(k in 1:9) {
  data <- TruncComp:::simTruncData(50, 0, 0.3, 0.85, 0.9)
  m <- truncComp.default(data$Y, data$A, data$R, method="SPLRT")
  a <- jointContrastCI(m)
}

#########################
rm(list=ls())

library(TruncComp)
library(parallel)
library(sp)

simCoverage <- function(nsim, n, scenario, method, ncores, seed=12345) {
  set.seed(seed)

  out <- parallel::mclapply(1:nsim, function(sim) {
    #cat(sim, "\n")

    if(scenario == 1) {
      data <- TruncComp:::simTruncData(n, 3, 4, 0.35, 0.35)
      muDeltaTrue <- 1
      orDeltaTrue <- (0.35/(1-0.35)) / (0.35/(1-0.35))
    } else if(scenario == 2) {
      data <- TruncComp:::simTruncData(n, 3.5, 3.5, 0.4, 0.3)
      muDeltaTrue <- 0
      orDeltaTrue <- (0.3/(1-0.3)) / (0.4/(1-0.4))
    } else if(scenario == 3) {
      data <- TruncComp:::simTruncData(n, 3, 4, 0.4, 0.3)
      muDeltaTrue <- 1
      orDeltaTrue <- (0.3/(1-0.3)) / (0.4/(1-0.4))
    } else if(scenatio == 4) {
      stop("not implemented yet")
    }

    m <- truncComp.default(data$Y, data$A, data$R, method=method)

    coverage.muDelta <- (m$muDeltaCI[1] < muDeltaTrue) & m$muDeltaCI[2] > muDeltaTrue
    coverage.OR <- (m$alphaDeltaCI[1] < orDeltaTrue) & (m$alphaDeltaCI[2] > orDeltaTrue)

    c(coverage.muDelta, coverage.OR)
  }, mc.cores=ncores)

  out <- do.call("rbind", out)

  cov.muDelta <- mean(out[,1] == TRUE, na.rm=TRUE)
  cov.OR <- mean(out[,2] == TRUE, na.rm=TRUE)

  c(cov.muDelta, cov.OR)
}

cov.1.50 <- simCoverage(nsim = 10^4, n = 50, scenario = 1, method="SPLRT", ncores=64)
cov.2.50 <- simCoverage(nsim = 10^4, n = 50, scenario = 2, method="SPLRT", ncores=64)
cov.3.50 <- simCoverage(nsim = 10^4, n = 50, scenario = 3, method="SPLRT", ncores=64)

cov.1.100 <- simCoverage(nsim = 10^4, n = 100, scenario = 1, method="SPLRT", ncores=64)
cov.2.100 <- simCoverage(nsim = 10^4, n = 100, scenario = 2, method="SPLRT", ncores=64)
cov.3.100 <- simCoverage(nsim = 10^4, n = 100, scenario = 3, method="SPLRT", ncores=64)

cov.1.200 <- simCoverage(nsim = 10^4, n = 200, scenario = 1, method="SPLRT", ncores=64)
cov.2.200 <- simCoverage(nsim = 10^4, n = 200, scenario = 2, method="SPLRT", ncores=64)
cov.3.200 <- simCoverage(nsim = 10^4, n = 200, scenario = 3, method="SPLRT", ncores=64)

cov.1.500 <- simCoverage(nsim = 10^4, n = 500, scenario = 1, method="SPLRT", ncores=64)
cov.2.500 <- simCoverage(nsim = 10^4, n = 500, scenario = 2, method="SPLRT", ncores=64)
cov.3.500 <- simCoverage(nsim = 10^4, n = 500, scenario = 3, method="SPLRT", ncores=64)


cov.1.50
cov.2.50
cov.3.50

cov.1.100
cov.2.100
cov.3.100

cov.1.200
cov.2.200
cov.3.200

cov.1.500
cov.2.500
cov.3.500



#########################
simJointCoverage <- function(nsim, n, ncores, seed=12345) {
  set.seed(seed)
  unlist(parallel::mclapply(1:nsim, function(sim) {
    cat(sim, "\n")
    data <- TruncComp:::simTruncData(n, 0, 0.3, 0.85, 0.9)
    muDeltaTrue <- 0.3
    logorDeltaTrue <- log((0.9/(1-0.9)) / (0.85/(1-0.85)))

    m <- truncComp.default(data$Y, data$A, data$R, method="SPLRT")
    a <- jointContrastCI(m, muDelta = seq(-5, 5, length.out = 100), logORdelta = seq(-10, 10, length.out = 100), plot=FALSE)
    cont <- contourLines(a$muDelta, a$logORdelta, a$surface, levels=stats::qchisq(0.95, 2))
    sp::point.in.polygon(muDeltaTrue, logorDeltaTrue, cont[[1]]$x, cont[[1]]$y)
  }, mc.cores=ncores, mc.preschedule = FALSE))
}

joint1 <- simJointCoverage(nsim = 500, n = 30, ncores=64)
joint2 <- simJointCoverage(nsim = 500, n = 50, ncores=64)
joint3 <- simJointCoverage(nsim = 500, n = 100, ncores=64)


mean(out == 1)

#########################
data <- TruncComp:::simTruncData(30, 0, 0.3, 0.85, 0.9)
m <- truncComp.default(data$Y, data$A, data$R, method="SPLRT")
a <- jointContrastCI(m, muDelta = seq(-5, 5, length.out = 100), logORdelta = seq(-10, 10, length.out = 100), plot=TRUE)

fields::image.plot(a$muDelta, a$logORdelta, a$surface)
contour(a$muDelta, a$logORdelta, a$surface, add=TRUE, levels=stats::qchisq(0.95, 2))

lines(c(m$muDelta, m$muDeltaCI[2]), rep(log(m$alphaDelta), 2), lwd=2)
lines(c(m$muDeltaCI[1], m$muDelta), rep(log(m$alphaDelta), 2), lwd=2)

lines(rep(m$muDelta, 2), c(log(m$alphaDelta), log(m$alphaDeltaCI[2])), lwd=2)
lines(rep(m$muDelta, 2), c(log(m$alphaDeltaCI[1]), log(m$alphaDelta)), lwd=2)

cont <- contourLines(a$muDelta, a$logORdelta, a$surface, levels=stats::qchisq(0.95, 2))
plot(sapply(1:25, function(i) TruncComp:::jointContrastLRT(data, cont[[1]]$x[i], cont[[1]]$y[i])))
abline(h=stats::qchisq(0.95, 2))
length(cont[[1]]$x)

points(m$muDeltaCI[1], log(m$alphaDelta))
points(m$muDeltaCI[2], log(m$alphaDelta))
points(m$muDelta, log(m$alphaDeltaCI[2]))
points(m$muDelta, log(m$alphaDeltaCI[1]))

TruncComp:::jointContrastLRT(data, m$muDeltaCI[1], log(m$alphaDelta))
TruncComp:::jointContrastLRT(data, m$muDeltaCI[2], log(m$alphaDelta))
TruncComp:::jointContrastLRT(data, m$muDelta, log(m$alphaDeltaCI[2]))
TruncComp:::jointContrastLRT(data, m$muDelta, log(m$alphaDeltaCI[1]))

points(m$muDeltaCI[1], log(m$alphaDelta))

deltaSeq <- seq(-1, 1, length.out=20)
plot(deltaSeq, sapply(deltaSeq, function(delta) TruncComp:::jointContrastLRT(data, delta, log(m$alphaDelta))))
abline(h=stats::qchisq(0.95, 2))


TruncComp:::jointContrastLRT(data, cont[[1]]$x[3], cont[[1]]$y[3])
points(cont[[1]]$x[3], cont[[1]]$y[3])

test <- function(data, muDelta) {
  yAlive1 <- data[data$R == 0 & data$A == 1, "Y"]
  yAlive2 <- data[data$R == 1 & data$A == 1, "Y"]
  EL::EL.means(yAlive2, yAlive1, mu = muDelta)$statistic
}

muDeltaSeq <- seq(-0.4, 1.2, length.out=100)
plot(muDeltaSeq, sapply(muDeltaSeq, function(muDelta) test(data, muDelta)), type="l")
abline(h = stats::qchisq(0.95, 1))
abline(v=m$muDeltaCI[1])
abline(v=m$muDeltaCI[2])

logORseq <- seq(-1, 4, length.out=100)
plot(logORseq, sapply(logORseq, function(logOR) TruncComp:::logit.LRT(data, logOR)), type="l")
abline(h = stats::qchisq(0.95, 1))
abline(v=log(m$alphaDeltaCI[1]))
abline(v=log(m$alphaDeltaCI[2]))




#########################
deltaSeq <- seq(-1, 3, length.out=50)
plot(deltaSeq, sapply(deltaSeq, function(delta) TruncComp:::logit.LRT(d, delta)), type="l")
abline(h = stats::qchisq(1-0.05, 1))
abline(v=log(m$alphaDeltaCI[1]))
abline(v=log(m$alphaDeltaCI[2]))

TruncComp:::logit.LRT(d, 0)
anova(stats::glm(A ~ 1, family=stats::binomial(), data=d),
      stats::glm(A ~ Z, family=stats::binomial(), data=d), test="LRT")$Deviance[2]

