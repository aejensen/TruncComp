rm(list=ls())

library(TruncComp)
library(parallel)
library(sp)

library(devtools)
install_github('aejensen/TruncComp')

simCoverage <- function(nsim, n, scenario, method, ncores, seed=12345) {
  set.seed(seed)

  out <- parallel::mclapply(1:nsim, function(sim) {
    #cat(sim, " ")

    if(scenario == 1) {
      f0 <- function(n) stats::rnorm(n, 3, 1)
      f1 <- function(n) stats::rnorm(n, 3.5, 1)
      pi0 <- 0.35
      pi1 <- 0.35
      muDeltaTrue <- 0.5
      orDeltaTrue <- (pi1/(1-pi1)) / (pi0/(1-pi0))
    } else if(scenario == 2) {
      f0 <- function(n) stats::rnorm(n, 3.5, 1)
      f1 <- function(n) stats::rnorm(n, 3.5, 1)
      pi0 <- 0.4
      pi1 <- 0.3
      muDeltaTrue <- 0
      orDeltaTrue <- (pi1/(1-pi1)) / (pi0/(1-pi0))
    } else if(scenario == 3) {
      f0 <- function(n) stats::rnorm(n, 3, 1)
      f1 <- function(n) stats::rnorm(n, 3.5, 1)
      pi0 <- 0.4
      pi1 <- 0.3
      muDeltaTrue <- 0.5
      orDeltaTrue <- (pi1/(1-pi1)) / (pi0/(1-pi0))
    } else if(scenario == 4) {
      f0 <- function(n) stats::rnorm(n, 5, 2)
      f1 <- function(n) rbeta(n, 0.5, 0.5) + 5
      pi0 <- 0.4
      pi1 <- 0.3
      orDeltaTrue <- (pi1/(1-pi1)) / (pi0/(1-pi0))
    }

    #Simulate data
    data <- TruncComp::simulateTruncatedData(n, f0, f1, pi0, pi1)

    #Fit model
    m <- TruncComp::truncComp.default(data$Y, data$A, data$R, method=method)


    #Get marginal inclusion indicators
    inclusion.muDelta <- (m$muDeltaCI[1] < muDeltaTrue) & m$muDeltaCI[2] > muDeltaTrue
    inclusion.OR <- (m$alphaDeltaCI[1] < orDeltaTrue) & (m$alphaDeltaCI[2] > orDeltaTrue)

    #Get joint surface (slow!)
    #jointSurface <- jointContrastCI(m, muDelta = seq(m$muDeltaCI[1] -0.4 , m$muDeltaCI[2] + 0.4, length.out = 40),
    #                                   logORdelta = seq(log(m$alphaDeltaCI[1]) - 0.4, log(m$alphaDeltaCI[2]) + 0.4, length.out = 40),
    #                                   plot=TRUE)

    #Get joint inclusion indicator
    #jointContour <- contourLines(jointSurface$muDelta, jointSurface$logORdelta, jointSurface$surface, levels=stats::qchisq(0.95, 2))
    #inclusion.joint <- sp::point.in.polygon(muDeltaTrue, log(orDeltaTrue), jointContour[[1]]$x, jointContour[[1]]$y)
    inclusion.joint <- 1

    c(m$muDelta, log(m$alphaDelta), inclusion.muDelta, inclusion.OR, as.numeric(inclusion.joint == 1))
  }, mc.cores=ncores, mc.preschedule = TRUE)

  out
}

nSim <- 10^4

#Setup 1
LRT.1.50 <- simCoverage(nsim = nSim, n = 50, scenario = 1, method="LRT", ncores=64)
SPLRT.1.50 <- simCoverage(nsim = nSim, n = 50, scenario = 1, method="SPLRT", ncores=64)

LRT.1.100 <- simCoverage(nsim = nSim, n = 100, scenario = 1, method="LRT", ncores=64)
SPLRT.1.100 <- simCoverage(nsim = nSim, n = 100, scenario = 1, method="SPLRT", ncores=64)

LRT.1.200 <- simCoverage(nsim = nSim, n = 200, scenario = 1, method="LRT", ncores=64)
SPLRT.1.200 <- simCoverage(nsim = nSim, n = 200, scenario = 1, method="SPLRT", ncores=64)

#Setup 2
LRT.2.50 <- simCoverage(nsim = nSim, n = 50, scenario = 2, method="LRT", ncores=64)
SPLRT.2.50 <- simCoverage(nsim = nSim, n = 50, scenario = 2, method="SPLRT", ncores=64)

LRT.2.100 <- simCoverage(nsim = nSim, n = 100, scenario = 2, method="LRT", ncores=64)
SPLRT.2.100 <- simCoverage(nsim = nSim, n = 100, scenario = 2, method="SPLRT", ncores=64)

LRT.2.200 <- simCoverage(nsim = nSim, n = 200, scenario = 2, method="LRT", ncores=64)
SPLRT.2.200 <- simCoverage(nsim = nSim, n = 200, scenario = 2, method="SPLRT", ncores=64)

#Setup 3
LRT.3.50 <- simCoverage(nsim = nSim, n = 50, scenario = 3, method="LRT", ncores=64)
SPLRT.3.50 <- simCoverage(nsim = nSim, n = 50, scenario = 3, method="SPLRT", ncores=64)

LRT.3.100 <- simCoverage(nsim = nSim, n = 100, scenario = 3, method="LRT", ncores=64)
SPLRT.3.100 <- simCoverage(nsim = nSim, n = 100, scenario = 3, method="SPLRT", ncores=64)

LRT.3.200 <- simCoverage(nsim = nSim, n = 200, scenario = 3, method="LRT", ncores=64)
SPLRT.3.200 <- simCoverage(nsim = nSim, n = 200, scenario = 3, method="SPLRT", ncores=64)

#Setup 4
LRT.4.50 <- simCoverage(nsim = nSim, n = 50, scenario = 4, method="LRT", ncores=64)
SPLRT.4.50 <- simCoverage(nsim = nSim, n = 50, scenario = 4, method="SPLRT", ncores=64)

LRT.4.100 <- simCoverage(nsim = nSim, n = 100, scenario = 4, method="LRT", ncores=64)
SPLRT.4.100 <- simCoverage(nsim = nSim, n = 100, scenario = 4, method="SPLRT", ncores=64)

LRT.4.200 <- simCoverage(nsim = nSim, n = 200, scenario = 4, method="LRT", ncores=64)
SPLRT.4.200 <- simCoverage(nsim = nSim, n = 200, scenario = 4, method="SPLRT", ncores=64)






getStats <- function(res, muDelta, logOR) {
  mat <- do.call("rbind", res)
  out <- c(biasMuDelta = mean(mat[,1] - muDelta, na.rm=TRUE),
           sdMuDelta = sd(mat[,1], na.rm=TRUE),
           biasLogOR = mean(mat[,2] - logOR),
           sdLogOR = sd(mat[,1]),
           covMuDelta = mean(mat[,3])*100,
           covLogOR = mean(mat[,4])*100,
           covJoint = mean(mat[,5])*100)
  round(out, 4)
}

getStats(LRT.1.50, 0.5, log(1))
getStats(SPLRT.1.50, 0.5, log(1))
getStats(LRT.1.100, 0.5, log(1))
getStats(SPLRT.1.100, 0.5, log(1))
getStats(LRT.1.200, 0.5, log(1))
getStats(SPLRT.1.200, 0.5, log(1))









cov.1.SPLRT.50 <- simCoverage(nsim = nSim, n = 250, scenario = 1, method="SPLRT", ncores=64)

cov.2.SPLRT.50 <- simCoverage(nsim = nSim, n = 250, scenario = 1, method="SPLRT", ncores=64)
cov.3.SPLRT.50 <- simCoverage(nsim = nSim, n = 250, scenario = 1, method="SPLRT", ncores=64)
cov.4.SPLRT.50 <- simCoverage(nsim = nSim, n = 250, scenario = 1, method="SPLRT", ncores=64)

getStats(cov.1.SPLRT.50, 0.5, log(1))

cov.1.LRT.50 <- simCoverage(nsim = nSim, n = 250, scenario = 1, method="LRT", ncores=64)
getStats(cov.1.LRT.50, 0.5, log(1))


nSim <- 1000

cov.4.SPLRT.200 <- simCoverage(nsim = nSim, n = 2000, scenario = 4, method="SPLRT", ncores=64)
getStats(cov.4.SPLRT.200, 0, log((0.3/(1-0.3)) / (0.4/(1-0.4))))

#nSim <- 2000

cov.1.SPLRT.50 <- simCoverage(nsim = nSim, n = 50, scenario = 1, method="SPLRT", ncores=64)
cov.2.SPLRT.50 <- simCoverage(nsim = nSim, n = 50, scenario = 2, method="SPLRT", ncores=64)
cov.3.SPLRT.50 <- simCoverage(nsim = nSim, n = 50, scenario = 3, method="SPLRT", ncores=64)
cov.4.SPLRT.50 <- simCoverage(nsim = nSim, n = 50, scenario = 4, method="SPLRT", ncores=64)

cov.1.SPLRT.100 <- simCoverage(nsim = nSim, n = 100, scenario = 1, method="SPLRT", ncores=64)
cov.2.SPLRT.100 <- simCoverage(nsim = nSim, n = 100, scenario = 2, method="SPLRT", ncores=64)
cov.3.SPLRT.100 <- simCoverage(nsim = nSim, n = 100, scenario = 3, method="SPLRT", ncores=64)

cov.1.SPLRT.200 <- simCoverage(nsim = nSim, n = 200, scenario = 1, method="SPLRT", ncores=64)
cov.2.SPLRT.200 <- simCoverage(nsim = nSim, n = 200, scenario = 2, method="SPLRT", ncores=64)
cov.3.SPLRT.200 <- simCoverage(nsim = nSim, n = 200, scenario = 3, method="SPLRT", ncores=64)
cov.4.SPLRT.200 <- simCoverage(nsim = nSim, n = 200, scenario = 4, method="SPLRT", ncores=64)



getStats(cov.1.SPLRT.50, 1, log(1))
getStats(cov.2.SPLRT.50, 0, log((0.3/(1-0.3)) / (0.4/(1-0.4))))
getStats(cov.3.SPLRT.50, 1, log((0.3/(1-0.3)) / (0.4/(1-0.4))))

getStats(cov.4.SPLRT.50, 1, log((0.3/(1-0.3)) / (0.4/(1-0.4))))

getStats(cov.1.SPLRT.100, 1, log(1))
getStats(cov.2.SPLRT.100, 0, log((0.3/(1-0.3)) / (0.4/(1-0.4))))
getStats(cov.3.SPLRT.100, 1, log((0.3/(1-0.3)) / (0.4/(1-0.4))))

getStats(cov.1.SPLRT.200, 1, log(1))
getStats(cov.2.SPLRT.200, 0, log((0.3/(1-0.3)) / (0.4/(1-0.4))))
getStats(cov.3.SPLRT.200, 1, log((0.3/(1-0.3)) / (0.4/(1-0.4))))
getStats(cov.4.SPLRT.200, 0, log((0.3/(1-0.3)) / (0.4/(1-0.4)))


save.image("coverageResults.RData")


####
test <- mclapply(1:1000, function(nsim) {
  f0 <- function() gamlss.dist::rSN1(1, mu = 3, sigma=1, nu=12)
  f1 <- function() gamlss.dist::rSN1(1, mu = 4, sigma=1, nu=12)

  d <- TruncComp:::simTruncData2(200, f0, f1, 0.3, 0.4)
  m <- truncComp.default(data$Y, data$A, data$R, method="SPLRT")
  m$muDelta
}, mc.cores=64, mc.preschedule=FALSE)
mean(unlist(test), na.rm=TRUE)


hist(d$Y[d$R == 1])



