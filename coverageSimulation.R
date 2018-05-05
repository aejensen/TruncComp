rm(list=ls())

library(TruncComp)
library(parallel)
library(sp)

library(devtools)
install_github('aejensen/TruncComp')

simCoverage <- function(nsim, n, scenario, method, ncores, seed=12345) {
  set.seed(seed)

  out <- parallel::mclapply(1:nsim, function(sim) {
    cat(sim, " ")

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
    } else if(scenario == 4) {
      data <- TruncComp:::simTruncData(n, 3, 4, 0.4, 0.3, dist="t-sq")
      muDeltaTrue <- 1
      orDeltaTrue <- (0.3/(1-0.3)) / (0.4/(1-0.4))
    }

    m <- truncComp.default(data$Y, data$A, data$R, method=method)

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
  }, mc.cores=ncores, mc.preschedule = FALSE)

  out
}

getStats <- function(res, muDelta, logOR) {
  mat <- do.call("rbind", res)
  out <- c(biasMuDelta = mean(mat[,1] - muDelta, na.rm=TRUE),
           sdMuDelta = sd(mat[,1], na.rm=TRUE),
           biasLogOR = mean(mat[,2] - logOR),
           sdLogOR = sd(mat[,1]),
           covMuDelta = mean(mat[,3]),
           covLogOR = mean(mat[,4]),
           covJoint = mean(mat[,5]))
  round(out, 4)
}

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
data <- TruncComp:::simTruncData(20000, 3, 4, 1, 1, dist="t-sq")
mean(data$Y[data$R == 1]) - mean(data$Y[data$R == 0])
