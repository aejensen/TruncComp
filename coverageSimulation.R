rm(list=ls())

library(TruncComp)
library(parallel)
library(sp)

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
      stop("not implemented yet")
    }

    m <- truncComp.default(data$Y, data$A, data$R, method=method)

    #Get marginal inclusion indicators
    inclusion.muDelta <- (m$muDeltaCI[1] < muDeltaTrue) & m$muDeltaCI[2] > muDeltaTrue
    inclusion.OR <- (m$alphaDeltaCI[1] < orDeltaTrue) & (m$alphaDeltaCI[2] > orDeltaTrue)

    #Get joint surface (slow!)
    jointSurface <- jointContrastCI(m, muDelta = seq(m$muDeltaCI[1] -0.4 , m$muDeltaCI[2] + 0.4, length.out = 40),
                                       logORdelta = seq(log(m$alphaDeltaCI[1]) - 0.4, log(m$alphaDeltaCI[2]) + 0.4, length.out = 40),
                                       plot=TRUE)

    #Get joint inclusion indicator
    jointContour <- contourLines(jointSurface$muDelta, jointSurface$logORdelta, jointSurface$surface, levels=stats::qchisq(0.95, 2))
    inclusion.joint <- sp::point.in.polygon(muDeltaTrue, log(orDeltaTrue), jointContour[[1]]$x, jointContour[[1]]$y)

    c(inclusion.muDelta, inclusion.OR, as.numeric(inclusion.joint == 1))
  }, mc.cores=ncores, mc.preschedule = FALSE)

  out <- do.call("rbind", out)

  cov.muDelta <- mean(out[,1] == 1, na.rm=TRUE)
  cov.OR <- mean(out[,2] == 1, na.rm=TRUE)
  cov.simult <- mean(out[,3] == 1, na.rm=TRUE)

  c(cov.muDelta, cov.OR, cov.simult)
}

nSim <- 2000

cov.1.SPLRT.50 <- simCoverage(nsim = nSim, n = 50, scenario = 1, method="SPLRT", ncores=64)
cov.2.SPLRT.50 <- simCoverage(nsim = nSim, n = 50, scenario = 2, method="SPLRT", ncores=64)
cov.3.SPLRT.50 <- simCoverage(nsim = nSim, n = 50, scenario = 3, method="SPLRT", ncores=64)

cov.1.SPLRT.100 <- simCoverage(nsim = nSim, n = 100, scenario = 1, method="SPLRT", ncores=64)
cov.2.SPLRT.100 <- simCoverage(nsim = nSim, n = 100, scenario = 2, method="SPLRT", ncores=64)
cov.3.SPLRT.100 <- simCoverage(nsim = nSim, n = 100, scenario = 3, method="SPLRT", ncores=64)

cov.1.SPLRT.200 <- simCoverage(nsim = nSim, n = 200, scenario = 1, method="SPLRT", ncores=64)
cov.2.SPLRT.200 <- simCoverage(nsim = nSim, n = 200, scenario = 2, method="SPLRT", ncores=64)
cov.3.SPLRT.200 <- simCoverage(nsim = nSim, n = 200, scenario = 3, method="SPLRT", ncores=64)



cov.1.SPLRT.50
cov.2.SPLRT.50
cov.3.SPLRT.50

cov.1.SPLRT.100
cov.2.SPLRT.100
cov.3.SPLRT.100

cov.1.SPLRT.200
cov.2.SPLRT.200
cov.3.SPLRT.200

