rm(list=ls())

library(devtools)
install_github('aejensen/TruncComp')

library(TruncComp)
library(parallel)
library(sp)

############################################################
# POINT ESTIMATION SIMULATION FOR THE MANUSCRIPT
############################################################
pointEstimateSim <- function(nsim, n, scenario, method, ncores, seed=12345) {
  set.seed(seed)

  out <- parallel::mclapply(1:nsim, function(sim) {
    #cat(sim, " ")

    if(scenario == 1) {
      f0 <- function(n) stats::rnorm(n, 3, 1)
      f1 <- function(n) stats::rnorm(n, 3.5, 1)
      pi0 <- 0.35
      pi1 <- 0.35
      muDeltaTrue <- 0.5
    } else if(scenario == 2) {
      f0 <- function(n) stats::rnorm(n, 3.5, 1)
      f1 <- function(n) stats::rnorm(n, 3.5, 1)
      pi0 <- 0.4
      pi1 <- 0.3
      muDeltaTrue <- 0
    } else if(scenario == 3) {
      f0 <- function(n) stats::rnorm(n, 3, 1)
      f1 <- function(n) stats::rnorm(n, 3.5, 1)
      pi0 <- 0.4
      pi1 <- 0.3
      muDeltaTrue <- 0.5
    } else if(scenario == 4) {
      f0 <- function(n) stats::rnorm(n, 5, 2)
      f1 <- function(n) stats::rbeta(n, 0.5, 0.5) + 5
      pi0 <- 0.4
      pi1 <- 0.3
      muDeltaTrue <- 0.5
    }
    orDeltaTrue <- (pi1/(1-pi1)) / (pi0/(1-pi0))

    #Simulate data
    data <- TruncComp::simulateTruncatedData(n, f0, f1, pi0, pi1)

    #Fit model
    m <- TruncComp::truncComp(Y ~ R, atom = 0, data = data, method = method)

    muDelta.width <-  abs(m$muDeltaCI[2] - m$muDeltaCI[1])
    logOR.width <- abs(m$alphaDeltaCI[2] - m$alphaDeltaCI[1])

    #Get marginal inclusion indicators
    inclusion.muDelta <- (m$muDeltaCI[1] < muDeltaTrue) & m$muDeltaCI[2] > muDeltaTrue
    inclusion.OR <- (m$alphaDeltaCI[1] < orDeltaTrue) & (m$alphaDeltaCI[2] > orDeltaTrue)

    if(method == "SPLRT") {
      #jSurf <- suppressMessages(confint(m, type="simultaneous", plot=FALSE, resolution = 40)) #40
      #cont <- contourLines(jSurf$muDelta, jSurf$logORdelta, jSurf$surface, levels=stats::qchisq(0.95, 2))
      #inclusion.joint <- sp::point.in.polygon(muDeltaTrue, log(orDeltaTrue), cont[[1]]$x, cont[[1]]$y) == 1
      inclusion.joint <- NA
    } else {
      inclusion.joint <- NA #Not implemented for parametric likelihood
    }

    c(m$muDelta, log(m$alphaDelta), inclusion.muDelta, inclusion.OR,
      inclusion.muDelta, inclusion.OR, inclusion.joint)
  }, mc.cores=ncores, mc.preschedule = TRUE)

  out
}

nSim <- 25 * 10^3

#Setup 1
LRT.1.50 <- pointEstimateSim(nsim = nSim, n = 50, scenario = 1, method="LRT", ncores=64)
SPLRT.1.50 <- pointEstimateSim(nsim = nSim, n = 50, scenario = 1, method="SPLRT", ncores=64)

LRT.1.100 <- pointEstimateSim(nsim = nSim, n = 100, scenario = 1, method="LRT", ncores=64)
SPLRT.1.100 <- pointEstimateSim(nsim = nSim, n = 100, scenario = 1, method="SPLRT", ncores=64)

LRT.1.200 <- pointEstimateSim(nsim = nSim, n = 200, scenario = 1, method="LRT", ncores=64)
SPLRT.1.200 <- pointEstimateSim(nsim = nSim, n = 200, scenario = 1, method="SPLRT", ncores=64)

#Setup 2
LRT.2.50 <- pointEstimateSim(nsim = nSim, n = 50, scenario = 2, method="LRT", ncores=64)
SPLRT.2.50 <- pointEstimateSim(nsim = nSim, n = 50, scenario = 2, method="SPLRT", ncores=64)

LRT.2.100 <- pointEstimateSim(nsim = nSim, n = 100, scenario = 2, method="LRT", ncores=64)
SPLRT.2.100 <- pointEstimateSim(nsim = nSim, n = 100, scenario = 2, method="SPLRT", ncores=64)

LRT.2.200 <- pointEstimateSim(nsim = nSim, n = 200, scenario = 2, method="LRT", ncores=64)
SPLRT.2.200 <- pointEstimateSim(nsim = nSim, n = 200, scenario = 2, method="SPLRT", ncores=64)

#Setup 3
LRT.3.50 <- pointEstimateSim(nsim = nSim, n = 50, scenario = 3, method="LRT", ncores=64)
SPLRT.3.50 <- pointEstimateSim(nsim = nSim, n = 50, scenario = 3, method="SPLRT", ncores=64)

LRT.3.100 <- pointEstimateSim(nsim = nSim, n = 100, scenario = 3, method="LRT", ncores=64)
SPLRT.3.100 <- pointEstimateSim(nsim = nSim, n = 100, scenario = 3, method="SPLRT", ncores=64)

LRT.3.200 <- pointEstimateSim(nsim = nSim, n = 200, scenario = 3, method="LRT", ncores=64)
SPLRT.3.200 <- pointEstimateSim(nsim = nSim, n = 200, scenario = 3, method="SPLRT", ncores=64)

#Setup 4
LRT.4.50 <- pointEstimateSim(nsim = nSim, n = 50, scenario = 4, method="LRT", ncores=64)
SPLRT.4.50 <- pointEstimateSim(nsim = nSim, n = 50, scenario = 4, method="SPLRT", ncores=64)

LRT.4.100 <- pointEstimateSim(nsim = nSim, n = 100, scenario = 4, method="LRT", ncores=64)
SPLRT.4.100 <- pointEstimateSim(nsim = nSim, n = 100, scenario = 4, method="SPLRT", ncores=64)

LRT.4.200 <- pointEstimateSim(nsim = nSim, n = 200, scenario = 4, method="LRT", ncores=64)
SPLRT.4.200 <- pointEstimateSim(nsim = nSim, n = 200, scenario = 4, method="SPLRT", ncores=64)

save(LRT.1.50, SPLRT.1.50, LRT.1.100, SPLRT.1.100, LRT.1.200, SPLRT.1.200,
     file = "coverageSimTest.RData")


save(LRT.1.50, SPLRT.1.50, LRT.1.100, SPLRT.1.100, LRT.1.200, SPLRT.1.200,
     LRT.2.50, SPLRT.2.50, LRT.2.100, SPLRT.2.100, LRT.2.200, SPLRT.2.200,
     LRT.3.50, SPLRT.3.50, LRT.3.100, SPLRT.3.100, LRT.3.200, SPLRT.3.200,
     LRT.4.50, SPLRT.4.50, LRT.4.100, SPLRT.4.100, LRT.4.200, SPLRT.4.200,
     file = "coverageSimulationNoJoint.RData")


############################################################
# GATHER THE RESULTS
############################################################
load("coverageSimulationNoJoint.RData")

getStats <- function(res, muDelta, logOR) {
  mat <- do.call("rbind", res)
  out <- c(biasMuDelta = mean(mat[,1] - muDelta, na.rm=TRUE),
           sdMuDelta = sd(mat[,1], na.rm=TRUE),
           biasLogOR = mean(mat[,2] - logOR),
           sdLogOR = sd(mat[,1]),
           widthMuDelta = mean(mat[,3]),
           widthLogOR = mean(mat[,4]),
           covMuDelta = mean(mat[,5])*100,
           covLogOR = mean(mat[,6])*100,
           covJoint = mean(mat[,7])*100)
  round(out, 4)
}

rbind(LRT = getStats(LRT.1.50, 0.5, 0), SPLRT = getStats(SPLRT.1.50, 0.5, 0))
rbind(LRT = getStats(LRT.1.100, 0.5, 0), SPLRT = getStats(SPLRT.1.100, 0.5, 0))
rbind(LRT = getStats(LRT.1.200, 0.5, 0), SPLRT = getStats(SPLRT.1.200, 0.5, 0))

rbind(LRT = getStats(LRT.2.50, 0, -0.4418328), SPLRT = getStats(SPLRT.2.50, 0, -0.4418328))
rbind(LRT = getStats(LRT.2.100, 0, -0.4418328), SPLRT = getStats(SPLRT.2.100, 0, -0.4418328))
rbind(LRT = getStats(LRT.2.200, 0, -0.4418328), SPLRT = getStats(SPLRT.2.200, 0, -0.4418328))

rbind(LRT = getStats(LRT.3.50, 0.5, -0.4418328), SPLRT = getStats(SPLRT.3.50, 0.5, -0.4418328))
rbind(LRT = getStats(LRT.3.100, 0.5, -0.4418328), SPLRT = getStats(SPLRT.3.100, 0.5, -0.4418328))
rbind(LRT = getStats(LRT.3.200, 0.5, -0.4418328), SPLRT = getStats(SPLRT.3.200, 0.5, -0.4418328))

rbind(LRT = getStats(LRT.4.50, 0.5, -0.4418328), SPLRT = getStats(SPLRT.4.50, 0.5, -0.4418328))
rbind(LRT = getStats(LRT.4.100, 0.5, -0.4418328), SPLRT = getStats(SPLRT.4.100, 0.5, -0.4418328))
rbind(LRT = getStats(LRT.4.200, 0.5, -0.4418328), SPLRT = getStats(SPLRT.4.200, 0.5, -0.4418328))


