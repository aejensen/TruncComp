rm(list=ls())

library(TruncComp)
library(parallel)

library(devtools)
install_github('aejensen/TruncComp')

powerSim <- function(n, scenario, nSim = 500, alpha = 0.05, ncores=64, seed=12345) {
  set.seed(12345)

  out <- parallel::mclapply(1:nSim, function(sim) {
    #Define the four scenarios
    if(scenario == 1) {
      f0 <- function(n) stats::rnorm(n, 3, 1)
      f1 <- function(n) stats::rnorm(n, 4, 1)
      pi0 <- 0.35
      pi1 <- 0.35
    } else if(scenario == 2) {
      f0 <- function(n) stats::rnorm(n, 3.5, 1)
      f1 <- function(n) stats::rnorm(n, 3.5, 1)
      pi0 <- 0.4
      pi1 <- 0.3
    } else if(scenario == 3) {
      f0 <- function(n) stats::rnorm(n, 3, 1)
      f1 <- function(n) stats::rnorm(n, 4, 1)
      pi0 <- 0.4
      pi1 <- 0.3
    } else if(scenario == 4) {
      f0 <- function() stats::rnorm(n, 8, 1)
      f1 <- function() stats::rexp(n, 0.1)
      pi0 <- 0.4
      pi1 <- 0.3
    }

    #Simulate data
    data <- TruncComp::simulateTruncatedData(n, f0, f1, pi0, pi1)

    #Fit models
    m.SPLRT <- TruncComp::truncComp.default(data$Y, data$A, data$R, method="SPLRT")
    m.LRT <- TruncComp::truncComp.default(data$Y, data$A, data$R, method="LRT")
    m.wil <- stats::wilcox.test(data$Y ~ data$R)

    #Get H0 rejection or not
    c(m.wil$p.value < alpha, m.LRT$p < alpha, m.SPLRT$p < alpha)
  }, mc.cores=ncores)
  out <- do.call("rbind", out)
  apply(out, 2, mean, na.rm=TRUE)
}

nSim <- 1000 #Number of simulations
nSeq <- seq(40, 250, length.out=7)

power1 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 1, nSim=nSim, ncores=64))
power2 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 2, nSim=nSim, ncores=64))
power3 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 3, nSim=nSim, ncores=64))
power4 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 4, nSim=nSim, ncores=64))



plot(nSeq, power1[1,], ylim=c(0, 1), type="l", col="green", bty="n")
lines(nSeq, power1[2,], col="black")
lines(nSeq, power1[3,], col="red")
