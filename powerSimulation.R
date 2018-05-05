rm(list=ls())


library(devtools)
install_github('aejensen/TruncComp')

library(TruncComp)
library(parallel)

powerSim <- function(n, scenario, nSim = 500, alpha = 0.05, ncores=64, seed=12345) {
  set.seed(seed)

  out <- parallel::mclapply(1:nSim, function(sim) {
    #Define the four scenarios
    if(scenario == 1) {
      f0 <- function(n) stats::rnorm(n, 3, 1)
      f1 <- function(n) stats::rnorm(n, 3.5, 1)
      pi0 <- 0.35
      pi1 <- 0.35
    } else if(scenario == 2) {
      f0 <- function(n) stats::rnorm(n, 3.5, 1)
      f1 <- function(n) stats::rnorm(n, 3.5, 1)
      pi0 <- 0.4
      pi1 <- 0.3
    } else if(scenario == 3) {
      f0 <- function(n) stats::rnorm(n, 3, 1)
      f1 <- function(n) stats::rnorm(n, 3.5, 1)
      pi0 <- 0.4
      pi1 <- 0.3
    } else if(scenario == 4) {
      #f0 <- function(n) stats::rnorm(n, 8, 1)
      #f1 <- function(n) stats::rexp(n, 0.1)
      f0 <- function(n) stats::rnorm(n, 5, 2)
      f1 <- function(n) rbeta(n, 0.5, 0.5) + 5
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
    c(wilcox = m.wil$p.value < alpha,
      LRT = m.LRT$p < alpha,
      SPLRT = m.SPLRT$p < alpha)
  }, mc.cores=ncores, mc.preschedule = TRUE)
  out <- do.call("rbind", out)
  apply(out, 2, mean, na.rm=TRUE)
}

nSim <- 10^4                        #number of simulations
nSeq <- seq(50, 350, length.out=16) #different sample sizes

power1 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 1, nSim=nSim, ncores=64))
power2 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 2, nSim=nSim, ncores=64))
power3 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 3, nSim=nSim, ncores=64))
power4 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 4, nSim=nSim, ncores=64))

save(nSeq, power1, power2, power3, power4, file="powerSimResults.RData")



#####
par(mfrow=c(2,2), mgp=c(2,1,0), mar=c(3.5, 3, 1, 0.5))

plot(nSeq, power1[1,], ylim=c(0, 1), type="l", col="cornflowerblue", bty="n",
     ylab="Power", xlab="Sample size per group", main="Setup #1", lwd=2)
lines(nSeq, power1[2,], col="firebrick1", lwd=2)
lines(nSeq, power1[3,], col="black", lwd=2)
legend("topleft", c("Parametric", "Semi-parametric", "Wilcoxon"),
       lwd=2, lty=1, col=c("firebrick1", "black", "cornflowerblue"), bty="n")

plot(nSeq, power2[1,], ylim=c(0, 1), type="l", col="cornflowerblue", bty="n",
     ylab="Power", xlab="Sample size per group", main="Setup #2", lwd=2)
lines(nSeq, power2[2,], col="firebrick1", lwd=2)
lines(nSeq, power2[3,], col="black", lwd=2)

plot(nSeq, power3[1,], ylim=c(0, 1), type="l", col="cornflowerblue", bty="n",
    ylab="Power", xlab="Sample size per group", main="Setup #3", lwd=2)
lines(nSeq, power3[2,], col="firebrick1", lwd=2)
lines(nSeq, power3[3,], col="black", lwd=2)

plot(nSeq, power4[1,], ylim=c(0, 1), type="l", col="cornflowerblue", bty="n",
     ylab="Power", xlab="Sample size per group", main="Setup #4", lwd=2)
lines(nSeq, power4[2,], col="firebrick1", lwd=2)
lines(nSeq, power4[3,], col="black", lwd=2)
