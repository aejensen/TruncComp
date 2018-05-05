rm(list=ls())

library(devtools)
install_github('aejensen/TruncComp')

library(TruncComp)
library(parallel)

############################################################
# POWER SIMULATIONS FOR THE MANUSCRIPT
############################################################
powerSim <- function(n, scenario, nSim, alpha = 0.05, ncores=64, seed=12345) {
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
      f1 <- function(n) stats::rbeta(n, 0.5, 0.5) + 5
      pi0 <- 0.4
      pi1 <- 0.3
    }

    #Simulate data
    data <- TruncComp::simulateTruncatedData(n, f0, f1, pi0, pi1)

    #Fit models
    m.SPLRT <- TruncComp::truncComp(Y ~ R, atom = 0, data = data, method = "SPLRT")
    m.LRT <- TruncComp::truncComp(Y ~ R, atom = 0, data = data, method = "LRT")
    m.wil <- stats::wilcox.test(Y ~ R, data = data)

    #Get H0 rejection or not
    c(wilcox = m.wil$p.value < alpha,
      LRT = m.LRT$p < alpha,
      SPLRT = m.SPLRT$p < alpha)
  }, mc.cores=ncores, mc.preschedule = TRUE)

  out <- do.call("rbind", out)
  apply(out, 2, mean, na.rm=TRUE)
}

nSim <- 25 * 10^3                   #number of simulations
nSeq <- seq(50, 350, length.out=16) #different sample sizes

power1 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 1, nSim=nSim, ncores=64))
power2 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 2, nSim=nSim, ncores=64))
power3 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 3, nSim=nSim, ncores=64))
power4 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 4, nSim=nSim, ncores=64))

save(nSeq, power1, power2, power3, power4, file="powerSimResults.RData")



############################################################
# PLOT THE RESULTS
############################################################
load("powerSimResults.RData")

par(mfrow=c(2,2), mgp=c(2,1,0), mar=c(3.5, 3, 1, 0.5))
plot(nSeq, power1[1,], ylim=c(0, 1), type="l", col="firebrick1", bty="n",
     ylab="Power", xlab="Sample size per group", main="Setup 1", lwd=2)
lines(nSeq, power1[2,], col="cornflowerblue", lwd=2)
lines(nSeq, power1[3,], col="black", lwd=2)
legend("topleft", c("Parametric", "Semi-parametric", "Wilcoxon"),
       lwd=2, lty=1, col=c("cornflowerblue", "black", "firebrick1"), bty="n")

plot(nSeq, power2[1,], ylim=c(0, 1), type="l", col="firebrick1", bty="n",
     ylab="Power", xlab="Sample size per group", main="Setup 2", lwd=2)
lines(nSeq, power2[2,], col="cornflowerblue", lwd=2)
lines(nSeq, power2[3,], col="black", lwd=2)

plot(nSeq, power3[1,], ylim=c(0, 1), type="l", col="firebrick1", bty="n",
    ylab="Power", xlab="Sample size per group", main="Setup 3", lwd=2)
lines(nSeq, power3[2,], col="cornflowerblue", lwd=2)
lines(nSeq, power3[3,], col="black", lwd=2)

plot(nSeq, power4[1,], ylim=c(0, 1), type="l", col="firebrick1", bty="n",
     ylab="Power", xlab="Sample size per group", main="Setup 4", lwd=2)
lines(nSeq, power4[2,], col="cornflowerblue", lwd=2)
lines(nSeq, power4[3,], col="black", lwd=2)


par(mgp=c(2.1,1,0))
plot(nSeq, (power1[3,] - power1[2,])/power1[2,] * 100, type="l", bty="n", ylim=c(0,30),
     col="black", lwd=2, xlab="Sample size per group",
     ylab="Percentage relative difference in power")
lines(nSeq, (power2[3,] - power2[2,])/power2[2,] * 100, col="cornflowerblue", lwd=2)
lines(nSeq, (power3[3,] - power3[2,])/power3[2,] * 100, col="darkgoldenrod1", lwd=2)
lines(nSeq, (power4[3,] - power4[2,])/power4[2,] * 100, col="firebrick1", lwd=2)
legend("topleft", c("Setup 1", "Setup 2", "Setup 3", "Setup 4"),
       lwd=2, col=c("black", "cornflowerblue", "darkgoldenrod1", "firebrick1"),
       bty="n")

