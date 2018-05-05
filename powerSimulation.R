rm(list=ls())

library(TruncComp)
library(parallel)

library(devtools)
install_github('aejensen/TruncComp')

powerSim <- function(n, scenario, nSim = 500, ncores=6) {
  out <- parallel::mclapply(1:nSim, function(sim) {
    if(scenario == 1) {
      data <- TruncComp:::simTruncData(n, 3, 4, 0.35, 0.35, dist="norm")
    } else if(scenario == 2) {
      data <- TruncComp:::simTruncData(n, 3.5, 3.5, 0.4, 0.3, dist="norm")
    } else if(scenario == 3) {
      data <- TruncComp:::simTruncData(n, 3, 4, 0.4, 0.3, dist="norm")
    } else if(scenario == 4) {
      data <- TruncComp:::simTruncData(n, 3, 4, 0.4, 0.3, dist="t-sq", df=2)
    }
    m.SPLRT <- truncComp.default(data$Y, data$A, data$R, method="SPLRT")
    m.LRT <- truncComp.default(data$Y, data$A, data$R, method="LRT")
    m.wil <- stats::wilcox.test(data$Y ~ data$R)

    c(m.wil$p.value < 0.05, m.LRT$p < 0.05, m.SPLRT$p < 0.05)
  }, mc.cores=ncores)
  out <- do.call("rbind", out)
  apply(out, 2, mean, na.rm=TRUE)
}

nSim <- 500
nSeq <- seq(40, 250, length.out=7)

power4 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 4, nSim=nSim, ncores=64))



power1 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 1, nSim=nSim, ncores=6))
power2 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 2, nSim=nSim, ncores=6))
power3 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 3, nSim=nSim, ncores=6))
power4 <- sapply(nSeq, function(n) powerSim(n = n, scenario = 4, nSim=nSim, ncores=6))
out


plot(nSeq, power1[1,], ylim=c(0, 1), type="l", col="green", bty="n")
lines(nSeq, power1[2,], col="black")
lines(nSeq, power1[3,], col="red")
