mu0 <- 5
pi0 <- 0.7
pi1 <- 0.5
muDelta <- mu0*((pi0 / pi1) - 1)

f0 <- function(n) stats::rnorm(n, mu0, 1)
f1 <- function(n) stats::rnorm(n, mu0 + muDelta, 1)

d <- TruncComp::simulateTruncatedData(25, f0, f1, pi0, pi1)
d$A <- NULL

save(TruncCompExample, file="TruncCompExample.RData")
