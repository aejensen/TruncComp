#curve(gamlss.dist::dSN2(x, mu = -4, sigma=1, nu=3), -12, 12)
#curve(gamlss.dist::dSN2(x, mu = 3, sigma=1, nu=0.2), -12, 12, add=TRUE, col=2)
#f0 <- function() gamlss.dist::rSN2(1, mu = -4, sigma=1, nu=3)
#f1 <- function() gamlss.dist::rSN2(1, mu = 3, sigma=1, nu=0.2)
#mean(sapply(1:10^4, function(q) f1())) - mean(sapply(1:10^4, function(q) f0()))

#nu <- 2
#tau <- 1

#mD <- uniroot(function(m) {
#  i1 <- integrate(function(x) x * gamlss.dist::dBCT(x, mu = m, sigma=2, nu=nu, tau=tau), 0, Inf)$value
#  i2 <- integrate(function(x) x * gamlss.dist::dBCT(x, mu = 3, sigma=2, nu=nu, tau=tau), 0, Inf)$value
#  i1 - i2 - 1
#}, interval = c(2.5, 4))$root


#curve(gamlss.dist::dBCT(x, mu = 3, sigma=2, nu=nu, tau=tau), 0, 24)
#curve(gamlss.dist::dBCT(x, mu = mD, sigma=2, nu=nu, tau=tau), 0, 24, add=TRUE, col=2)

#f0 <- function() gamlss.dist::rBCT(1, mu = 3, sigma=2, nu=nu, tau=tau)
#f1 <- function() gamlss.dist::rBCT(1, mu = mD, sigma=2, nu=nu, tau=tau)

#integrate(function(x) x * gamlss.dist::dBCT(x, mu = mD, sigma=2, nu=nu, tau=tau), 0, Inf)$value - integrate(function(x) x * gamlss.dist::dBCT(x, mu = 3, sigma=2, nu=nu, tau=tau), 0, Inf)$value

#hist(sapply(1:10^4, function(q) f0()))
#hist(sapply(1:10^4, function(q) f1()))

f0 <- function() rnorm(n, 8, 1)
f1 <- function() rexp(n, 0.1)

mean(sapply(1:10^4, function(q) f1())) - mean(sapply(1:10^4, function(q) f0()))
hist(sapply(1:10^4, function(q) f0()))
hist(sapply(1:10^4, function(q) f1()))

plot(density(sapply(1:10^4, function(q) f0())))
lines(density(sapply(1:10^4, function(q) f1())), col=2)

test <- mclapply(1:10000, function(nsim) {
  data <- TruncComp:::simTruncData2(25, f0, f1, 0.4, 0.3)
  m <- truncComp.default(data$Y, data$A, data$R, method="SPLRT")
  m$muDelta
}, mc.cores=64, mc.preschedule=FALSE)
mean(unlist(test), na.rm=TRUE)

mean(unlist(test[sapply(test, is.numeric)]), na.rm=TRUE)


test <- mclapply(1:5000, function(hej) {
  mean(sapply(1:25, function(q) f1())) - mean(sapply(1:25, function(q) f0()))
}, mc.cores=64)
mean(unlist(test))
mTrue <-11137
mean(unlist(test)) - mTrue
((mean(unlist(test)) - mTrue) / mTrue) * 100



powerSim <- function(n, nSim = 500, ncores=6) {
  out <- parallel::mclapply(1:nSim, function(sim) {
    data <- TruncComp:::simTruncData2(n, f0, f1, 0.4, 0.3)

    m.SPLRT <- truncComp.default(data$Y, data$A, data$R, method="SPLRT")
    m.LRT <- truncComp.default(data$Y, data$A, data$R, method="LRT")
    m.wil <- stats::wilcox.test(data$Y ~ data$R)

    c(wilcox = m.wil$p.value < 0.05, LRT = m.LRT$p < 0.05, SPLRT = m.SPLRT$p < 0.05)
  }, mc.cores=ncores)
  out <- do.call("rbind", out)
  apply(out, 2, mean, na.rm=TRUE)
}

nSeq <- seq(40, 250, length.out=7)
power4 <- sapply(nSeq, function(n) powerSim(n = n, nSim=2000, ncores=64))
power4

#a <- c(0.129, 0.212, 0.317, 0.399, 0.449, 0.551, 0.572)
#b <- c(0.203, 0.284, 0.366, 0.444, 0.495, 0.590, 0.604)
#plot(nSeq, a, ylim=c(0,1), type="l")
#lines(nSeq, b)



simTest <- function(n, nsim=5000) {
  out <- mclapply(1:nsim, function(sim) {
    #y1Mix <- sample(0:1, n, replace=TRUE, prob=c(0.5, 0.5))
    #y1 <- rnorm(n, 3, 1)*y1Mix + rnorm(n, 8, 1)*(1-y1Mix)
    #y1 <- rnorm(n, 4, 1)
    y1 <- rnorm(n, 8, 1)

    #y2Mix <- sample(0:1, n, replace=TRUE, prob=c(0.5, 0.5))
    #y2 <- rnorm(n, 6, 1)*y2Mix + rexp(n, 1)*(1-y2Mix)
    y2 <- rexp(n, 0.1)
    mean(y2) - mean(y1)

    ELRT <- EL::EL.means(y1, y2)
    TT <- t.test(y1, y2)
    c(EL = ELRT$p.value < 0.05, TT = TT$p.value < 0.05)
  }, mc.cores=64)
  out <- do.call("rbind", out)
  colMeans(out)
}

simTest(25)
simTest(50)
simTest(100)
simTest(250)
hist(y2)



#3.25, 1
nSeq <- seq(25, 250, length.out=10)
out <- sapply(nSeq, function(n) simTest(n))

