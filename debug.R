#########################
rm(list=ls())

par(mfrow=c(3,3))
for(k in 1:9) {
  data <- TruncComp:::simTruncData(50, 0, 0.3, 0.85, 0.9)
  m <- truncComp.default(data$Y, data$A, data$R, method="SPLRT")
  a <- jointContrastCI(m)
}

simCoverage <- function(nsim, n, method, seed=12345) {
  set.seed(seed)
  coverageOut <- t(sapply(1:nsim, function(sim) {
    #cat(sim, "\n")
    data <- TruncComp:::simTruncData(n, 0, 0.3, 0.85, 0.9)
    m <- truncComp.default(data$Y, data$A, data$R, method=method)

    muDeltaTrue <- 0.3
    orDeltaTrue <- (0.9/(1-0.9)) / (0.85/(1-0.85))

    coverage.muDelta <- (m$muDeltaCI[1] < muDeltaTrue) & m$muDeltaCI[2] > muDeltaTrue
    coverage.OR <- (m$alphaDeltaCI[1] < orDeltaTrue) & (m$alphaDeltaCI[2] > orDeltaTrue)

    c(coverage.muDelta, coverage.OR)
  }))
  c(mean(coverageOut[,1] == TRUE, na.rm=TRUE), mean(coverageOut[,2] == TRUE, na.rm=TRUE))
}

simCoverage(nsim = 500, n = 30, method="SPLRT")
simCoverage(nsim = 500, n = 50, method="SPLRT")
simCoverage(nsim = 500, n = 100, method="SPLRT")



data <- TruncComp:::simTruncData(100000, 0, 0.3, 0.6, 0.9)
m <- truncComp.default(data$Y, data$A, data$R, method="SPLRT")
m$alphaDelta


image.plot(muDeltaSeq, piDeltaSeq, matOut)
points(m$muDelta, log(m$alphaDelta), pch=19, cex=4)
points(0, 0, cex=5, pch=19)

contour(muDeltaSeq, piDeltaSeq, matOut, add=FALSE, levels=stats::qchisq(1-0.05, 2))
contour(muDeltaSeq, piDeltaSeq, matOut, add=TRUE, levels=stats::qchisq(1-0.2, 2))




lines(c(m$muDelta, m$muDelta), c(0.331035, -0.07047824), lwd=2)
0.7249986

#######

deltaSeq <- seq(-1, 3, length.out=50)
plot(deltaSeq, sapply(deltaSeq, function(delta) TruncComp:::logit.LRT(d, delta)), type="l")
abline(h = stats::qchisq(1-0.05, 1))
abline(v=log(m$alphaDeltaCI[1]))
abline(v=log(m$alphaDeltaCI[2]))

TruncComp:::logit.LRT(d, 0)
anova(stats::glm(A ~ 1, family=stats::binomial(), data=d),
      stats::glm(A ~ Z, family=stats::binomial(), data=d), test="LRT")$Deviance[2]



#######################
d <- simTruncData(50, 1, 2, 0.4, 0.6)
m <- truncComp(d$Y, d$A, d$Z, method="SPLRT")

m$ELRT
m$muW
m$ELRT$statistic

yAlive1 <- d[d$Z == 0 & d$A == 1, "Y"]
yAlive2 <- d[d$Z == 1 & d$A == 1, "Y"]

m$yAlive1 == yAlive1
m$yAlive2 == yAlive2

el <- EL::EL.means(yAlive2, yAlive1, mu=0)
el$statistic
m$ELRT$statistic

EL::EL.means(yAlive2, yAlive1, mu = 0, conf.level = 0.95)$statistic
m$ELRT$statistic


?EL::EL.means

head(d)


#######################

simOut <- sapply(1:1000, function(q) {
  d <- simTruncData(50, 1, 2, 0.4, 0.6)
  test <- truncComp(d$Y, d$A, d$Z, method="SPLRT")
  c(test$muDelta, test$alphaDelta)
})

plot(t(simOut))

d <- simTruncData(50, 1, 2, 0.4, 0.6)
alive <- d[d$A == 1,]

EL::EL.means(alive$Y[alive$Z == 1], alive$Y[alive$Z == 0])
truncComp(d$Y, d$A, d$Z, method="SPLRT")



deltaSeq <- seq(0.5, 2, length.out=100)
out <- sapply(deltaSeq, function(d) EL::EL.means(alive$Y[alive$Z == 1], alive$Y[alive$Z == 0], mu=d)$statistic)
plot(deltaSeq, out)
abline(v=1.115)
abline(h=3.84)
abline(v=0.7244)
abline(v=1.55)

stats::qchisq(1-0.05, 2)



hejhej <- function(data, muDelta, alphaDelta) {
  yAlive1 <- data[data$R == 0 & data$A == 1, "Y"]
  yAlive2 <- data[data$R == 1 & data$A == 1, "Y"]
  ELRT <- EL::EL.means(yAlive2, yAlive1, mu = muDelta)

  binom <- TruncComp:::logit.LRT(data, alphaDelta)

  as.numeric(ELRT$statistic + binom)
}
