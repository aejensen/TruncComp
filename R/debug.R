#########################
rm(list=ls())

d <- simTruncData(50, 1, 2, 0.4, 0.6)
m <- truncComp(d$Y, d$A, d$Z, method="SPLRT")
m


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

