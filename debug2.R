rm(list=ls())

n <- 10^4

d <- TruncComp:::simTruncData(n, 3, 4, 0.6, 0.8)

m1 <- mean(d[d$A==1 & d$R == 1,]$Y) #treatment mean among alive
m0 <- mean(d[d$A==1 & d$R == 0,]$Y) #control mean among live


p1.mid.1 <- mean(d[d$R == 1,]$A == 1)
p0.mid.1 <- mean(d[d$R == 0,]$A == 1)
p1.mid.0 <- mean(d[d$R == 1,]$A == 0)
p0.mid.0 <- mean(d[d$R == 0,]$A == 0)

q1 <- m1 * p1.mid.1 + 0 * p1.mid.0
q2 <- m0 * p0.mid.1 + 0 * p0.mid.0

q1 - q2

mean(d[d$R == 1,]$Y) - mean(d[d$R == 0,]$Y)




(m1-m0 ) * (p0.mid.1 - p0.mid.1)

mean(d[d$R == 1,]$Y) - mean(d[d$R == 0,]$Y)
m1 * p1.mid.1 - m0 * p0.mid.1



d <- TruncComp:::simTruncData(100, 3, 4, 0.6, 0.8)
sub1 <- subset(d, d$R == 1)


mean(sub1$Y)
mean(subset(sub1, A == 1)$Y)*mean(sub1$A == 1) + mean(subset(sub1, A == 0)$Y)*mean(sub1$A == 0)
mean(subset(sub1, A == 1)$Y)*mean(sub1$A == 1) + mean(subset(sub1, A == 0)$Y)*(1 - mean(sub1$A == 1))

mean(sub1$A == 0)





mean(subset(sub1, A == 1)$Y) * (mean(d$A == 1))
