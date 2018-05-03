jointContrastCI <- function(data, muDelta, alphaDelta) {
  yAlive1 <- data[data$R == 0 & data$A == 1, "Y"]
  yAlive2 <- data[data$R == 1 & data$A == 1, "Y"]
  ELRT <- EL::EL.means(yAlive2, yAlive1, mu = muDelta)

  binom <- TruncComp:::logit.LRT(data, alphaDelta)

  as.numeric(ELRT$statistic + binom)
}



