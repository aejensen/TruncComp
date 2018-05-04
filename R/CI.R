jointContrastLRT <- function(data, muDelta, alphaDelta) {
  yAlive1 <- data[data$R == 0 & data$A == 1, "Y"]
  yAlive2 <- data[data$R == 1 & data$A == 1, "Y"]
  ELRT <- EL::EL.means(yAlive2, yAlive1, mu = muDelta)

  binom <- TruncComp:::logit.LRT(data, alphaDelta)

  as.numeric(ELRT$statistic + binom)
}

jointContrastCI <- function(model, muDelta = seq(-1, 1, length.out = 20), piDelta = seq(-3, 3, length.out = 20), plot=TRUE, conf.level = 0.95) {
  data <- model$data

  matOut <- matrix(NA, length(muDelta), length(piDelta))
  for(a in 1:length(muDelta)) {
    for(b in 1:length(piDelta)) {
      matOut[a,b] <- jointContrastLRT(data, muDelta[a], piDelta[b])
    }
  }

  if(plot) {
    fields::image.plot(muDelta, piDelta, matOut, xlab="Difference in mean", ylab="log OR")
    points(model$muDelta, log(model$alphaDelta), pch=19, cex=3)
    points(0, 0, cex=3, pch=3)
    contour(muDelta, piDelta, matOut, add=TRUE, levels=stats::qchisq(conf.level, 2), lwd=2)
  }
  matOut
}



