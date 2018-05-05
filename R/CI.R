jointContrastLRT <- function(data, muDelta, alphaDelta) {
  yAlive1 <- data[data$R == 0 & data$A == 1, "Y"]
  yAlive2 <- data[data$R == 1 & data$A == 1, "Y"]
  ELRT <- EL::EL.means(yAlive2, yAlive1, mu = muDelta)

  binom <- TruncComp:::logit.LRT(data, alphaDelta)

  as.numeric(ELRT$statistic + binom)
}

jointContrastCI <- function(m, muDelta = NULL, logORdelta = NULL,
                            plot=TRUE, conf.level = 0.95) {

  if(!("TruncComp" %in% class(m))) {
    stop("m must be an object of type TruncComp")
  }

  if(is.null(muDelta)) {
    muDelta <- seq(m$muDeltaCI[1] -0.4 , m$muDeltaCI[2] + 0.4, length.out = 40)
  }

  if(is.null(logORdelta)) {
    logORdelta = seq(log(m$alphaDeltaCI[1]) - 0.4, log(m$alphaDeltaCI[2]) + 0.4, length.out = 40)
  }

  matOut <- matrix(NA, length(muDelta), length(logORdelta))
  for(a in 1:length(muDelta)) {
    for(b in 1:length(logORdelta)) {
      matOut[a,b] <- jointContrastLRT(m$data, muDelta[a], logORdelta[b])
    }
  }

  if(plot) {
    fields::image.plot(muDelta, logORdelta, matOut, xlab="Difference in mean", ylab="log OR")
    points(m$muDelta, log(m$alphaDelta), pch=19, cex=3)
    points(0, 0, cex=3, pch=3)
    contour(muDelta, logORdelta, matOut, add=TRUE, levels=stats::qchisq(conf.level, 2), lwd=2)
  }

  list(muDelta = muDelta, logORdelta = logORdelta, surface = matOut)
}



