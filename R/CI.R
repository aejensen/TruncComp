confint.TruncComp <- function(object, type = "marginal", muDelta = NULL, logORdelta = NULL, conf.level = 0.95, plot=TRUE, offset = 1, resolution = 35,  ...) {
  if(!(type == "marginal" | type == "simultaneous")) {
    stop("Type of confidence interval must be either marginal or simultaneous.")
  }

  if(!object$success) {
    stop("Estimation failed. Cannot display confidence intervals.")
  }

  if(object$method == "Parametric Likelihood Ratio Test" & type == "simultaneous") {
    stop("Simultaneous confidence regions is only implemented for the semi-parametric likelihood ratio model.")
  }

  if(length(conf.level) > 1) {
    message("More than one confidence level given. Using only the first")
    conf.level <- conf.level[1]
  }


  if (type == "marginal") {
    if(conf.level != object$conf.level) {
      #Fix this to be done automatically
      stop("Please refit model with the chosen confidence level and call again.")
    }

    cMat <- matrix(NA, 4, 2)
    cMat[1,] <- object$muDeltaCI
    cMat[2,] <- object$alphaDeltaCI
    cMat[3,] <- log(object$alphaDeltaCI)
    cMat[4,] <- object$DeltaCI

    a <- (1 - object$conf.level)/2
    a <- c(a, 1 - a)
    pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")

    colnames(cMat) <- pct
    rownames(cMat) <- c("Mean difference:", "Odds ratio:", "log Odds ratio:", "Delta")

    print.default(cMat)
  } else {
    message("Calculating joint likelihood surface. This may take some time depending on the resolution.")
    joint <- jointContrastCI(object, muDelta, logORdelta, conf.level, plot, offset, resolution)
  }

  if(type == "marginal") {
    invisible(cMat)
  } else {
    invisible(joint)
  }
}

jointContrastLRT <- function(data, muDelta, alphaDelta) {
  yAlive1 <- data[data$R == 0 & data$A == 1, "Y"]
  yAlive2 <- data[data$R == 1 & data$A == 1, "Y"]
  ELRT <- EL::EL.means(yAlive2, yAlive1, mu = muDelta)

  binom <- TruncComp:::logit.LRT(data, alphaDelta)

  as.numeric(ELRT$statistic + binom)
}

jointContrastCI <- function(m, muDelta = NULL, logORdelta = NULL, conf.level = 0.95, plot=TRUE, offset, resolution) {
  if(!("TruncComp" %in% class(m))) {
    stop("m must be an object of type TruncComp")
  }

  if(is.null(muDelta)) {
    muDelta <- seq(m$muDeltaCI[1] - offset , m$muDeltaCI[2] + offset, length.out = resolution)
  }

  if(is.null(logORdelta)) {
    logORdelta = seq(log(m$alphaDeltaCI[1]) - offset, log(m$alphaDeltaCI[2]) + offset, length.out = resolution)
  }

  matOut <- matrix(NA, length(muDelta), length(logORdelta))
  for(a in 1:length(muDelta)) {
    for(b in 1:length(logORdelta)) {
      matOut[a,b] <- jointContrastLRT(m$data, muDelta[a], logORdelta[b])
    }
  }

  if(plot) {
    fields::image.plot(muDelta, logORdelta, matOut, xlab="Difference in mean", ylab="log Odds ratio")
    points(m$muDelta, log(m$alphaDelta), pch=19, cex=3)
    points(0, 0, cex=3, pch=3)
    contour(muDelta, logORdelta, matOut, add=TRUE, levels=stats::qchisq(conf.level, 2), lwd=2)
  }

  #Also add the contour itself to the return
  list(muDelta = muDelta, logORdelta = logORdelta, surface = matOut)
}



