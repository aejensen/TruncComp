jointContrastGrid <- function(interval = NULL, center = NULL, offset = 1, resolution = 35) {
  if(length(offset) != 1 || !is.numeric(offset) || !is.finite(offset) || offset < 0) {
    stop("offset must be a single finite non-negative number.")
  }

  if(length(resolution) != 1 || !is.numeric(resolution) || !is.finite(resolution) || resolution < 1) {
    stop("resolution must be a single positive integer.")
  }
  resolution <- as.integer(resolution)

  if(!is.null(interval) && length(interval) >= 2 && all(is.finite(interval[1:2]))) {
    return(seq(interval[1] - offset, interval[2] + offset, length.out = resolution))
  }

  if(length(center) == 1 && is.finite(center)) {
    return(seq(center - offset, center + offset, length.out = resolution))
  }

  seq(-offset, offset, length.out = resolution)
}

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
    rownames(cMat) <- c("Difference in means among the observed:",
                        "Odds ratio of being observed:",
                        "log Odds ratio of being observed:",
                        "Delta")

    print.default(cMat)
  } else {
    message("Calculating joint likelihood surface.\nThis may take some time depending on the resolution.")
    joint <- jointContrastCI(object, muDelta, logORdelta, conf.level, plot, offset, resolution)
  }

  if(type == "marginal") {
    invisible(cMat)
  } else {
    invisible(joint)
  }
}

jointContrastLRT.cached <- function(yAlive1, yAlive2, muDelta, logORdelta, logitReference) {
  muW <- el_mean_diff_statistic(yAlive2, yAlive1, muDelta)
  binom <- logit.LRT.prepared(logitReference, logORdelta)

  as.numeric(muW + binom)
}

jointContrastLRT <- function(data, muDelta, logORdelta) {
  yAlive1 <- data[data$R == 0 & data$A == 1, "Y"]
  yAlive2 <- data[data$R == 1 & data$A == 1, "Y"]
  logitReference <- logit.prepare(data)

  jointContrastLRT.cached(yAlive1, yAlive2, muDelta, logORdelta, logitReference)
}

jointContrastCI <- function(m, muDelta = NULL, logORdelta = NULL, conf.level = 0.95, plot=TRUE, offset, resolution) {
  if(!inherits(m, "TruncComp")) {
    stop("m must be an object of type TruncComp")
  }

  if(!isTRUE(m$success)) {
    stop("Estimation failed. Cannot calculate simultaneous confidence regions.")
  }

  if(!identical(m$method, "Semi-empirical Likelihood Ratio Test")) {
    stop("Simultaneous confidence regions are only implemented for the semi-parametric likelihood ratio model.")
  }

  yAlive1 <- m$data[m$data$R == 0 & m$data$A == 1, "Y"]
  yAlive2 <- m$data[m$data$R == 1 & m$data$A == 1, "Y"]
  logitReference <- logit.prepare(m$data)

  if(is.null(muDelta)) {
    muDelta <- jointContrastGrid(m$muDeltaCI, m$muDelta, offset = offset, resolution = resolution)
  }

  if(is.null(logORdelta)) {
    logAlphaCI <- suppressWarnings(log(m$alphaDeltaCI))
    logAlpha <- suppressWarnings(log(as.numeric(m$alphaDelta)))
    logORdelta <- jointContrastGrid(logAlphaCI, logAlpha, offset = offset, resolution = resolution)
  }

  matOut <- matrix(NA, length(muDelta), length(logORdelta))
  for(a in seq_along(muDelta)) {
    for(b in seq_along(logORdelta)) {
      matOut[a,b] <- jointContrastLRT.cached(yAlive1, yAlive2, muDelta[a], logORdelta[b], logitReference)
    }
  }

  if(plot) {
    #fields::image.plot(muDelta, logORdelta, matOut, useRaster = TRUE,
    #                   xlab="Mean difference among the observed",
    #                   ylab="log OR of being observed")
    image(muDelta, logORdelta, matOut,
          xlab="Difference in means among the observed",
          ylab="log OR of being observed",
          col = rev(fields::tim.colors(128)), useRaster = TRUE)
    logAlphaEstimate <- suppressWarnings(log(as.numeric(m$alphaDelta)))
    if(is.finite(m$muDelta) && is.finite(logAlphaEstimate)) {
      points(m$muDelta, logAlphaEstimate, pch=19, cex=1)
    }
    #points(0, 0, cex=3, pch=1)
    contour(muDelta, logORdelta, matOut, add=TRUE,
            levels=stats::qchisq(conf.level, 2), lwd=1, labels=conf.level)
  }

  #Also add the contour itself to the return
  list(muDelta = muDelta, logORdelta = logORdelta, surface = matOut)
}
