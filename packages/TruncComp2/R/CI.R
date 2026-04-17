utils::globalVariables(c("muDelta", "logORdelta", "statistic"))

validateJointContrastResolution <- function(resolution) {
  if(length(resolution) != 1 || !is.numeric(resolution) || !is.finite(resolution) || resolution < 1) {
    stop("resolution must be a single positive integer.")
  }

  as.integer(resolution)
}

jointContrastAxisBounds <- function(interval = NULL, center = NULL, offset) {
  if(length(offset) != 1 || !is.numeric(offset) || !is.finite(offset) || offset < 0) {
    stop("offset must be a single finite non-negative number.")
  }

  if(!is.null(interval) && length(interval) >= 2 && all(is.finite(interval[1:2]))) {
    return(as.numeric(c(interval[1] - offset, interval[2] + offset)))
  }

  if(length(center) == 1 && is.finite(center)) {
    return(as.numeric(c(center - offset, center + offset)))
  }

  as.numeric(c(-offset, offset))
}

jointContrastGrid <- function(interval = NULL, center = NULL, offset, resolution = 35) {
  resolution <- validateJointContrastResolution(resolution)
  bounds <- jointContrastAxisBounds(interval = interval, center = center, offset = offset)
  seq(bounds[1], bounds[2], length.out = resolution)
}

jointContrastMuFallbackOffset <- function(m) {
  data <- m$data
  y0 <- data$Y[data$R == 0 & data$A == 1]
  y1 <- data$Y[data$R == 1 & data$A == 1]

  v0 <- if(length(y0) > 1) stats::var(y0) else 0
  v1 <- if(length(y1) > 1) stats::var(y1) else 0
  se <- sqrt(v0 / max(length(y0), 1) + v1 / max(length(y1), 1))
  if(is.finite(se) && se > 0) {
    return(2 * se)
  }

  data_scale <- max(abs(data$Y[is.finite(data$Y)]), na.rm = TRUE)
  if(is.finite(data_scale) && data_scale > 0) {
    return(data_scale / 4)
  }

  1 / sqrt(max(nrow(data), 1))
}

jointContrastLogORFallbackOffset <- function(m) {
  data <- m$data
  n0 <- sum(data$R == 0)
  n1 <- sum(data$R == 1)
  k0 <- sum(data$R == 0 & data$A == 1)
  k1 <- sum(data$R == 1 & data$A == 1)

  se <- sqrt(
    1 / (k0 + 0.5) +
      1 / (n0 - k0 + 0.5) +
      1 / (k1 + 0.5) +
      1 / (n1 - k1 + 0.5)
  )
  if(is.finite(se) && se > 0) {
    return(2 * se)
  }

  1 / sqrt(max(nrow(data), 1))
}

jointContrastAxisOffset <- function(interval = NULL, center = NULL, fallback = NULL) {
  if(!is.null(interval) && length(interval) >= 2 && all(is.finite(interval[1:2]))) {
    width <- diff(interval[1:2])
    if(is.finite(width) && width > 0) {
      return(as.numeric(width / 2))
    }
  }

  if(length(fallback) == 1 && is.finite(fallback) && fallback > 0) {
    return(as.numeric(fallback))
  }

  if(length(center) == 1 && is.finite(center) && abs(center) > 0) {
    return(as.numeric(abs(center)))
  }

  1
}

jointContrastDefaultOffsets <- function(m) {
  logAlphaCI <- suppressWarnings(log(m$alphaDeltaCI))
  logAlpha <- suppressWarnings(log(as.numeric(m$alphaDelta)))

  c(
    muDelta = jointContrastAxisOffset(
      interval = m$muDeltaCI,
      center = m$muDelta,
      fallback = jointContrastMuFallbackOffset(m)
    ),
    logORdelta = jointContrastAxisOffset(
      interval = logAlphaCI,
      center = logAlpha,
      fallback = jointContrastLogORFallbackOffset(m)
    )
  )
}

normalizeJointContrastOffsets <- function(m, offset = NULL) {
  if(is.null(offset)) {
    return(jointContrastDefaultOffsets(m))
  }

  if(!is.numeric(offset) || any(!is.finite(offset)) || any(offset < 0) || !(length(offset) %in% c(1, 2))) {
    stop("offset must be NULL, a single finite non-negative number, or a length-2 numeric vector.")
  }

  if(length(offset) == 1) {
    return(rep(as.numeric(offset), 2))
  }

  as.numeric(offset)
}

jointContrastDefaultBounds <- function(m, offset = NULL) {
  offsets <- normalizeJointContrastOffsets(m, offset)
  logAlphaCI <- suppressWarnings(log(m$alphaDeltaCI))
  logAlpha <- suppressWarnings(log(as.numeric(m$alphaDelta)))

  list(
    muDelta = jointContrastAxisBounds(m$muDeltaCI, m$muDelta, offsets[1]),
    logORdelta = jointContrastAxisBounds(logAlphaCI, logAlpha, offsets[2])
  )
}

jointContrastPlot <- function(muDelta, logORdelta, surface, m, conf.level) {
  plot_data <- expand.grid(
    muDelta = muDelta,
    logORdelta = logORdelta
  )
  plot_data$statistic <- as.vector(surface)

  plot_obj <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = muDelta, y = logORdelta, fill = statistic)
  ) +
    ggplot2::geom_raster() +
    ggplot2::geom_contour(
      ggplot2::aes(z = statistic),
      breaks = stats::qchisq(conf.level, 2),
      color = "black",
      linewidth = 0.5
    ) +
    ggplot2::scale_fill_gradientn(colors = rev(grDevices::hcl.colors(128, palette = "YlOrRd", rev = FALSE))) +
    ggplot2::labs(
      x = "Difference in means among the observed",
      y = "log OR of being observed",
      fill = "Joint LR"
    ) +
    ggplot2::theme_minimal()

  logAlphaEstimate <- suppressWarnings(log(as.numeric(m$alphaDelta)))
  if(is.finite(m$muDelta) && is.finite(logAlphaEstimate)) {
    plot_obj <- plot_obj + ggplot2::geom_point(
      data = data.frame(muDelta = m$muDelta, logORdelta = logAlphaEstimate),
      ggplot2::aes(x = muDelta, y = logORdelta),
      inherit.aes = FALSE,
      shape = 16,
      size = 2
    )
  }

  plot_obj
}

validateConfidenceLevel <- function(conf.level) {
  if(!(length(conf.level) == 1 &&
       is.numeric(conf.level) &&
       is.finite(conf.level) &&
       conf.level > 0 &&
       conf.level < 1)) {
    stop("conf.level must be a single number strictly between 0 and 1.")
  }

  conf.level
}

buildMarginalCIMatrix <- function(object) {
  buildComponentCIMatrix(object, c("muDelta", "alphaDelta"))
}

ciColumnLabels <- function(conf.level) {
  a <- (1 - conf.level) / 2
  a <- c(a, 1 - a)
  paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
}

componentCIInfo <- function(object) {
  list(
    muDelta = list(
      label = "Difference in means among the observed:",
      interval = object$muDeltaCI
    ),
    alphaDelta = list(
      label = "Odds ratio of being observed:",
      interval = object$alphaDeltaCI
    )
  )
}

buildComponentCIMatrix <- function(object, parameter) {
  info <- componentCIInfo(object)
  rows <- lapply(parameter, function(name) info[[name]]$interval)
  cMat <- do.call(rbind, rows)
  rownames(cMat) <- unname(vapply(parameter, function(name) info[[name]]$label, character(1)))
  colnames(cMat) <- ciColumnLabels(object$conf.level)
  cMat
}

deltaIntervalLabel <- function(method) {
  switch(
    method,
    welch = "Delta (welch)",
    profile = "Delta (profile)",
    projected = "Delta (projected)"
  )
}

buildSingleIntervalMatrix <- function(interval, label, conf.level) {
  matrix(
    interval,
    nrow = 1,
    dimnames = list(label, ciColumnLabels(conf.level))
  )
}

#' Confidence intervals for a TruncComp2 fit
#'
#' Computes component confidence intervals, joint confidence-region surfaces,
#' and on-demand `Delta` intervals for a fitted `"TruncComp2"` object.
#'
#' @param object A `"TruncComp2"` object returned by [truncComp()].
#' @param parameter Parameter selection for the requested interval. Use
#'   `"muDelta"` and/or `"alphaDelta"` for the stored component intervals,
#'   `"Delta"` for the derived combined-outcome contrast, or `"joint"` for the
#'   two-parameter likelihood-ratio surface.
#' @param method Interval construction for `parameter = "Delta"`. One of
#'   `"welch"`, `"profile"`, or `"projected"`.
#' @param muDelta Optional grid values for the mean-difference axis when
#'   `parameter = "joint"`.
#' @param logORdelta Optional grid values for the log-odds-ratio axis when
#'   `parameter = "joint"`.
#' @param conf.level Confidence level for the interval or contour threshold.
#' @param plot Logical; if `TRUE`, plot the joint confidence surface.
#' @param offset Optional simultaneous-grid expansion. If omitted, a
#'   data-adaptive default is derived from the fitted marginal intervals or from
#'   fallback data scales. A single number is applied to both axes; a length-2
#'   vector supplies separate expansions for `muDelta` and `logORdelta`.
#' @param resolution Number of grid points per axis for the surface-based
#'   surface-based calculations.
#' @param algorithm For `parameter = "Delta"` with `method = "projected"` or
#'   `method = "profile"`, whether to use the default grid-based approximation
#'   (`"grid"`) or the slower direct optimization alternative (`"optimize"`).
#' @param ... Unused additional arguments.
#' @return Invisibly returns a printed matrix for the component and `Delta`
#'   intervals, or a list with the evaluated joint surface for
#'   `parameter = "joint"`.
#' @details
#' Adjusted fits support only the stored component confidence intervals.
#' Joint regions and `Delta` intervals are available only for successful
#' unadjusted `LRT` and `SPLRT` fits.
#' @examples
#' library(TruncComp2)
#' d <- loadTruncComp2Example()
#' fit <- truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
#' confint(fit)
#' confint(fit, parameter = "joint", plot = FALSE, resolution = 10)
#' confint(fit, parameter = "Delta", method = "profile")
#' @rdname confint.TruncComp
#' @export
confint.TruncComp2 <- function(object, parameter = c("muDelta", "alphaDelta"),
                              method = "welch", muDelta = NULL, logORdelta = NULL,
                              conf.level = object$conf.level, plot = TRUE,
                              offset = NULL, resolution = 35, algorithm = c("grid", "optimize"),
                              ...) {
  parameter <- unique(match.arg(
    parameter,
    choices = c("muDelta", "alphaDelta", "Delta", "joint"),
    several.ok = TRUE
  ))
  conf.level <- validateConfidenceLevel(conf.level)

  if(!isTRUE(object$success)) {
    stop("Estimation failed. Cannot display confidence intervals.")
  }

  if(length(parameter) > 1 && !all(parameter %in% c("muDelta", "alphaDelta"))) {
    stop("Multiple parameters are only supported for c(\"muDelta\", \"alphaDelta\").")
  }

  if(length(parameter) == 1 && identical(parameter, "joint")) {
    if(!is.null(object$adjust)) {
      stop("Joint confidence regions are not implemented for adjusted fits.")
    }

    message("Calculating joint likelihood surface.\nThis may take some time depending on the resolution.")
    joint <- jointContrastCI(object, muDelta, logORdelta, conf.level, plot, offset, resolution)
    return(invisible(joint))
  }

  if(length(parameter) == 1 && identical(parameter, "Delta")) {
    if(!is.null(object$adjust)) {
      stop("Delta intervals are not implemented for adjusted fits.")
    }

    method <- match.arg(method, c("welch", "profile", "projected"))
    algorithm <- validateDeltaIntervalAlgorithm(algorithm)

    interval <- switch(
      method,
      welch = delta_welch_interval(object$data$Y, object$data$R, conf.level),
      profile = delta_profile_interval(
        object,
        conf.level = conf.level,
        offset = offset,
        resolution = resolution,
        algorithm = algorithm
      ),
      projected = delta_projected_interval(
        object,
        conf.level = conf.level,
        offset = offset,
        resolution = resolution,
        algorithm = algorithm
      )
    )
    interval_mat <- buildSingleIntervalMatrix(interval, deltaIntervalLabel(method), conf.level)
    print.default(interval_mat)
    return(invisible(interval_mat))
  }

  if(!isTRUE(all.equal(conf.level, object$conf.level, tolerance = sqrt(.Machine$double.eps)))) {
    stop(paste0("Component confidence intervals are stored at the fitted confidence level (",
                format(object$conf.level),
                "). Refit the model with the desired conf.level to change them."))
  }

  cMat <- buildComponentCIMatrix(object, parameter)
  print.default(cMat)
  invisible(cMat)
}

parametricJointReference <- function(data, atom = 0) {
  prepareParametricJointReference(data, atom = atom)
}

jointContrastLRT.parametric.cached <- function(parametricReference, muDelta, logORdelta,
                                               tol = 1e-12) {
  parametricJointCandidate(parametricReference, muDelta, logORdelta, tol = tol)$statistic
}

jointContrastLRT.parametric <- function(data, muDelta, logORdelta) {
  parametricReference <- parametricJointReference(data, atom = resolveDefaultAtom(data$Y, data$A))
  jointContrastLRT.parametric.cached(parametricReference, muDelta, logORdelta)
}

jointContrastLRT.cached <- function(yAlive1, yAlive2, muDelta, logORdelta, logitReference,
                                    atom = 0) {
  splrtReference <- list(
    yAlive0 = yAlive1,
    yAlive1 = yAlive2,
    logitReference = logitReference,
    atom = atom
  )

  splrtJointCandidate(splrtReference, muDelta, logORdelta)$statistic
}

jointContrastLRT <- function(data, muDelta, logORdelta) {
  yAlive1 <- data[data$R == 0 & data$A == 1, "Y"]
  yAlive2 <- data[data$R == 1 & data$A == 1, "Y"]
  logitReference <- logit.prepare(data)
  atom <- resolveDefaultAtom(data$Y, data$A)

  jointContrastLRT.cached(yAlive1, yAlive2, muDelta, logORdelta, logitReference, atom = atom)
}

jointContrastSurfaceData <- function(m, muDelta = NULL, logORdelta = NULL,
                                     conf.level = m$conf.level,
                                     plot = TRUE, offset = NULL, resolution = 35,
                                     include_delta = FALSE) {
  if(!inherits(m, "TruncComp2")) {
    stop("m must be an object of type TruncComp2")
  }

  conf.level <- validateConfidenceLevel(conf.level)

  if(!isTRUE(m$success)) {
    stop("Estimation failed. Cannot calculate simultaneous confidence regions.")
  }

  if(!(identical(m$method, "Semi-empirical Likelihood Ratio Test") ||
       identical(m$method, "Parametric Likelihood Ratio Test"))) {
    stop("Simultaneous confidence regions are only implemented for the parametric and semi-parametric likelihood ratio models.")
  }

  if(!is.null(m$adjust)) {
    stop("Simultaneous confidence regions are not implemented for adjusted fits.")
  }

  resolution <- validateJointContrastResolution(resolution)
  offsets <- normalizeJointContrastOffsets(m, offset)

  if(is.null(muDelta)) {
    muDelta <- jointContrastGrid(m$muDeltaCI, m$muDelta, offset = offsets[1], resolution = resolution)
  }

  if(is.null(logORdelta)) {
    logAlphaCI <- suppressWarnings(log(m$alphaDeltaCI))
    logAlpha <- suppressWarnings(log(as.numeric(m$alphaDelta)))
    logORdelta <- jointContrastGrid(logAlphaCI, logAlpha, offset = offsets[2], resolution = resolution)
  }

  matOut <- matrix(NA, length(muDelta), length(logORdelta))
  deltaOut <- if(isTRUE(include_delta)) matrix(NA, length(muDelta), length(logORdelta)) else NULL
  if(identical(m$method, "Semi-empirical Likelihood Ratio Test")) {
    splrtReference <- prepareSPLRTJointReference(m$data, atom = m$atom)

    for(a in seq_along(muDelta)) {
      for(b in seq_along(logORdelta)) {
        candidate <- splrtJointCandidate(splrtReference, muDelta[a], logORdelta[b])
        matOut[a,b] <- candidate$statistic
        if(isTRUE(include_delta)) {
          deltaOut[a,b] <- candidate$Delta
        }
      }
    }
  } else {
    parametricReference <- parametricJointReference(m$data, atom = m$atom)
    for(a in seq_along(muDelta)) {
      for(b in seq_along(logORdelta)) {
        candidate <- parametricJointCandidate(
          parametricReference,
          muDelta[a],
          logORdelta[b]
        )
        matOut[a,b] <- candidate$statistic
        if(isTRUE(include_delta)) {
          deltaOut[a,b] <- candidate$Delta
        }
      }
    }
  }

  if(plot) {
    print(jointContrastPlot(muDelta, logORdelta, matOut, m, conf.level))
  }

  out <- list(muDelta = muDelta, logORdelta = logORdelta, surface = matOut)
  if(isTRUE(include_delta)) {
    out$deltaSurface <- deltaOut
  }
  out
}

#' Joint confidence-region surface for an unadjusted fit
#'
#' Evaluates the simultaneous confidence-region surface for a fitted
#' `"TruncComp2"` object from the unadjusted parametric or semi-parametric
#' likelihood-ratio method.
#'
#' @param m A successful unadjusted `"TruncComp2"` object fitted with
#'   `method = "LRT"` or `method = "SPLRT"`.
#' @param muDelta Optional grid values for the mean-difference axis.
#' @param logORdelta Optional grid values for the log-odds-ratio axis.
#' @param conf.level Confidence level used for the plotted contour. Defaults to
#'   `m$conf.level`.
#' @param plot Logical; if `TRUE`, plot the surface and contour.
#' @param offset Optional grid expansion. If omitted, a data-adaptive default is
#'   derived from the fitted marginal intervals or fallback data scales. A
#'   single number is applied to both axes; a length-2 vector supplies separate
#'   expansions for `muDelta` and `logORdelta`.
#' @param resolution Number of grid points per axis.
#' @return A list with components `muDelta`, `logORdelta`, and `surface`,
#'   containing the evaluated grid and the corresponding joint likelihood-ratio
#'   statistics.
#' @details
#' This helper is exported for convenience and compatibility. Internally it is a
#' thin wrapper around [jointContrastSurfaceData()] with `include_delta = FALSE`.
#' @examples
#' library(TruncComp2)
#' d <- loadTruncComp2Example()
#' fit <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
#' surface <- jointContrastCI(fit, plot = FALSE, resolution = 10)
#' str(surface)
#' @export
jointContrastCI <- function(m, muDelta = NULL, logORdelta = NULL,
                            conf.level = m$conf.level, plot = TRUE,
                            offset = NULL, resolution = 35) {
  surface <- jointContrastSurfaceData(
    m,
    muDelta = muDelta,
    logORdelta = logORdelta,
    conf.level = conf.level,
    plot = plot,
    offset = offset,
    resolution = resolution,
    include_delta = FALSE
  )

  list(muDelta = surface$muDelta, logORdelta = surface$logORdelta, surface = surface$surface)
}
