utils::globalVariables(c("mu_delta", "log_or_delta", "statistic"))

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
  logAlphaCI <- suppressWarnings(log(m$alpha_delta_ci))
  logAlpha <- suppressWarnings(log(as.numeric(m$alpha_delta)))

  mu_offset <- jointContrastAxisOffset(
    interval = m$mu_delta_ci,
    center = m$mu_delta,
    fallback = jointContrastMuFallbackOffset(m)
  )
  log_or_offset <- jointContrastAxisOffset(
    interval = logAlphaCI,
    center = logAlpha,
    fallback = jointContrastLogORFallbackOffset(m)
  )

  c(
    mu_delta = mu_offset,
    log_or_delta = log_or_offset
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
  logAlphaCI <- suppressWarnings(log(m$alpha_delta_ci))
  logAlpha <- suppressWarnings(log(as.numeric(m$alpha_delta)))

  list(
    mu_delta = jointContrastAxisBounds(m$mu_delta_ci, m$mu_delta, offsets[1]),
    log_or_delta = jointContrastAxisBounds(logAlphaCI, logAlpha, offsets[2]),
    muDelta = jointContrastAxisBounds(m$mu_delta_ci, m$mu_delta, offsets[1]),
    logORdelta = jointContrastAxisBounds(logAlphaCI, logAlpha, offsets[2])
  )
}

jointContrastPlot <- function(mu_delta, log_or_delta, surface, m, conf.level) {
  plot_data <- expand.grid(
    mu_delta = mu_delta,
    log_or_delta = log_or_delta
  )
  plot_data$statistic <- as.vector(surface)

  plot_obj <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = mu_delta, y = log_or_delta, fill = statistic)
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

  logAlphaEstimate <- suppressWarnings(log(as.numeric(m$alpha_delta)))
  if(is.finite(m$mu_delta) && is.finite(logAlphaEstimate)) {
    plot_obj <- plot_obj + ggplot2::geom_point(
      data = data.frame(mu_delta = m$mu_delta, log_or_delta = logAlphaEstimate),
      ggplot2::aes(x = mu_delta, y = log_or_delta),
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
  buildComponentCIMatrix(object, c("mu_delta", "alpha_delta"))
}

ciColumnLabels <- function(conf.level) {
  a <- (1 - conf.level) / 2
  a <- c(a, 1 - a)
  paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
}

componentCIInfo <- function(object) {
  list(
    mu_delta = list(
      label = "mu_delta",
      interval = object$mu_delta_ci
    ),
    alpha_delta = list(
      label = "alpha_delta",
      interval = object$alpha_delta_ci
    )
  )
}

buildComponentCIMatrix <- function(object, parameter) {
  info <- componentCIInfo(object)
  rows <- lapply(parameter, function(name) info[[name]]$interval)
  cMat <- do.call(rbind, rows)
  rownames(cMat) <- unname(vapply(parameter, function(name) info[[name]]$label, character(1)))
  colnames(cMat) <- ciColumnLabels(object$conf_level)
  cMat
}

deltaIntervalLabel <- function(method) {
  switch(
    method,
    welch = "delta (welch)",
    profile = "delta (profile)",
    projected = "delta (projected)"
  )
}

buildSingleIntervalMatrix <- function(interval, label, conf.level) {
  matrix(
    interval,
    nrow = 1,
    dimnames = list(label, ciColumnLabels(conf.level))
  )
}

#' Confidence intervals for a truncated-comparison fit
#'
#' Computes component confidence intervals, joint confidence-region surfaces,
#' and on-demand `delta` intervals for a fitted `"trunc_comp_fit"` object.
#'
#' @param object A `"trunc_comp_fit"` object returned by [truncComp()].
#' @param parameter Parameter selection for the requested interval. Use
#'   `"mu_delta"` and/or `"alpha_delta"` for the stored component intervals,
#'   `"delta"` for the derived combined-outcome contrast, or `"joint"` for the
#'   two-parameter likelihood-ratio surface.
#' @param method Interval construction for `parameter = "delta"`. One of
#'   `"welch"`, `"profile"`, or `"projected"`.
#' @param mu_delta Optional grid values for the mean-difference axis when
#'   `parameter = "joint"`.
#' @param log_or_delta Optional grid values for the log-odds-ratio axis when
#'   `parameter = "joint"`.
#' @param conf_level Confidence level for the interval or contour threshold.
#' @param plot Logical; if `TRUE`, plot the joint confidence surface.
#' @param offset Optional simultaneous-grid expansion. If omitted, a
#'   data-adaptive default is derived from the fitted marginal intervals or from
#'   fallback data scales. A single number is applied to both axes; a length-2
#'   vector supplies separate expansions for `mu_delta` and `log_or_delta`.
#' @param resolution Number of grid points per axis for the surface-based
#'   surface-based calculations.
#' @param algorithm For `parameter = "delta"` with `method = "projected"` or
#'   `method = "profile"`, whether to use the default grid-based approximation
#'   (`"grid"`) or the slower direct optimization alternative (`"optimize"`).
#' @param ... Unused additional arguments.
#' @return Invisibly returns a printed matrix for the component and `delta`
#'   intervals, or a list with the evaluated joint surface for
#'   `parameter = "joint"`.
#' @details
#' Adjusted fits support only the stored component confidence intervals.
#' Joint regions and `delta` intervals are available only for successful
#' unadjusted `LRT` and `SPLRT` fits.
#' @examples
#' library(TruncComp2)
#' d <- load_trunc_comp2_example()
#' fit <- truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
#' confint(fit)
#' confint(fit, parameter = "joint", plot = FALSE, resolution = 10)
#' confint(fit, parameter = "delta", method = "profile")
#' @rdname confint
#' @export
confint.trunc_comp_fit <- function(object, parameter = c("mu_delta", "alpha_delta"),
                                   method = "welch", mu_delta = NULL, log_or_delta = NULL,
                                   conf_level = object$conf_level, plot = TRUE,
                                   offset = NULL, resolution = 35,
                                   algorithm = c("grid", "optimize"), ...) {
  parameter <- unique(gsub("^muDelta$", "mu_delta", parameter))
  parameter <- unique(gsub("^alphaDelta$", "alpha_delta", parameter))
  parameter <- unique(gsub("^Delta$", "delta", parameter))
  parameter <- unique(match.arg(
    parameter,
    choices = c("mu_delta", "alpha_delta", "delta", "joint"),
    several.ok = TRUE
  ))
  conf.level <- validateConfidenceLevel(conf_level)

  if(!isTRUE(object$success)) {
    stop("Estimation failed. Cannot display confidence intervals.")
  }

  if(length(parameter) > 1 && !all(parameter %in% c("mu_delta", "alpha_delta"))) {
    stop("Multiple parameters are only supported for c(\"mu_delta\", \"alpha_delta\").")
  }

  if(length(parameter) == 1 && identical(parameter, "joint")) {
    if(!is.null(object$adjust)) {
      stop("Joint confidence regions are not implemented for adjusted fits.")
    }

    message("Calculating joint likelihood surface.\nThis may take some time depending on the resolution.")
    joint <- joint_contrast_ci(object, mu_delta, log_or_delta, conf.level, plot, offset, resolution)
    return(invisible(joint))
  }

  if(length(parameter) == 1 && identical(parameter, "delta")) {
    if(!is.null(object$adjust)) {
      stop("delta intervals are not implemented for adjusted fits.")
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

  if(!isTRUE(all.equal(conf.level, object$conf_level, tolerance = sqrt(.Machine$double.eps)))) {
    stop(paste0("Component confidence intervals are stored at the fitted confidence level (",
                format(object$conf_level),
                "). Refit the model with the desired conf_level to change them."))
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
                                     conf.level = m$conf_level,
                                     plot = TRUE, offset = NULL, resolution = 35,
                                     include_delta = FALSE) {
  if(!inherits(m, "trunc_comp_fit")) {
    stop("m must be an object of type trunc_comp_fit")
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
    muDelta <- jointContrastGrid(m$mu_delta_ci, m$mu_delta, offset = offsets[1], resolution = resolution)
  }

  if(is.null(logORdelta)) {
    logAlphaCI <- suppressWarnings(log(m$alpha_delta_ci))
    logAlpha <- suppressWarnings(log(as.numeric(m$alpha_delta)))
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

  out <- list(
    mu_delta = muDelta,
    log_or_delta = logORdelta,
    muDelta = muDelta,
    logORdelta = logORdelta,
    surface = matOut
  )
  if(isTRUE(include_delta)) {
    out$delta_surface <- deltaOut
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
#' @param m A successful unadjusted `"trunc_comp_fit"` object fitted with
#'   `method = "LRT"` or `method = "SPLRT"`.
#' @param mu_delta Optional grid values for the mean-difference axis.
#' @param log_or_delta Optional grid values for the log-odds-ratio axis.
#' @param conf_level Confidence level used for the plotted contour. Defaults to
#'   `m$conf_level`.
#' @param plot Logical; if `TRUE`, plot the surface and contour.
#' @param offset Optional grid expansion. If omitted, a data-adaptive default is
#'   derived from the fitted marginal intervals or fallback data scales. A
#'   single number is applied to both axes; a length-2 vector supplies separate
#'   expansions for `mu_delta` and `log_or_delta`.
#' @param resolution Number of grid points per axis.
#' @return A list with components `mu_delta`, `log_or_delta`, and `surface`,
#'   containing the evaluated grid and the corresponding joint likelihood-ratio
#'   statistics.
#' @details
#' This helper is exported for convenience. Internally it is a
#' thin wrapper around [jointContrastSurfaceData()] with `include_delta = FALSE`.
#' @examples
#' library(TruncComp2)
#' d <- load_trunc_comp2_example()
#' fit <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
#' surface <- joint_contrast_ci(fit, plot = FALSE, resolution = 10)
#' str(surface)
#' @export
joint_contrast_ci <- function(m, mu_delta = NULL, log_or_delta = NULL,
                              conf_level = m$conf_level, plot = TRUE,
                              offset = NULL, resolution = 35) {
  surface <- jointContrastSurfaceData(
    m,
    muDelta = mu_delta,
    logORdelta = log_or_delta,
    conf.level = conf_level,
    plot = plot,
    offset = offset,
    resolution = resolution,
    include_delta = FALSE
  )

  list(
    mu_delta = surface$mu_delta,
    log_or_delta = surface$log_or_delta,
    muDelta = surface$mu_delta,
    logORdelta = surface$log_or_delta,
    surface = surface$surface
  )
}

jointContrastCI <- joint_contrast_ci
