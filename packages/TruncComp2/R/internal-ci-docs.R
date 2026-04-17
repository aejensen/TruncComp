#' Internal confidence-surface helpers
#'
#' Developer-facing documentation for the helpers that build marginal interval
#' displays, simultaneous confidence grids, and joint likelihood surfaces.
#'
#' @name trunccomp2-ci-helpers
#' @title Internal confidence-surface helpers
#' @description
#' These helpers validate plotting and grid arguments, choose adaptive grid
#' bounds, evaluate joint likelihood-ratio surfaces, and provide compatibility
#' wrappers around the parametric and semi-parametric candidate calculators.
#' @aliases validateJointContrastResolution
#' @aliases jointContrastAxisBounds
#' @aliases jointContrastGrid
#' @aliases jointContrastMuFallbackOffset
#' @aliases jointContrastLogORFallbackOffset
#' @aliases jointContrastAxisOffset
#' @aliases jointContrastDefaultOffsets
#' @aliases normalizeJointContrastOffsets
#' @aliases jointContrastDefaultBounds
#' @aliases jointContrastPlot
#' @aliases validateConfidenceLevel
#' @aliases buildMarginalCIMatrix
#' @aliases parametricJointReference
#' @aliases jointContrastLRT.parametric.cached
#' @aliases jointContrastLRT.parametric
#' @aliases jointContrastLRT.cached
#' @aliases jointContrastLRT
#' @aliases jointContrastSurfaceData
#' @usage
#' validateJointContrastResolution(resolution)
#' jointContrastAxisBounds(interval = NULL, center = NULL, offset)
#' jointContrastGrid(interval = NULL, center = NULL, offset, resolution = 35)
#' jointContrastMuFallbackOffset(m)
#' jointContrastLogORFallbackOffset(m)
#' jointContrastAxisOffset(interval = NULL, center = NULL, fallback = NULL)
#' jointContrastDefaultOffsets(m)
#' normalizeJointContrastOffsets(m, offset = NULL)
#' jointContrastDefaultBounds(m, offset = NULL)
#' jointContrastPlot(muDelta, logORdelta, surface, m, conf.level)
#' validateConfidenceLevel(conf.level)
#' buildMarginalCIMatrix(object)
#' parametricJointReference(data, atom = 0)
#' jointContrastLRT.parametric.cached(parametricReference, muDelta, logORdelta, tol = 1e-12)
#' jointContrastLRT.parametric(data, muDelta, logORdelta)
#' jointContrastLRT.cached(yAlive1, yAlive2, muDelta, logORdelta, logitReference, atom = 0)
#' jointContrastLRT(data, muDelta, logORdelta)
#' jointContrastSurfaceData(
#'   m, muDelta = NULL, logORdelta = NULL, conf.level = m$conf.level,
#'   plot = TRUE, offset = NULL, resolution = 35, include_delta = FALSE
#' )
#' @details
#' ### `validateJointContrastResolution(resolution)`
#'
#' Validates the requested simultaneous-surface resolution. The helper returns
#' the value coerced to integer and errors unless `resolution` is a single
#' positive finite number. Its role is to keep grid generation and the `Delta`
#' surface approximations from silently using malformed resolutions.
#'
#' ### `jointContrastAxisBounds(interval = NULL, center = NULL, offset)`
#'
#' Derives plotting bounds for one joint-inference axis. `interval` is an
#' optional marginal interval, `center` is an optional point estimate, and
#' `offset` is the expansion amount. The helper returns a length-two numeric
#' vector. It errors when `offset` is not scalar, finite, and non-negative. Its
#' role is to translate marginal inference output into an initial simultaneous
#' search box.
#'
#' ### `jointContrastGrid(interval = NULL, center = NULL, offset,
#' resolution = 35)`
#'
#' Builds an equally spaced grid for one axis by combining
#' `validateJointContrastResolution()` with `jointContrastAxisBounds()`. The
#' helper returns the numeric grid vector. Its role is to keep both axes of the
#' simultaneous region synchronized with the same resolution logic.
#'
#' ### `jointContrastMuFallbackOffset(m)` and
#' `jointContrastLogORFallbackOffset(m)`
#'
#' Provide fallback expansions for the mean-difference and log-odds axes when
#' marginal intervals are unavailable or degenerate. `m` is a fitted
#' `"TruncComp2"` object. Each helper returns a scalar positive offset derived
#' from observed-data variability or, failing that, a scale based on the raw
#' data and sample size. Their role is to make simultaneous-region evaluation
#' possible even when the component intervals are `NA`.
#'
#' ### `jointContrastAxisOffset(interval = NULL, center = NULL, fallback = NULL)`
#'
#' Chooses a single axis expansion from the preferred available source: half the
#' width of `interval`, otherwise `fallback`, otherwise the absolute `center`,
#' otherwise `1`. The helper returns a scalar numeric offset. Its role is to
#' encode the package's priority order for default simultaneous-region windows.
#'
#' ### `jointContrastDefaultOffsets(m)` and
#' `normalizeJointContrastOffsets(m, offset = NULL)`
#'
#' `jointContrastDefaultOffsets()` returns the named default axis offsets for the
#' fitted object `m`. `normalizeJointContrastOffsets()` either returns those
#' defaults or validates a user-supplied scalar or length-two override. The
#' second helper errors on negative, non-finite, or wrong-length overrides.
#' Their role is to let public and internal callers share the same windowing
#' rules.
#'
#' ### `jointContrastDefaultBounds(m, offset = NULL)`
#'
#' Combines the fitted marginal intervals, transformed odds-ratio scale, and the
#' normalized offsets into a named list with `muDelta` and `logORdelta` bounds.
#' The helper returns that list directly. Its role is to provide the initial box
#' used by both simultaneous regions and grid-based `Delta` intervals.
#'
#' ### `jointContrastPlot(muDelta, logORdelta, surface, m, conf.level)`
#'
#' Creates the `ggplot2` raster-and-contour visualization for the joint
#' likelihood surface. The grid vectors `muDelta` and `logORdelta`, the matrix
#' `surface`, the fitted object `m`, and the contour level `conf.level` fully
#' determine the result. The helper returns a `ggplot` object. It quietly omits
#' the fitted-point overlay when the estimate is not finite. Its role is to keep
#' all plotting logic separate from the numerical surface evaluation.
#'
#' ### `validateConfidenceLevel(conf.level)`
#'
#' Shared scalar confidence-level validator used across confidence-interval and
#' surface helpers. It returns the supplied value unchanged when it lies strictly
#' between `0` and `1`, otherwise it errors. Its role is to make confidence-level
#' handling consistent across the package.
#'
#' ### `buildMarginalCIMatrix(object)`
#'
#' Converts a successful `"TruncComp2"` object into the printed component
#' confidence-interval matrix used by `confint()` when `parameter` selects the
#' stored component intervals. The helper returns a numeric matrix with
#' human-readable row and column names. Its role is to centralize the display
#' formatting of the stored one-parameter intervals.
#'
#' ### `parametricJointReference(data, atom = 0)`
#'
#' Compatibility wrapper around `prepareParametricJointReference()`. It returns
#' the cached parametric joint-reference list used by the surface and `Delta`
#' helpers. Its role is to keep older helper names available inside the package
#' while the actual implementation lives in `delta.R`.
#'
#' ### `jointContrastLRT.parametric.cached(parametricReference, muDelta,
#' logORdelta, tol = 1e-12)` and
#' `jointContrastLRT.parametric(data, muDelta, logORdelta)`
#'
#' Return the parametric joint likelihood-ratio statistic at fixed candidate
#' contrasts, either from a cached reference object or from raw standardized
#' data. The helpers return a scalar statistic and treat infeasible constrained
#' fits as `Inf`. Their role is to expose the parametric joint objective in a
#' minimal form for earlier interfaces and validation utilities.
#'
#' ### `jointContrastLRT.cached(yAlive1, yAlive2, muDelta, logORdelta,
#' logitReference, atom = 0)` and `jointContrastLRT(data, muDelta, logORdelta)`
#'
#' Provide the analogous semi-parametric joint statistic, again either from
#' cached components or from raw data. The helpers return a scalar joint
#' likelihood-ratio statistic. Their role is to bridge the empirical-likelihood
#' continuous component and the profiled logistic component for simultaneous
#' confidence-region evaluation.
#'
#' ### `jointContrastSurfaceData(m, muDelta = NULL, logORdelta = NULL,
#' conf.level = m$conf.level, plot = TRUE, offset = NULL, resolution = 35,
#' include_delta = FALSE)`
#'
#' Core surface evaluator used by [jointContrastCI()] and the grid-based
#' `Delta` intervals. `m` must be a successful unadjusted `"TruncComp2"` fit;
#' the remaining arguments define or derive the evaluation grid, whether to plot,
#' and whether the derived `Delta` values should also be stored. The helper
#' returns a list containing the axis vectors, the surface matrix, and optionally
#' `deltaSurface`. It errors on failed fits, adjusted fits, unsupported methods,
#' or invalid confidence levels. Its role is to centralize all joint-surface
#' evaluation logic so the public helpers and the internal `Delta` routines stay
#' aligned.
#'
#' @seealso [confint.TruncComp2()], [jointContrastCI()], [delta_surface_for_inference()]
#' @keywords internal
NULL
