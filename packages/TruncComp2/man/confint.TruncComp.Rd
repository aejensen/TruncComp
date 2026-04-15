\name{confint.TruncComp2}
\alias{confint.TruncComp2}

\title{Confidence intervals for a TruncComp2 fit}

\usage{
\method{confint}{TruncComp2}(object, type = "marginal", muDelta = NULL,
  logORdelta = NULL, conf.level = object$conf.level, plot = TRUE,
  offset = NULL, resolution = 35, algorithm = c("grid", "optimize"), ...)
}

\description{
Computes marginal confidence intervals or a simultaneous confidence-region
surface for a fitted \code{TruncComp2} object. For successful unadjusted fits,
the method also supports on-demand projected and profile-style intervals for
the derived combined-outcome contrast \code{Delta}.
}

\arguments{
\item{object}{A \code{TruncComp2} object returned by \code{\link{truncComp}}.}
\item{type}{One of \code{"marginal"}, \code{"simultaneous"},
\code{"delta_projected"}, or \code{"delta_profile"}.}
\item{muDelta}{Optional grid values for the mean-difference axis when
\code{type = "simultaneous"}.}
\item{logORdelta}{Optional grid values for the log-odds-ratio axis when
\code{type = "simultaneous"}.}
\item{conf.level}{Confidence level for the interval or contour threshold.}
\item{plot}{Logical; if \code{TRUE}, plot the simultaneous confidence surface.}
\item{offset}{Optional simultaneous-grid expansion. If omitted, a
data-adaptive default is derived from the fitted marginal intervals or
fallback data scale. A single number is applied to both axes; a length-2
vector supplies separate expansions for \code{muDelta} and
\code{logORdelta}.}
\item{resolution}{Number of grid points per axis for the simultaneous surface.}
\item{algorithm}{For \code{type = "delta_projected"} and
\code{type = "delta_profile"}, whether to use the default grid-based
approximation (\code{"grid"}) or the slower direct constrained optimization
alternative (\code{"optimize"}).}
\item{...}{Unused additional arguments.}
}

\value{
If \code{type = "marginal"}, a matrix of reported marginal confidence
intervals, returned invisibly. For successful unadjusted fits, the matrix also
includes \code{Delta (marginal)} and \code{Delta (profile likelihood)} rows.

If \code{type = "simultaneous"}, a list with components \code{muDelta},
\code{logORdelta}, and \code{surface}, returned invisibly.

If \code{type = "delta_projected"}, a one-row matrix for
\code{Delta (projected)}, returned invisibly.

If \code{type = "delta_profile"}, a one-row matrix for
\code{Delta (profile likelihood)}, returned invisibly.
}

\details{
Simultaneous confidence regions and projected \code{Delta} intervals are
implemented for successful unadjusted fits obtained with \code{method = "LRT"}
or \code{method = "SPLRT"}. Adjusted fits support only \code{type = "marginal"}
in the current version.

Marginal confidence intervals are stored on the fitted object. For
\code{type = "marginal"}, the requested \code{conf.level} must therefore match
the fitted \code{object$conf.level}; otherwise the model must be refit.

The package distinguishes four different kinds of confidence output:

\itemize{
  \item{\code{type = "marginal"} returns one-dimensional intervals for the
  reported treatment-effect components. For successful unadjusted fits, the
  returned matrix also includes \code{Delta (marginal)} and
  \code{Delta (profile likelihood)}.}
  \item{\code{type = "simultaneous"} returns the two-dimensional joint surface
  for \code{(muDelta, logORdelta)} and uses the \eqn{\chi^2_2} cutoff for the
  simultaneous region. This is the right object for joint inference on the two
  primary treatment-effect components.}
  \item{\code{type = "delta_projected"} returns the projection of that same
  simultaneous region onto the atom-dependent combined-outcome contrast
  \code{Delta}. Because it inherits the full two-parameter simultaneous
  statement, it is usually more conservative than the one-dimensional profile
  interval \code{DeltaProfileCI}. By default this uses a grid-based projection;
  \code{algorithm = "optimize"} requests the slower direct constrained
  optimization version.}
  \item{\code{type = "delta_profile"} returns the one-dimensional profile-style
  interval for \code{Delta}. By default this uses the same evaluated joint
  surface with a \eqn{\chi^2_1} cutoff; \code{algorithm = "optimize"}
  requests the slower direct constrained profiling version.}
}

For successful unadjusted fits, \code{DeltaProfileCI} stored on the fitted
object is the default grid-based one-dimensional profile approximation for
\code{Delta}. Its target is defined through

\deqn{
W_\Delta(d) = \inf_{\eta : \Delta(\eta) = d} W(\eta)
}

and

\deqn{
CI_{\Delta,\mathrm{profile}} = \{ d : W_\Delta(d) \le q_{\chi^2_1}(\mathrm{conf.level}) \}.
}

By contrast, \code{type = "delta_projected"} returns

\deqn{
CI_{\Delta,\mathrm{proj}} =
\left[
\inf_{(\delta_\mu,\delta_\psi)\in C_{\mathrm{joint}}}\Delta(\delta_\mu,\delta_\psi),
\sup_{(\delta_\mu,\delta_\psi)\in C_{\mathrm{joint}}}\Delta(\delta_\mu,\delta_\psi)
\right],
}

where \eqn{C_{\mathrm{joint}}} is the simultaneous region in
\code{(muDelta, logORdelta)}.

The default \code{"grid"} algorithm for both \code{type = "delta_projected"}
and \code{type = "delta_profile"} reads these intervals off an adaptively
expanded evaluated joint surface. The alternative \code{"optimize"} algorithm
instead solves the corresponding constrained optimization problem directly.

When \code{type = "simultaneous"} and \code{offset} is omitted, the grid
expansion is chosen adaptively from the fitted data. This keeps the default
surface window responsive to the scale of the fitted intervals rather than
using one fixed expansion for every dataset.

When \code{plot = TRUE}, the simultaneous surface is drawn with
\pkg{ggplot2} using \code{theme_minimal()}.
}

\examples{
library(TruncComp2)
f0 <- function(n) stats::rnorm(n, 3, 1)
f1 <- function(n) stats::rnorm(n, 3.5, 1)
d <- simulateTruncatedData(n = 40, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)
fit <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
fit_lrt <- truncComp(Y ~ R, atom = 0, data = d, method = "LRT")

confint(fit, type = "marginal")
confint(fit, type = "delta_projected", plot = FALSE)
confint(fit, type = "delta_profile", plot = FALSE)
confint(fit, type = "delta_profile", algorithm = "optimize", plot = FALSE)
confint(fit_lrt, type = "simultaneous", plot = FALSE, resolution = 5)
confint(fit, type = "simultaneous", plot = FALSE, resolution = 5)
}
