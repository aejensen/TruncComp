\name{jointContrastCI}
\alias{jointContrastCI}

\title{Joint confidence-region surface for an unadjusted fit}

\usage{
jointContrastCI(m, muDelta = NULL, logORdelta = NULL,
  conf.level = m$conf.level, plot = TRUE, offset = NULL, resolution = 35)
}

\description{
Evaluates the simultaneous confidence-region surface for a fitted
\code{TruncComp2} object from the unadjusted parametric or semi-parametric
likelihood-ratio method.
}

\arguments{
\item{m}{A successful unadjusted \code{TruncComp2} object fitted with
\code{method = "LRT"} or \code{method = "SPLRT"}.}
\item{muDelta}{Optional grid values for the mean-difference axis.}
\item{logORdelta}{Optional grid values for the log-odds-ratio axis.}
\item{conf.level}{Confidence level used for the plotted contour. Defaults to
the fitted \code{m$conf.level}.}
\item{plot}{Logical; if \code{TRUE}, plot the surface and contour.}
\item{offset}{Optional grid expansion. If omitted, a data-adaptive default is
derived from the fitted marginal intervals or fallback data scale. A single
number is applied to both axes; a length-2 vector supplies separate
expansions for \code{muDelta} and \code{logORdelta}.}
\item{resolution}{Number of grid points per axis.}
}

\value{
A list with components:
\item{muDelta}{Grid values on the mean-difference axis.}
\item{logORdelta}{Grid values on the log-odds-ratio axis.}
\item{surface}{Matrix of joint likelihood-ratio statistics evaluated on the grid.}
}

\details{
This helper is exported for compatibility, but it is only defined for
successful unadjusted fits. Failed fits and adjusted fits are rejected.

Mathematically, the returned surface is the joint statistic

\deqn{
W(\delta_\mu, \delta_\psi)
}

evaluated over a grid in \eqn{\delta_\mu = muDelta} and
\eqn{\delta_\psi = logORdelta}. The simultaneous confidence region is

\deqn{
C_{\mathrm{joint}}
= \{ (\delta_\mu, \delta_\psi) :
W(\delta_\mu, \delta_\psi) \le q_{\chi^2_2}(\mathrm{conf.level}) \}.
}

This is a two-parameter region. It is different from the one-dimensional
\code{DeltaProfileCI}, whose stored default version is the grid-based
profile approximation on the \code{Delta} scale. It is also different from
\code{DeltaProjectedCI}, which projects this two-dimensional region onto the
\code{Delta} scale.

When \code{offset} is omitted, the default grid window is chosen from the
fitted data rather than from one fixed hard-coded expansion.

When \code{plot = TRUE}, the surface is rendered with \pkg{ggplot2} using
\code{theme_minimal()}.
}

\examples{
library(TruncComp2)
d <- loadTruncComp2Example()
fit <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
fit_lrt <- truncComp(Y ~ R, atom = 0, data = d, method = "LRT")

joint <- jointContrastCI(fit, plot = FALSE, offset = 1, resolution = 5)
joint_lrt <- jointContrastCI(fit_lrt, plot = FALSE, offset = 1, resolution = 5)
str(joint)
str(joint_lrt)
}
