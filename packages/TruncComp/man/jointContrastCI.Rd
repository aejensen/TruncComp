\name{jointContrastCI}
\alias{jointContrastCI}

\title{Joint confidence-region surface for the semi-parametric fit}

\usage{
jointContrastCI(m, muDelta = NULL, logORdelta = NULL,
  conf.level = m$conf.level, plot = TRUE, offset = 1, resolution = 35)
}

\description{
Evaluates the simultaneous confidence-region surface for a fitted
\code{TruncComp} object from the semi-parametric likelihood-ratio method.
}

\arguments{
\item{m}{A successful \code{TruncComp} object fitted with
\code{method = "SPLRT"}.}
\item{muDelta}{Optional grid values for the mean-difference axis.}
\item{logORdelta}{Optional grid values for the log-odds-ratio axis.}
\item{conf.level}{Confidence level used for the plotted contour. Defaults to
the fitted \code{m$conf.level}.}
\item{plot}{Logical; if \code{TRUE}, plot the surface and contour.}
\item{offset}{Amount by which the default grid is expanded beyond the marginal
interval or centered estimate.}
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
successful semi-parametric fits. Parametric fits and failed fits are rejected.
}

\examples{
library(TruncComp)
data("TruncCompExample")
fit <- truncComp(Y ~ R, atom = 0, data = TruncCompExample, method = "SPLRT")

joint <- jointContrastCI(fit, plot = FALSE, offset = 1, resolution = 5)
str(joint)
}
