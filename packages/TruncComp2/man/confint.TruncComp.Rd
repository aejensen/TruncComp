\name{confint.TruncComp2}
\alias{confint.TruncComp2}

\title{Confidence intervals for a TruncComp2 fit}

\usage{
\method{confint}{TruncComp2}(object, type = "marginal", muDelta = NULL,
  logORdelta = NULL, conf.level = object$conf.level, plot = TRUE,
  offset = 1, resolution = 35, ...)
}

\description{
Computes marginal confidence intervals or a simultaneous confidence-region
surface for a fitted \code{TruncComp2} object.
}

\arguments{
\item{object}{A \code{TruncComp2} object returned by \code{\link{truncComp}}.}
\item{type}{Either \code{"marginal"} or \code{"simultaneous"}.}
\item{muDelta}{Optional grid values for the mean-difference axis when
\code{type = "simultaneous"}.}
\item{logORdelta}{Optional grid values for the log-odds-ratio axis when
\code{type = "simultaneous"}.}
\item{conf.level}{Confidence level for the interval or contour threshold.}
\item{plot}{Logical; if \code{TRUE}, plot the simultaneous confidence surface.}
\item{offset}{Amount by which the default simultaneous-confidence grid is
expanded beyond the marginal interval or centered estimate.}
\item{resolution}{Number of grid points per axis for the simultaneous surface.}
\item{...}{Unused additional arguments.}
}

\value{
If \code{type = "marginal"}, a matrix of reported marginal confidence
intervals, returned invisibly. The matrix includes only intervals that are
currently implemented, so the unimplemented \code{DeltaCI} placeholder is not
printed.

If \code{type = "simultaneous"}, a list with components \code{muDelta},
\code{logORdelta}, and \code{surface}, returned invisibly.
}

\details{
Simultaneous confidence regions are only implemented for fits obtained with
\code{method = "SPLRT"} and no covariate adjustment. Adjusted semi-parametric
fits support only \code{type = "marginal"} in the current version.

Marginal confidence intervals are stored on the fitted object. For
\code{type = "marginal"}, the requested \code{conf.level} must therefore match
the fitted \code{object$conf.level}; otherwise the model must be refit.
}

\examples{
library(TruncComp2)
f0 <- function(n) stats::rnorm(n, 3, 1)
f1 <- function(n) stats::rnorm(n, 3.5, 1)
d <- simulateTruncatedData(n = 40, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)
fit <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")

confint(fit, type = "marginal")
confint(fit, type = "simultaneous", plot = FALSE, resolution = 5)
}
