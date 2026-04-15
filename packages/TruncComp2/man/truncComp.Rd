\name{truncComp}
\alias{truncComp}
\alias{truncComp.default}
\title{Two-sample comparison of truncated continuous outcome}
\usage{
truncComp(formula, atom, data, method, conf.level = 0.95, init = NULL,
  adjust = NULL)

truncComp.default(y, a, r, method, conf.level = 0.95, init = NULL,
  adjust = NULL, atom = NULL)
}
\description{
Two-sample comparison of a continuous outcome with a distinguished atom value
that represents an unobserved or undefined outcome. The function implements
either a parametric likelihood-ratio test or a semi-parametric likelihood-ratio
test.
}
\arguments{
\item{formula}{A formula with a continuous outcome on the left-hand side and a
single binary treatment indicator on the right-hand side.}
\item{atom}{A single numeric value used for the special atom outcome. The
observation indicator is reconstructed internally as \code{Y != atom}. For the
default interface, \code{atom} may be omitted only when \code{y[a == 0]} has
exactly one unique finite value, in which case it is inferred and stored on the
fitted object.}
\item{data}{A data frame containing the variables in \code{formula}.}
\item{y}{Outcome vector for the default interface.}
\item{a}{Binary indicator for whether the continuous outcome is observed.}
\item{r}{Binary treatment indicator for the default interface.}
\item{method}{\code{LRT} for the parametric likelihood-ratio method and
\code{SPLRT} for the semi-parametric likelihood-ratio method.}
\item{conf.level}{The level for the confidence intervals (default = 0.95)}
\item{init}{An optional object accepted for compatibility with older
parametric fits. The current model-backed \code{LRT} implementation stores it on
the returned object if supplied, but does not use it for estimation.}
\item{adjust}{Optional covariate adjustment. For the formula interface, this
must be \code{NULL} or a one-sided formula such as \code{~ age + sex}. For the
default interface, this must be \code{NULL}, a data frame, or a matrix of
baseline covariates. The same additive adjustment specification is used in both
the observed-outcome and observation components. For adjusted \code{SPLRT}
fits, only marginal intervals are implemented; simultaneous confidence regions
remain unavailable. The same is true for adjusted \code{LRT} fits.}
}

\value{
Returns an S3 object of class \code{TruncComp2}. The object includes the
following elements:
\item{muDelta}{Estimated difference in means among the observed outcomes.}
\item{muDeltaCI}{Confidence interval for \code{muDelta}.}
\item{alphaDelta}{Estimated odds ratio of being observed.}
\item{alphaDeltaCI}{Confidence interval for \code{alphaDelta}.}
\item{Delta}{Combined-outcome mean difference at the fitted \code{atom} value.}
\item{DeltaMarginalCI}{Welch interval for \code{Delta} in successful
unadjusted fits.}
\item{DeltaProjectedCI}{Projection of the simultaneous confidence region onto
\code{Delta} in successful unadjusted fits. This interval is available on
demand through \code{\link{confint.TruncComp2}} with
\code{type = "delta_projected"} and is not stored on the fitted object by
default. The default returned version is grid-based, and
\code{algorithm = "optimize"} requests the slower direct constrained
optimization alternative.}
\item{DeltaProfileCI}{Profile interval for \code{Delta} in successful
unadjusted fits. The stored version is the default grid-based profile
approximation; \code{\link{confint.TruncComp2}} with
\code{type = "delta_profile", algorithm = "optimize"} requests the slower
direct optimization-based version.}
\item{DeltaCI}{Backward-compatible alias of \code{DeltaProfileCI}. For adjusted
fits, \code{Delta} and all \code{Delta} intervals are returned as \code{NA}
because the reported treatment effects are conditional on the supplied
covariates.}
\item{W}{Joint likelihood-ratio test statistic.}
\item{p}{Joint p-value for no treatment effect on either component.}
\item{method}{Description of the estimation method used.}
\item{conf.level}{Confidence level used for reported intervals.}
\item{success}{Logical indicator for whether estimation succeeded.}
\item{error}{Error message when estimation fails.}
\item{init}{Initial values used for the parametric method, if applicable.}
\item{data}{Internal standardized data frame with columns \code{Y}, \code{A},
and \code{R}, plus any supplied adjustment covariates.}
\item{adjust}{Human-readable adjustment specification stored on the fitted
object, or \code{NULL} for unadjusted fits.}
\item{atom}{The fitted atom value used to define the combined outcome.}
}

\details{
For successful unadjusted fits, the derived combined-outcome contrast is

\deqn{
\Delta =
[p_1 \mu_1 + (1 - p_1)\,atom] -
[p_0 \mu_0 + (1 - p_0)\,atom].
}

The package reports three different confidence intervals for this quantity:

\itemize{
  \item{\code{DeltaMarginalCI}: Welch interval for the raw combined-outcome mean
  difference. This is the simplest descriptive interval.}
  \item{\code{DeltaProjectedCI}: projection of the simultaneous two-parameter
  region for \code{(muDelta, logORdelta)} onto the \code{Delta} scale. This is
  computed only on demand because it is more expensive. The default public
  version is grid-based, with a slower direct constrained-optimization
  alternative available through \code{\link{confint.TruncComp2}}.}
  \item{\code{DeltaProfileCI}: one-dimensional profile interval for
  \code{Delta}. The stored version is the default grid-based profile
  approximation, and \code{DeltaCI} is kept as its backward-compatible alias.
  A slower direct optimization-based alternative is available through
  \code{\link{confint.TruncComp2}} with \code{type = "delta_profile"}.}
}

These intervals share the same point estimate but answer different inferential
questions. \code{DeltaProjectedCI} inherits the full simultaneous two-parameter
statement and is therefore usually more conservative than
\code{DeltaProfileCI}.

For adjusted fits, \code{Delta} and all \code{Delta} intervals remain \code{NA}
because the implemented adjusted treatment effects are conditional regression
coefficients rather than a standardized marginal combined-outcome contrast.
}

\author{Andreas Kryger Jensen <aeje@sund.ku.dk>}

\seealso{
\code{\link{print.TruncComp2}}, \code{\link{summary.TruncComp2}},
\code{\link{confint.TruncComp2}}, \code{\link{jointContrastCI}}
}

\examples{
library(TruncComp2)
f0 <- function(n) stats::rnorm(n, 3, 1)
f1 <- function(n) stats::rnorm(n, 3.5, 1)
d <- simulateTruncatedData(n = 100, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)

# Formula interface
truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
d_adjusted <- loadTruncComp2AdjustedExample()
truncComp(Y ~ R, atom = 0, data = d_adjusted[, c("Y", "R")], method = "LRT")
truncComp(Y ~ R, atom = 0, data = d_adjusted, method = "LRT", adjust = ~ L)
truncComp(Y ~ R, atom = 0, data = d_adjusted, method = "SPLRT", adjust = ~ L)
truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")

# Default interface
truncComp.default(d$Y, d$A, d$R, method = "SPLRT")
truncComp.default(d$Y, d$A, d$R, method = "LRT", atom = 0)
}
