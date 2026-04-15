\name{truncComp}
\alias{truncComp}
\alias{truncComp.default}
\title{Two-sample comparison of truncated continuous outcome}
\usage{
truncComp(formula, atom, data, method, conf.level = 0.95, init = NULL)

truncComp.default(y, a, r, method, conf.level = 0.95, init = NULL)
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
observation indicator is reconstructed internally as \code{Y != atom}.}
\item{data}{A data frame containing the variables in \code{formula}.}
\item{y}{Outcome vector for the default interface.}
\item{a}{Binary indicator for whether the continuous outcome is observed.}
\item{r}{Binary treatment indicator for the default interface.}
\item{method}{\code{LRT} for the parametric likelihood-ratio method and
\code{SPLRT} for the semi-parametric likelihood-ratio method.}
\item{conf.level}{The level for the confidence intervals (default = 0.95)}
\item{init}{An optional named list of initial values for the parametric maximum
likelihood optimization. If omitted, the function attempts to construct
moment-based starting values.}
}

\value{
Returns an S3 object of class \code{TruncComp}. The object includes the
following elements:
\item{muDelta}{Estimated difference in means among the observed outcomes.}
\item{muDeltaCI}{Confidence interval for \code{muDelta}.}
\item{alphaDelta}{Estimated odds ratio of being observed.}
\item{alphaDeltaCI}{Confidence interval for \code{alphaDelta}.}
\item{Delta}{Empirical difference in means of the combined observed outcome.}
\item{DeltaCI}{Placeholder for a confidence interval for \code{Delta}; currently
returned as \code{c(NA, NA)}.}
\item{W}{Joint likelihood-ratio test statistic.}
\item{p}{Joint p-value for no treatment effect on either component.}
\item{method}{Description of the estimation method used.}
\item{conf.level}{Confidence level used for reported intervals.}
\item{success}{Logical indicator for whether estimation succeeded.}
\item{error}{Error message when estimation fails.}
\item{init}{Initial values used for the parametric method, if applicable.}
\item{data}{Internal standardized data frame with columns \code{Y}, \code{A},
and \code{R}.}
}

\author{Andreas Kryger Jensen <aeje@sund.ku.dk>}

\seealso{
\code{\link{print.TruncComp}}, \code{\link{summary.TruncComp}}
}

\examples{
library(TruncComp)
d <- simTruncData(n = 100, mu0 = 3, mu1 = 3.5, pi0 = 0.6, pi1 = 0.5)

# Formula interface
truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")

# Default interface
truncComp.default(d$Y, d$A, d$R, method = "SPLRT")
}
