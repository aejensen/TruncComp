# TruncComp
Development version of the R package TruncComp for two-sample comparison of truncated continuous outcomes.

The package implements:

- a parametric likelihood-ratio test (`method = "LRT"`)
- a semi-parametric likelihood-ratio test (`method = "SPLRT"`)

The current implementation keeps the same public API for both methods, but now computes both paths internally:

- the parametric method uses a closed-form likelihood-ratio calculation built from the Bernoulli and Normal sufficient statistics
- the semi-parametric method uses an internal pure-R empirical-likelihood implementation for the observed-outcome mean difference

As a result, the package no longer depends on `EL` or `bbmle` at runtime.

The currently supported scope is:

- one binary treatment indicator coded `0/1`
- one atom value representing the unobserved or undefined outcome
- no additional covariates in the public interface

To install the development version of TruncComp run the following commands from within R

```r
library(remotes)
install_github("aejensen/TruncComp", subdir = "packages/TruncComp")
```

# Documentation

- Statistical model specification: [MODEL.md](MODEL.md)
- Implementation walkthrough: [IMPLEMENTATION.md](IMPLEMENTATION.md)

# Main Interface

The primary entry point is:

```r
truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
```

The fitted object reports:

- `muDelta`: difference in means among the observed
- `alphaDelta`: odds ratio of being observed
- `W`: joint likelihood-ratio test statistic
- `p`: joint p-value

For `method = "SPLRT"`, simultaneous confidence-region surfaces are available through `confint(..., type = "simultaneous")`.

# Example
```r
library(TruncComp)

#Define the two distributions for the observed data
f0 <- function(n) stats::rnorm(n, 3, 1)
f1 <- function(n) stats::rnorm(n, 3.5, 1)

#Define probabilities of being observed
pi0 <- 0.35
pi1 <- 0.6

#Simulate data
d <- TruncComp::simulateTruncatedData(25, f0, f1, pi0, pi1)

#Estimate parameters using the parametric method
fit_lrt <- truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
summary(fit_lrt)
confint(fit_lrt, type = "marginal")

#Estimate parameters using the semi-parametric method
fit_splrt <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
summary(fit_splrt)

#Get simultaneous confidence region
confint(fit_splrt, type = "simultaneous", plot = TRUE, resolution = 10)
```

# Development

The repository includes a GitHub Actions workflow that runs package install,
the `testthat` suite, and `R CMD check --no-manual --no-build-vignettes` on
changes. That check combination should be treated as the release gate for
package changes.

From the monorepo root, local package verification is:

```sh
cd packages/TruncComp
R CMD build .
R CMD check --no-manual --no-build-vignettes TruncComp_*.tar.gz
```
