# TruncComp2
Development version of the R package TruncComp2 for two-sample comparison of truncated continuous outcomes.

The package implements:

- a parametric likelihood-ratio test (`method = "LRT"`)
- a semi-parametric likelihood-ratio test (`method = "SPLRT"`)

The current implementation keeps the same public API for both methods, but now computes both paths internally:

- the parametric method uses `glm` and observed-outcome `lm` fits with ML log-likelihood comparison, plus explicit fallbacks for singular boundary cases
- both methods can adjust for baseline covariates through a shared `adjust = ~ ...` specification, returning conditional treatment effects
- the semi-parametric method uses an internal pure-R empirical-likelihood implementation for both the unadjusted two-sample mean difference and the adjusted observed-outcome treatment coefficient

As a result, the package no longer depends on `EL` or `bbmle` at runtime.

The currently supported scope is:

- one binary treatment indicator coded `0/1`
- one atom value representing the unobserved or undefined outcome
- optional additive baseline-covariate adjustment for both `method = "LRT"` and `method = "SPLRT"` through `adjust = ~ ...`
- adjusted `SPLRT` provides fitted tests and marginal confidence intervals, but not simultaneous confidence regions

To install the development version of TruncComp2 run the following commands from within R

```r
library(remotes)
install_github("aejensen/TruncComp", subdir = "packages/TruncComp2")
```

# Documentation

- Statistical model specification: [MODEL.md](MODEL.md)
- Implementation walkthrough: [IMPLEMENTATION.md](IMPLEMENTATION.md)
- Adjusted semi-parametric methodology note: [ADJUSTED_SPLRT.md](ADJUSTED_SPLRT.md)
- Package-local development guide: [DEVELOPMENT.md](DEVELOPMENT.md)
- Packaged example data loader: `loadTruncComp2Example()`
- Packaged adjusted example data loader: `loadTruncComp2AdjustedExample()`

# Main Interface

The primary entry point is:

```r
truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
truncComp(Y ~ R, atom = 0, data = d, method = "LRT", adjust = ~ age + sex)
truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT", adjust = ~ age + sex)
```

The fitted object reports:

- `muDelta`: difference in means among the observed
- `alphaDelta`: odds ratio of being observed
- `Delta`: combined-outcome mean difference at the fitted atom value
- `DeltaMarginalCI`, `DeltaProfileCI`: the stored unadjusted confidence intervals for `Delta`, with `DeltaProfileCI` computed by the default grid-based profile approximation
- `DeltaProjectedCI`: available on demand through `confint(..., type = "delta_projected")`
- `W`: joint likelihood-ratio test statistic
- `p`: joint p-value

When `adjust` is supplied with either method, `muDelta` and `alphaDelta` are
conditional treatment effects from the adjusted observed-outcome and logistic
submodels. In that adjusted setting, `Delta` is not reported and remains `NA`.

For both `method = "LRT"` and `method = "SPLRT"`, simultaneous
confidence-region surfaces are available for unadjusted fits through
`confint(..., type = "simultaneous")`. If `offset` is omitted, the default
surface window is expanded adaptively from the fitted data, and the plot is
rendered with `ggplot2::theme_minimal()`. Adjusted fits currently support only
marginal intervals.

# Confidence Intervals

`TruncComp2` reports several different inferential objects, and they do not all
answer the same question.

For the two primary treatment-effect components:

- `muDeltaCI` is a one-dimensional confidence interval for the treatment effect
  on the observed-outcome mean.
- `alphaDeltaCI` is a one-dimensional confidence interval for the treatment
  effect on the odds of being observed.

These are componentwise intervals. They describe uncertainty about one
parameter at a time, not about the pair jointly.

For successful unadjusted fits, the package also supports a two-dimensional
simultaneous region in

```text
(muDelta, logORdelta)
```

defined by

```text
C_joint = { (mu, psi) : W(mu, psi) <= qchisq(conf.level, 2) }.
```

This is what `confint(fit, type = "simultaneous")` and `jointContrastCI(fit)`
compute.

For the derived combined-outcome contrast

```text
Delta = [p1 * mu1 + (1 - p1) * atom] - [p0 * mu0 + (1 - p0) * atom],
```

the package distinguishes three different intervals for successful unadjusted
fits:

1. `DeltaMarginalCI`
   This is a Welch interval for the raw combined-outcome mean difference. It is
   fast, model-light, and answers the question "what is the uncertainty in the
   ordinary two-sample mean difference on the combined outcome scale?"

2. `DeltaProjectedCI`
   This is the projection of the two-dimensional simultaneous region onto
   `Delta`:

   ```text
   [ min Delta(mu, psi), max Delta(mu, psi) ] over (mu, psi) in C_joint.
   ```

   By default it is computed numerically from a high-resolution evaluation of
   the same 2D joint surface used for the simultaneous region, then projecting
   the accepted surface points onto the `Delta` scale. It is available on
   demand through `confint(..., type = "delta_projected")`. If you want the
   slower direct constrained-optimization version instead, call
   `confint(..., type = "delta_projected", algorithm = "optimize")`. Because
   it inherits the full two-parameter simultaneous guarantee, it is usually
   more conservative than the profile interval.

3. `DeltaProfileCI`
   This is the one-dimensional profile interval for `Delta`:

   ```text
   C_profile = { d : W_Delta(d) <= qchisq(conf.level, 1) }
   ```

   By default, `TruncComp2` estimates this interval from the same evaluated 2D
   surface but with the `chi^2_1` cutoff, so the stored `DeltaProfileCI` and
   backward-compatible alias `DeltaCI` come from a grid-based profile
   approximation. If you want the slower direct optimization-based version,
   call `confint(..., type = "delta_profile", algorithm = "optimize")`, which
   profiles in `Delta` itself rather than reading the interval off the surface.

Practical guidance:

- Use `muDeltaCI` and `alphaDeltaCI` when the scientific question is about the
  two component effects separately.
- Use `confint(..., type = "simultaneous")` when you want joint inference on the
  pair `(muDelta, logORdelta)`.
- Use `DeltaMarginalCI` when you want a simple descriptive interval for the raw
  combined-outcome mean difference.
- Use `DeltaProfileCI` when `Delta` itself is the inferential target and you
  want the model-based one-dimensional interval. The stored version is the
  default grid-based approximation; use `confint(..., type = "delta_profile",
  algorithm = "optimize")` when you specifically want the direct
  optimization-based profile interval.
- Use `DeltaProjectedCI` when you want the range of `Delta` values compatible
  with the full joint simultaneous region. Use
  `algorithm = "optimize"` only when you explicitly want the slower direct
  optimization-based projection.

# Example
```r
library(TruncComp2)

#Define the two distributions for the observed data
f0 <- function(n) stats::rnorm(n, 3, 1)
f1 <- function(n) stats::rnorm(n, 3.5, 1)

#Define probabilities of being observed
pi0 <- 0.35
pi1 <- 0.6

#Simulate data
d <- TruncComp2::simulateTruncatedData(25, f0, f1, pi0, pi1, atom = 0)

#Estimate parameters using the parametric method
fit_lrt <- truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
summary(fit_lrt)
confint(fit_lrt, type = "marginal")
confint(fit_lrt, type = "delta_projected")
confint(fit_lrt, type = "delta_profile")
confint(fit_lrt, type = "delta_profile", algorithm = "optimize")

#Load the fixed adjusted example and compare unadjusted vs adjusted LRT
d_adjusted <- loadTruncComp2AdjustedExample()
fit_lrt_unadjusted_example <- truncComp(Y ~ R, atom = 0, data = d_adjusted[, c("Y", "R")], method = "LRT")
fit_lrt_adjusted <- truncComp(Y ~ R, atom = 0, data = d_adjusted, method = "LRT", adjust = ~ L)
summary(fit_lrt_unadjusted_example)
summary(fit_lrt_adjusted)

#The same adjusted example also works with the semi-parametric method
fit_splrt_unadjusted_example <- truncComp(Y ~ R, atom = 0, data = d_adjusted[, c("Y", "R")], method = "SPLRT")
fit_splrt_adjusted <- truncComp(Y ~ R, atom = 0, data = d_adjusted, method = "SPLRT", adjust = ~ L)
summary(fit_splrt_unadjusted_example)
summary(fit_splrt_adjusted)

#Estimate parameters using the semi-parametric method
fit_splrt <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
summary(fit_splrt)

#The default interface also accepts atom explicitly and can infer it when y[a == 0]
truncComp.default(d$Y, d$A, d$R, method = "LRT", atom = 0)

#Get simultaneous confidence region for an unadjusted fit
confint(fit_lrt, type = "simultaneous", plot = TRUE, resolution = 10)
confint(fit_splrt, type = "simultaneous", plot = TRUE, resolution = 10)
```

# Development

`TruncComp2` is intended to be developed from its own package subtree. For
ordinary work, no edits outside `packages/TruncComp2/` should be needed.

From the package root, the package-local verification entry point is:

```sh
Rscript tools/check-package.R
```

Manual verification is also available:

```sh
R CMD build .
R CMD check --no-manual --no-build-vignettes TruncComp2_*.tar.gz
```
