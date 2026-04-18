# TruncComp2
Development version of the R package TruncComp2 for two-sample comparison of truncated continuous outcomes.

The package implements:

- a parametric likelihood-ratio test (`method = "lrt"`)
- a semi-parametric likelihood-ratio test (`method = "splrt"`)
- an experimental Bayesian two-part Dirichlet process mixture path via
  `trunc_comp_bayes()`

The current implementation exposes a unified R-idiomatic API for both methods and computes both paths internally:

- the parametric method uses `glm` and observed-outcome `lm` fits with ML log-likelihood comparison, plus explicit fallbacks for singular boundary cases
- both methods can adjust for baseline covariates through a shared `adjust = ~ ...` specification, returning conditional treatment effects
- the semi-parametric method uses an internal pure-R empirical-likelihood implementation for both the unadjusted two-sample mean difference and the adjusted observed-outcome treatment coefficient

As a result, the package no longer depends on `EL` or `bbmle` at runtime.

The currently supported scope is:

- one binary treatment indicator coded `0/1`
- one atom value representing the unobserved or undefined outcome
- optional additive baseline-covariate adjustment for both `method = "lrt"` and `method = "splrt"` through `adjust = ~ ...`
- adjusted `SPLRT` provides fitted tests and component confidence intervals, but not joint confidence regions or `delta` intervals
- the Bayesian path is currently experimental, no-covariate only, and supports
  either real-line non-atom outcomes or strictly positive non-atom outcomes via
  `continuous_support = "positive_real"`

To install the development version of TruncComp2 run the following commands from within R

```r
library(remotes)
install_github("aejensen/TruncComp", subdir = "packages/TruncComp2")
```

# Documentation

- Source-level roxygen comments under `R/` are the source of truth for
  `NAMESPACE` and `man/`. Regenerate them from the package root with:

```sh
Rscript -e 'if(!requireNamespace("roxygen2", quietly = TRUE)) install.packages("roxygen2", repos = "https://cloud.r-project.org"); roxygen2::roxygenise(".")'
```

- Statistical model specification: [MODEL.md](MODEL.md)
- Implementation walkthrough: [IMPLEMENTATION.md](IMPLEMENTATION.md)
- Adjusted semi-parametric methodology note: [ADJUSTED_SPLRT.md](ADJUSTED_SPLRT.md)
- Package-local development guide: [DEVELOPMENT.md](DEVELOPMENT.md)
- Packaged example dataset: `trunc_comp_example`
- Packaged adjusted example dataset: `trunc_comp_adjusted_example`

# Experimental Bayesian Interface

`trunc_comp_bayes()` fits an explicit two-part Bayesian model with:

- one atom probability per treatment arm
- one truncated stick-breaking mixture per treatment arm for the non-atom
  outcomes

The current Bayesian implementation:

- uses packaged `rstan` models under `inst/stan/`
- supports only the no-covariate model
- supports `continuous_support = "real_line"` for Gaussian kernels and
  `continuous_support = "positive_real"` for Gamma kernels on the positive real
  line
- reports posterior summaries and diagnostics rather than p-values or
  likelihood-ratio statistics
- returns the combined-outcome contrast `delta` as the headline summary,
  alongside `delta_atom`, `mu_delta`, and `alpha_delta`
- includes `posterior_density_plot()` for the arm-specific posterior outcome
  densities implied by the fitted two-part DPMM
- includes `posterior_predictive_check()` for `bayesplot`-based posterior
  predictive checks of both the atom and continuous model parts

# Main Interface

The primary entry point is:

```r
trunc_comp(Y ~ R, atom = 0, data = d, method = "lrt")
trunc_comp(Y ~ R, atom = 0, data = d, method = "lrt", adjust = ~ age + sex)
trunc_comp(Y ~ R, atom = 0, data = d, method = "splrt")
trunc_comp(Y ~ R, atom = 0, data = d, method = "splrt", adjust = ~ age + sex)
trunc_comp_bayes(Y ~ R, atom = 0, data = d)
```

The fitted object reports:

- `mu_delta`: difference in means among the observed
- `alpha_delta`: odds ratio of being observed
- `delta`: combined-outcome mean difference at the fitted atom value
- `statistic`: joint likelihood-ratio test statistic
- `p.value`: joint p-value

The fitted object stores the component intervals `mu_delta_ci` and
`alpha_delta_ci`, but it does not store `delta` confidence intervals. Any
interval for `delta` is computed on demand through `confint()`.

When `adjust` is supplied with either method, `mu_delta` and `alpha_delta` are
conditional treatment effects from the adjusted observed-outcome and logistic
submodels. In that adjusted setting, `delta` is not reported and remains `NA`,
and `confint()` rejects `parameter = "delta"` and `parameter = "joint"`.

For `trunc_comp_bayes()`, the returned object instead stores posterior draws
and summaries for:

- `delta_atom`: difference in atom probability
- `mu_delta`: difference in mean among non-atom outcomes
- `alpha_delta`: odds ratio of being observed
- `delta`: combined-outcome contrast

along with Stan diagnostics such as divergences, `Rhat`, and effective sample
sizes.

The Bayesian path also provides:

- `posterior_density_plot(fit)` to visualize the arm-specific posterior mean
  outcome densities, with pointwise credible ribbons and the atom mass shown
  separately in each arm
- `posterior_predictive_check(fit)` to assess the atom and continuous parts of
  the Bayesian model using posterior predictive checks

For positive-support Bayesian fits, the continuous PPC is shown as a density
overlay on `log(Y)` so the visual diagnostic is not distorted by the boundary
at zero.

For both `method = "lrt"` and `method = "splrt"`, joint confidence-region
surfaces are available for unadjusted fits through
`confint(..., parameter = "joint")`. If `offset` is omitted, the default
surface window is expanded adaptively from the fitted data, and the plot is
rendered with `ggplot2::theme_minimal()`. Adjusted fits currently support only
the stored component intervals.

A compact Bayesian plotting workflow on the packaged example data is:

```r
data("trunc_comp_example", package = "TruncComp2")
fit_bayes <- trunc_comp_bayes(
  Y ~ R,
  atom = 0,
  data = trunc_comp_example,
  chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0
)

posterior_density_plot(fit_bayes)

ppc_plots <- posterior_predictive_check(fit_bayes, seed = 1)
ppc_plots$atom
ppc_plots$continuous
```

For strictly positive non-atom outcomes, switch the continuous kernel family:

```r
positive_data <- simulate_truncated_data(
  20,
  f0 = function(n) stats::rgamma(n, shape = 4, rate = 2.5),
  f1 = function(n) stats::rgamma(n, shape = 5, rate = 2.2),
  pi0 = 0.7,
  pi1 = 0.6,
  atom = 0
)

fit_bayes_positive <- trunc_comp_bayes(
  Y ~ R,
  atom = 0,
  data = positive_data,
  continuous_support = "positive_real",
  chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0
)
```

# Confidence Intervals

`TruncComp2` reports several different inferential objects, and they do not all
answer the same question.

For the two primary treatment-effect components, use:

- `confint(fit, parameter = "mu_delta")`
- `confint(fit, parameter = "alpha_delta")`
- `confint(fit)` or `confint(fit, parameter = c("mu_delta", "alpha_delta"))`

These are componentwise intervals. They describe uncertainty about one
parameter at a time, not about the pair jointly.

For successful unadjusted fits, the package also supports a two-dimensional
simultaneous region in

```text
(mu_delta, log_or_delta)
```

defined by

```text
C_joint = { (mu, psi) : W(mu, psi) <= qchisq(conf.level, 2) }.
```

This is what `confint(fit, parameter = "joint")` and
`joint_contrast_surface(fit)`
compute.

For the derived combined-outcome contrast

```text
Delta = [p1 * mu1 + (1 - p1) * atom] - [p0 * mu0 + (1 - p0) * atom],
```

the package supports three interval constructions for successful unadjusted
fits:

1. `confint(fit, parameter = "delta", method = "welch")`
   This is the Welch interval for the raw combined-outcome mean difference. It
   is fast, model-light, and descriptive.

2. `confint(fit, parameter = "delta", method = "projected")`
   This is the projection of the two-dimensional simultaneous region onto the
   `delta` scale:

   ```text
   [ min Delta(mu, psi), max Delta(mu, psi) ] over (mu, psi) in C_joint.
   ```

   The default implementation is grid-based. The slower direct constrained
   optimization alternative remains available through
   `algorithm = "optimize"`.

3. `confint(fit, parameter = "delta", method = "profile")`
   This is the one-dimensional profile interval for `delta`:

   ```text
   C_profile = { d : W_Delta(d) <= qchisq(conf.level, 1) }
   ```

   The default implementation is grid-based. The slower direct
   optimization-based alternative remains available through
   `algorithm = "optimize"`.

Practical guidance:

- Use `confint(fit)` when the scientific question is about the two component
  effects separately.
- Use `confint(fit, parameter = "joint")` when you want joint inference on the
  pair `(mu_delta, log_or_delta)`.
- Use `method = "welch"` when you want a fast descriptive interval for the
  combined outcome scale.
- Use `method = "profile"` when `delta` itself is the primary inferential
  target.
- Use `method = "projected"` when you want the range of `delta` values
  compatible with the full joint simultaneous region.

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
d <- TruncComp2::simulate_truncated_data(25, f0, f1, pi0, pi1, atom = 0)

#Estimate parameters using the parametric method
fit_lrt <- trunc_comp(Y ~ R, atom = 0, data = d, method = "lrt")
summary(fit_lrt)
confint(fit_lrt)
confint(fit_lrt, parameter = "delta", method = "welch")
confint(fit_lrt, parameter = "delta", method = "projected")
confint(fit_lrt, parameter = "delta", method = "profile")
confint(fit_lrt, parameter = "delta", method = "profile", algorithm = "optimize")

#Load the fixed adjusted example and compare unadjusted vs adjusted LRT
data("trunc_comp_adjusted_example", package = "TruncComp2")
fit_lrt_unadjusted_example <- trunc_comp(
  Y ~ R,
  atom = 0,
  data = trunc_comp_adjusted_example[, c("Y", "R")],
  method = "lrt"
)
fit_lrt_adjusted <- trunc_comp(
  Y ~ R,
  atom = 0,
  data = trunc_comp_adjusted_example,
  method = "lrt",
  adjust = ~ L
)
summary(fit_lrt_unadjusted_example)
summary(fit_lrt_adjusted)

#The same adjusted example also works with the semi-parametric method
fit_splrt_unadjusted_example <- trunc_comp(
  Y ~ R,
  atom = 0,
  data = trunc_comp_adjusted_example[, c("Y", "R")],
  method = "splrt"
)
fit_splrt_adjusted <- trunc_comp(
  Y ~ R,
  atom = 0,
  data = trunc_comp_adjusted_example,
  method = "splrt",
  adjust = ~ L
)
summary(fit_splrt_unadjusted_example)
summary(fit_splrt_adjusted)

#Estimate parameters using the semi-parametric method
fit_splrt <- trunc_comp(Y ~ R, atom = 0, data = d, method = "splrt")
summary(fit_splrt)

#The default interface also accepts atom explicitly and can infer it when y[a == 0]
trunc_comp(d$Y, d$A, d$R, method = "lrt", atom = 0)

#Get simultaneous confidence region for an unadjusted fit
confint(fit_lrt, parameter = "joint", plot = TRUE, resolution = 10)
confint(fit_splrt, parameter = "joint", plot = TRUE, resolution = 10)
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
