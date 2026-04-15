# Implementation Notes for TruncComp2

## Overview

This document describes how the current `TruncComp2` implementation is organized in code, with particular focus on the semi-parametric likelihood-ratio procedure and the internal empirical-likelihood engine that replaced the previous dependency on `EL::EL.means`.

The package has two main estimation paths:

- `method = "LRT"`: a fully parametric likelihood-ratio test implemented in [R/LRT.R](R/LRT.R)
- `method = "SPLRT"`: a semi-parametric likelihood-ratio test implemented in [R/SPLRT.R](R/SPLRT.R), using the internal empirical-likelihood helpers in [R/empiricalLikelihood.R](R/empiricalLikelihood.R)

The exported entry point for both methods is [R/truncComp.R](R/truncComp.R).

## End-to-End Control Flow

### 1. User interface

The public entry point is:

```r
truncComp(formula, atom, data, method, conf.level = 0.95, init = NULL)
```

The formula interface expects:

- a single outcome on the left-hand side
- one binary treatment indicator on the right-hand side
- one atom value indicating the unobserved or undefined outcome

`truncComp()` extracts the variables from the formula, checks that the treatment is binary, and reconstructs

```text
A = as.numeric(Y != atom)
```

before passing the data to `truncComp.default()`.

### 2. Internal standardized data representation

Inside `truncComp.default()` the data are stored in a standardized frame with columns:

- `Y`: combined observed outcome
- `A`: observation indicator
- `R`: treatment indicator

This representation is used by both the parametric and semi-parametric paths, as well as by the confidence-region code.

### 3. Data validation

There are two layers of validation:

- `truncComp()` ensures that the formula is valid, the treatment variable is binary, the atom is numeric, and that there is at least one observed outcome in each arm.
- `isDataOkay()` in [R/utility.R](R/utility.R) requires at least two observed outcomes in each treatment arm.

If the internal data check fails, `truncComp.default()` now returns an error object immediately rather than continuing into estimation.

Failed fits use the same top-level S3 class as successful fits, so
\code{print()}, \code{summary()}, and \code{confint()} dispatch consistently
even when estimation does not succeed.

For marginal confidence intervals, the fitted object stores the intervals at
its original `conf.level`. The `confint()` method now uses that fitted level by
default and treats a different requested marginal level as a refit request
rather than silently pretending it can recompute the stored intervals.

## Parametric Path

The parametric likelihood-ratio implementation in [R/LRT.R](R/LRT.R) now uses standard model fits as its primary engine:

- `glm(A ~ 1, family = binomial())` and `glm(A ~ R, family = binomial())` for the observation component
- `lm(Y ~ 1)` and `lm(Y ~ R)` on the observed-outcome subset for the Normal component

The test statistic is then assembled from ML log-likelihood differences, with explicit fallback logic for singular boundary cases.

The model is unchanged:

- `A | R` is Bernoulli with one probability per arm under `HA` and a common probability under `H0`
- `Y | A = 1, R` is Normal with one mean per arm under `HA`, a common mean under `H0`, and a common variance in both models

### Bernoulli component

For the observation indicator, the code fits:

- `glm(A ~ 1, family = binomial(), data = data)`
- `glm(A ~ R, family = binomial(), data = data)`

In regular cases it computes `W_A` from the fitted log-likelihoods and derives:

- `alphaDelta = exp(coef(m1)["R"])`
- a Wald interval on the log-odds scale from the fitted coefficient variance

When the observation table hits a boundary, the code preserves exact `0` / `Inf` odds-ratio behavior by falling back to empirical group proportions and returns `NA` intervals.

### Normal component

For the observed outcomes, the code fits:

- `lm(Y ~ 1, data = observed_data)`
- `lm(Y ~ R, data = observed_data)`

In regular cases it computes `W_Y` from `logLik(..., REML = FALSE)` and derives:

- `muDelta` from the treatment coefficient in the observed-outcome model
- a Wald interval from the fitted coefficient variance

When the observed outcomes have zero residual variance, the code falls back to explicit residual-sum-of-squares rules so exact-fit cases still return the intended `0` or `Inf` LR contribution and `NA` intervals.

### Final statistic

The parametric joint test is assembled additively:

```math
W = W_A + W_Y,
```

and the p-value is computed from `chisq(df = 2)`.

The treatment effects returned to the user are:

- `muDelta`: difference in means among observed outcomes
- `alphaDelta`: odds ratio of being observed

The public `init` argument is still accepted for compatibility, but the current model-backed parametric path does not use it for estimation.

## Semi-Parametric Path

The semi-parametric method in [R/SPLRT.R](R/SPLRT.R) decomposes the joint test into two components:

- a mean-difference empirical likelihood ratio for the observed outcomes
- a logistic likelihood-ratio test for the observation indicator

The returned joint statistic is

```math
W = W_\mu + W_\alpha,
```

and the p-value is computed from `chisq(df = 2)`.

### Observed-outcome component

For the continuous part, `SPLRT()` extracts

- `yAlive1`: observed outcomes in the control group
- `yAlive2`: observed outcomes in the treatment group

and passes them to:

```r
el_mean_diff_fit(yAlive2, yAlive1, conf.level = conf.level)
```

This returns:

- `estimate`: the mean difference among observed outcomes
- `conf.int`: the empirical-likelihood confidence interval
- `statistic`: the observed-outcome likelihood-ratio statistic

The result is then mapped into the package-level outputs:

- `muDelta`
- `muDeltaCI`
- `muW`

### Observation component

For the binary part, `SPLRT()` fits:

```r
m0 <- glm(A ~ 1, family = binomial(), data = data)
m1 <- glm(A ~ R, family = binomial(), data = data)
```

and derives:

- `alphaDelta = exp(coef(m1)["R"])`
- `alphaDeltaCI` from `confint(m1)`
- `alphaW` from the logistic likelihood-ratio test `anova(m0, m1, test = "LRT")`

### Final semi-parametric result

The package combines the two components as:

```r
W <- as.numeric(muW + alphaW)
p <- 1 - stats::pchisq(W, 2)
```

The returned object still has the same external shape as before, so the public API did not change when the empirical-likelihood dependency was internalized.

## Internal Empirical-Likelihood Engine

The file [R/empiricalLikelihood.R](R/empiricalLikelihood.R) contains a focused implementation of the two-sample empirical likelihood for a difference in means. It intentionally implements only the slice of functionality that `TruncComp2` needs.

### Design goals

The replacement was designed to:

- remove the runtime dependency on `EL`
- keep the user-facing `SPLRT` results stable
- preserve the `htest`-like shape needed by the current code
- provide a fast statistic-only path for simultaneous confidence-region calculations

### Main internal helpers

The engine is built from the following helpers.

#### `el_validate_mean_diff_inputs()`

Checks that:

- `x` and `y` are numeric
- all values are finite
- both groups contain at least two observations
- `conf.level` is strictly between `0` and `1`

This function is intentionally stricter than the previous external dependency because the package already constrains the supported input shape.

#### `el_one_sample_lambda()`

For a one-sample empirical likelihood problem with residuals `r_i`, this function solves

```math
\sum_i \frac{r_i}{1 + \lambda r_i} = 0.
```

Behavior:

- if all residuals are zero, it returns `lambda = 0`, `statistic = 0`, and `gradient = 0`
- if the residuals do not span zero, the mean constraint is infeasible and the function returns `feasible = FALSE`
- otherwise, it brackets the root using the admissible interval implied by `1 + \lambda r_i > 0` and solves with `stats::uniroot()`

The returned statistic contribution is:

```math
2 \sum_i \log(1 + \lambda r_i).
```

#### `el_mean_diff_delta_range()`

Computes the feasible interval for the mean difference:

```math
[\min(x) - \max(y),\ \max(x) - \min(y)].
```

Outside this interval the profiled empirical-likelihood statistic is treated as `Inf`.

#### `el_mean_diff_theta_range()`

For a fixed mean difference `delta`, the nuisance mean `theta` must satisfy both support constraints:

```math
\theta \in [\max(\min(x), \min(y) + \delta),\ \min(\max(x), \max(y) + \delta)].
```

This range is computed directly from the observed samples.

#### `el_mean_diff_theta_eval()`

Evaluates the profiled statistic at a fixed `theta` and `delta` by forming residuals:

- `x - theta`
- `y - theta + delta`

and solving the two one-sample lambda problems separately.

This returns:

- feasibility
- profiled statistic value
- the derivative of the profiled criterion with respect to `theta`

#### `el_mean_diff_theta()`

Finds the nuisance mean `theta` for a fixed `delta`.

The implementation uses:

- `stats::uniroot()` when the `theta` gradient changes sign over the feasible interval
- `stats::optimize()` on the profiled statistic as a deterministic fallback when the gradient is numerically awkward or does not provide a clean bracket

If the feasible interval collapses to a single point, that point is evaluated directly. This is important for exact point-mass edge cases such as constant separated samples.

#### `el_mean_diff_statistic()`

This is the fast path used by the simultaneous-confidence-region code. It:

- validates the inputs
- returns `0` exactly at the empirical mean difference `mean(x) - mean(y)`
- returns `Inf` outside the feasible delta range
- otherwise profiles out `theta` and returns the `-2 log ELR` statistic

This path avoids constructing confidence intervals when only the statistic is needed.

#### `el_mean_diff_confint()`

Constructs the empirical-likelihood confidence interval by finding the set of `delta` values such that

```math
W_\mu(\delta) \leq \chi^2_1(1 - \alpha).
```

The interval bounds are found by bisection from:

- the estimate toward the lower feasible delta bound
- the estimate toward the upper feasible delta bound

Because the feasible interval is not artificially shrunk, degenerate point-interval cases such as `[0, 0]` and `[-1, -1]` are preserved exactly.

#### `el_mean_diff_fit()`

This is the user-facing internal wrapper used by `SPLRT()`. It returns an `htest`-like object with:

- `estimate`
- `conf.int`
- `statistic`
- `p.value`
- `method`
- `null.value`
- `data.name`

The estimate is deliberately set to the exact sample mean difference:

```r
mean(x) - mean(y)
```

rather than being obtained by a numerical optimizer.

## Confidence-Region Implementation

The simultaneous-confidence-region code lives in [R/CI.R](R/CI.R).

### `jointContrastLRT()`

This function now combines:

- `el_mean_diff_statistic()` for the observed-outcome component
- `logit.LRT()` for the observation component

for any candidate pair:

- `muDelta`
- `logORdelta`

The local variable was renamed from `alphaDelta` to `logORdelta` internally because the function works on the log-odds-ratio scale for the logistic component.

### `jointContrastCI()`

This function evaluates `jointContrastLRT()` on a rectangular grid and returns:

- `muDelta`
- `logORdelta`
- `surface`

Before building the surface, it now checks that the supplied object:

- inherits from `TruncComp2`
- represents a successful fit
- was fitted with `method = "SPLRT"`

The default grids are also guarded so that non-finite odds-ratio confidence
limits fall back to a finite centered grid rather than propagating `log(0)` or
`log(Inf)` into plotting and surface evaluation.

The exported helper now has its own default `conf.level`, `offset`, and
`resolution`, so it can be called directly without routing through
`confint.TruncComp2()`.

When `plot = TRUE`, it renders the heat map and contour using the joint `chisq(df = 2)` threshold.

To avoid repeated work in the inner loop, `jointContrastCI()` caches the
unconstrained logistic fit and the alive-only outcome splits before evaluating
the grid.

For the marginal path, `confint.TruncComp2()` now prints only the implemented
intervals. The stored `DeltaCI` placeholder remains on the fit object for API
compatibility, but it is omitted from the printed matrix until a real interval
is implemented.

## Numerical Stability Choices

The current implementation makes the following explicit stability choices:

- exact feasible intervals are kept whenever possible
- epsilon adjustments are used only for root bracketing, not to globally shrink the parameter space
- the empirical-likelihood statistic is clipped to `0` when tiny negative/positive roundoff would otherwise appear at the fitted point
- infeasible mean constraints return `Inf` rather than generating ad hoc approximations

These choices are what allow the implementation to match the previous `EL`-based behavior closely while also handling constant-sample edge cases cleanly.

## Test Coverage

The current test suite in [tests/testthat](tests/testthat) covers five layers.

### 1. Low-level numerical tests

These tests validate:

- lambda solving
- infeasible residual configurations
- non-negativity and feasibility behavior of the statistic
- exact point-mass cases

### 2. Frozen upstream parity fixtures

A reference fixture generated from CRAN `EL` 1.3 is stored in:

[tests/testthat/fixtures/el_means_reference.rds](tests/testthat/fixtures/el_means_reference.rds)

and can be regenerated with:

[tools/generate-el-fixture.R](tools/generate-el-fixture.R)

This lets the tests compare against upstream behavior without requiring `EL` at runtime.

### 3. Parametric model-fit and optimization cross-checks

The parametric tests now compare the package results against:

- direct `glm` / `lm` references on regular interior datasets
- a separate numerical optimization of the full joint likelihood

This checks both the model-backed implementation details and the underlying LR statistic without depending on a frozen fixture from an older implementation.

### 4. TruncComp2 integration tests

These verify that:

- `LRT` matches direct `glm` / `lm` references on regular cases
- `LRT` boundary cases still return a usable test statistic even when Wald intervals are undefined
- `SPLRT` still reproduces the package example outputs
- marginal confidence intervals still print correctly
- the joint confidence-region surface is finite and well-shaped
- the null and fitted-point joint statistics behave as expected

### 5. Failure-mode and property tests

These cover:

- invalid inputs
- symmetry under swapping the two groups
- containment of the estimate inside the reported confidence interval
- consistency between the reported test statistic and p-value

## Source Files

The main files involved in the current implementation are:

- [R/truncComp.R](R/truncComp.R): public entry point and dispatch
- [R/LRT.R](R/LRT.R): parametric likelihood-ratio implementation
- [R/SPLRT.R](R/SPLRT.R): semi-parametric likelihood-ratio implementation
- [R/empiricalLikelihood.R](R/empiricalLikelihood.R): internal empirical-likelihood engine
- [R/CI.R](R/CI.R): marginal and simultaneous confidence procedures
- [R/logitFunctions.R](R/logitFunctions.R): logistic profile-likelihood helper used for the confidence surface
- [R/utility.R](R/utility.R): data validation and error-object construction

Together, these files now give a complete in-package implementation of both testing procedures without relying on the external `EL` package at runtime.
