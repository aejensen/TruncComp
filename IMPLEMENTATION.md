# Implementation Notes for TruncComp

## Overview

This document describes how the current `TruncComp` implementation is organized in code, with particular focus on the semi-parametric likelihood-ratio procedure and the internal empirical-likelihood engine that replaced the previous dependency on `EL::EL.means`.

The package has two main estimation paths:

- `method = "LRT"`: a fully parametric likelihood-ratio test implemented in [R/LRT.R](/Users/czv146/Documents/GitHub/TruncComp/R/LRT.R)
- `method = "SPLRT"`: a semi-parametric likelihood-ratio test implemented in [R/SPLRT.R](/Users/czv146/Documents/GitHub/TruncComp/R/SPLRT.R), using the internal empirical-likelihood helpers in [R/empiricalLikelihood.R](/Users/czv146/Documents/GitHub/TruncComp/R/empiricalLikelihood.R)

The exported entry point for both methods is [R/truncComp.R](/Users/czv146/Documents/GitHub/TruncComp/R/truncComp.R).

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
- `isDataOkay()` in [R/utility.R](/Users/czv146/Documents/GitHub/TruncComp/R/utility.R) requires at least two observed outcomes in each treatment arm.

If the internal data check fails, `truncComp.default()` now returns an error object immediately rather than continuing into estimation.

## Parametric Path

The parametric likelihood-ratio implementation in [R/LRT.R](/Users/czv146/Documents/GitHub/TruncComp/R/LRT.R) is organized around two nested likelihoods:

- `logLikH0()`: common observation probability and common observed-outcome mean across groups
- `logLikHA()`: treatment effect allowed in both components

The helper `dQoL()` computes the observation-level contribution

```math
f_Y(Y_i \mid \mu_i, \sigma)^{A_i}\pi_i^{A_i}(1-\pi_i)^{1-A_i}.
```

The models are fitted with `bbmle::mle2()`, and the nested-model comparison from `bbmle::anova()` produces:

- the joint statistic `W`
- the p-value `p`

The treatment effects returned to the user are:

- `muDelta`: difference in means among observed outcomes
- `alphaDelta`: odds ratio of being observed

## Semi-Parametric Path

The semi-parametric method in [R/SPLRT.R](/Users/czv146/Documents/GitHub/TruncComp/R/SPLRT.R) decomposes the joint test into two components:

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

The file [R/empiricalLikelihood.R](/Users/czv146/Documents/GitHub/TruncComp/R/empiricalLikelihood.R) contains a focused implementation of the two-sample empirical likelihood for a difference in means. It intentionally implements only the slice of functionality that `TruncComp` needs.

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

The simultaneous-confidence-region code lives in [R/CI.R](/Users/czv146/Documents/GitHub/TruncComp/R/CI.R).

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

When `plot = TRUE`, it renders the heat map and contour using the joint `chisq(df = 2)` threshold.

## Numerical Stability Choices

The current implementation makes the following explicit stability choices:

- exact feasible intervals are kept whenever possible
- epsilon adjustments are used only for root bracketing, not to globally shrink the parameter space
- the empirical-likelihood statistic is clipped to `0` when tiny negative/positive roundoff would otherwise appear at the fitted point
- infeasible mean constraints return `Inf` rather than generating ad hoc approximations

These choices are what allow the implementation to match the previous `EL`-based behavior closely while also handling constant-sample edge cases cleanly.

## Test Coverage

The current test suite in [tests/testthat](/Users/czv146/Documents/GitHub/TruncComp/tests/testthat) covers four layers.

### 1. Low-level numerical tests

These tests validate:

- lambda solving
- infeasible residual configurations
- non-negativity and feasibility behavior of the statistic
- exact point-mass cases

### 2. Frozen upstream parity fixtures

A reference fixture generated from CRAN `EL` 1.3 is stored in:

[tests/testthat/fixtures/el_means_reference.rds](/Users/czv146/Documents/GitHub/TruncComp/tests/testthat/fixtures/el_means_reference.rds)

and can be regenerated with:

[tools/generate-el-fixture.R](/Users/czv146/Documents/GitHub/TruncComp/tools/generate-el-fixture.R)

This lets the tests compare against upstream behavior without requiring `EL` at runtime.

### 3. TruncComp integration tests

These verify that:

- `SPLRT` still reproduces the package example outputs
- marginal confidence intervals still print correctly
- the joint confidence-region surface is finite and well-shaped
- the null and fitted-point joint statistics behave as expected

### 4. Failure-mode and property tests

These cover:

- invalid inputs
- symmetry under swapping the two groups
- containment of the estimate inside the reported confidence interval
- consistency between the reported test statistic and p-value

## Source Files

The main files involved in the current implementation are:

- [R/truncComp.R](/Users/czv146/Documents/GitHub/TruncComp/R/truncComp.R): public entry point and dispatch
- [R/LRT.R](/Users/czv146/Documents/GitHub/TruncComp/R/LRT.R): parametric likelihood-ratio implementation
- [R/SPLRT.R](/Users/czv146/Documents/GitHub/TruncComp/R/SPLRT.R): semi-parametric likelihood-ratio implementation
- [R/empiricalLikelihood.R](/Users/czv146/Documents/GitHub/TruncComp/R/empiricalLikelihood.R): internal empirical-likelihood engine
- [R/CI.R](/Users/czv146/Documents/GitHub/TruncComp/R/CI.R): marginal and simultaneous confidence procedures
- [R/logitFunctions.R](/Users/czv146/Documents/GitHub/TruncComp/R/logitFunctions.R): logistic profile-likelihood helper used for the confidence surface
- [R/utility.R](/Users/czv146/Documents/GitHub/TruncComp/R/utility.R): data validation and error-object construction

Together, these files now give a complete in-package implementation of both testing procedures without relying on the external `EL` package at runtime.
