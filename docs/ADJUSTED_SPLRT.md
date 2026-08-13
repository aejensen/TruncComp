# Adjusted semiparametric likelihood ratio inference

This note describes the covariate-adjusted `method = "splrt"` implementation.

## Target

The package uses additive models without treatment-by-covariate interactions.

```math
\logit P(A = 1 \mid R, L)
= \beta_0 + \beta_\delta R + s_\beta(L)
```

```math
E[Y \mid A = 1, R, L]
= \mu_0 + \mu_\delta R + s_\mu(L)
```

The fitted object reports:

- `mu_delta`, the conditional treatment coefficient \(\mu_\delta\)
- `alpha_delta`, the conditional odds ratio \(\exp(\beta_\delta)\) for being
  outside the atom

An adjusted fit stores `delta = NA` because these conditional coefficients do
not by themselves define a standardized marginal coded endpoint contrast.
Consequently, `confint()` rejects `parameter = "delta"`. Joint inference for
`(mu_delta, log(alpha_delta))` remains available through
`confint(fit, parameter = "joint")` and `joint_contrast_surface()`.

## Joint test

The joint null is

```math
H_0:\ \mu_\delta = 0 \text{ and } \beta_\delta = 0.
```

The statistic separates into the profiled empirical likelihood contribution
for the conditional non-atom mean and the profiled logistic likelihood
contribution for the probability of being outside the atom.

```math
W = W_\mu + W_A.
```

Under the regularity conditions stated in the manuscript, each contribution
has a one-degree-of-freedom limit, their limiting efficient scores have zero
cross covariance, and \(W\) has a \(\chi^2_2\) null limit.

## Conditional non-atom mean component

Among observations with `A = 1`, let

```math
X_i = (1, R_i, h(L_i)^T)^T
```

and define

```math
g_i(\gamma) = X_i\left(Y_i - X_i^T\gamma\right).
```

For a fixed candidate treatment coefficient \(m\), the package fixes the
coefficient of `R` at \(m\), profiles the remaining coefficients, and computes
the empirical likelihood ratio under

```math
E[g_i(\gamma) \mid A_i = 1] = 0.
```

For a fixed coefficient vector, the dual multiplier solves

```math
\sum_{i:A_i=1}
\frac{g_i(\gamma)}{1 + \lambda^T g_i(\gamma)} = 0
```

subject to

```math
1 + \lambda^T g_i(\gamma) > 0
\quad\text{for every }i\text{ with }A_i=1.
```

The profiled component statistic is

```math
W_\mu(m)
= \inf_{\gamma:\,\gamma_R=m}
2\sum_{i:A_i=1}\log\left(1 + \lambda^T g_i(\gamma)\right).
```

The point estimate `mu_delta` is the unrestricted least-squares treatment
coefficient. Its marginal interval inverts

```math
W_\mu(m) \leq \chi^2_{1,\texttt{conf.level}}.
```

## Probability component

The logistic component compares `A ~ L` with `A ~ R + L`. It reports
`alpha_delta = exp(beta_delta)` and obtains its marginal interval by inverting
the profile likelihood ratio on the log odds ratio scale.

## Joint confidence region

At each candidate pair `(m, b)`, the package fixes the treatment coefficient
in the conditional mean component at `m` and the logistic treatment coefficient
at `b`. It profiles the adjustment coefficients separately in the two
components and adds the resulting profile statistics. The simultaneous region
is

```math
\{(m,b): W_\mu(m) + W_A(b)
\leq \chi^2_{2,\texttt{conf.level}}\}.
```

## Regularity policy

The adjusted fit fails when the conditional non-atom design is rank deficient,
too few non-atom observations remain for the fitted design, the logistic model
is nonregular or separated, or the empirical likelihood profile has no regular
feasible solution.
