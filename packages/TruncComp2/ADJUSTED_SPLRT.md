# Adjusted Semi-Parametric Extension in TruncComp2

This note documents the covariate-adjusted `method = "SPLRT"` extension in
`TruncComp2`. This goes beyond the original manuscript and should be read as a
package-level methodological extension rather than a restatement of the
published theory.

## Target

For adjusted semi-parametric fits, the package targets conditional treatment
effects under an additive no-interaction specification:

```math
\logit P(A = 1 \mid R, L) = \beta_0 + \beta_\delta R + s_\beta(L)
```

```math
E[Y \mid A = 1, R, L] = \mu_0 + \mu_\delta R + s_\mu(L)
```

The reported quantities are:

- `muDelta = \mu_\delta`
- `alphaDelta = exp(\beta_\delta)`

These are conditional treatment effects. For adjusted `SPLRT`, the package sets
`Delta = NA` and `DeltaCI = c(NA, NA)`.

## Joint Test

The adjusted semi-parametric joint null is

```math
H_0:\ \mu_\delta = 0 \text{ and } \beta_\delta = 0.
```

The package constructs the joint statistic additively:

```math
W = W_\mu + W_\alpha.
```

Under regularity, `W_\mu` is treated as asymptotically `\chi^2_1`, `W_\alpha`
is the usual adjusted logistic likelihood-ratio statistic, and the combined
statistic `W` is compared to `\chi^2_2`.

## Continuous Component

Among subjects with `A = 1`, let

```math
X_i = (1, R_i, h(L_i)^T)
```

and let

```math
g_i(\beta) = X_i \left(Y_i - X_i^T \beta\right),
\quad
\beta = (\beta_0, \beta_R, \beta_L).
```

For a fixed candidate treatment effect `\delta`, the package constrains
`\beta_R = \delta`, profiles over the nuisance coefficients
`\beta_{-R} = (\beta_0, \beta_L)`, and computes the empirical-likelihood ratio
from the moment equations

```math
E[g_i(\beta)] = 0.
```

For a fixed coefficient vector `\beta`, the dual problem is solved for
`\lambda` such that

```math
\sum_i \frac{g_i(\beta)}{1 + \lambda^T g_i(\beta)} = 0
```

with the feasibility condition

```math
1 + \lambda^T g_i(\beta) > 0 \quad \text{for all } i.
```

The continuous-component statistic is then

```math
W_\mu(\delta) = 2 \sum_i \log\left(1 + \lambda^T g_i(\beta)\right),
```

after profiling over the nuisance coefficients.

The point estimate `muDelta` is the unconstrained OLS treatment coefficient.
The marginal interval `muDeltaCI` is obtained by inverting

```math
W_\mu(\delta) \le q_{\chi^2_1}(\texttt{conf.level}).
```

## Binary Component

The observation component remains logistic:

- null: `A ~ L`
- alternative: `A ~ R + L`

The package computes

```math
W_\alpha = 2(\ell_\text{alt} - \ell_\text{null}),
```

reports `alphaDelta = exp(\hat\beta_\delta)`, and forms a Wald interval on the
log-odds scale.

## Regularity Policy

This first version prefers clean failure over approximate fallback behavior.
Adjusted `SPLRT` returns a failed fit when:

- the observed-outcome design is rank deficient
- there are too few observed outcomes relative to the adjusted design
- the adjusted logistic model is non-regular or effectively separated
- the empirical-likelihood profile cannot find a feasible constrained solution

Simultaneous confidence regions are not implemented for adjusted `SPLRT` in
this version.
