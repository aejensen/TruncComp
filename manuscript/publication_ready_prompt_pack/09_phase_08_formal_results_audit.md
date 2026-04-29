# Phase 8 - Formal Results Audit

```text
Audit every proposition, lemma, formal definition, displayed null hypothesis,
and theorem-like claim in manuscript/manuscrip.tex.

For each formal result, check:
1. Exact role in the manuscript.
2. Whether all objects are defined before the result.
3. Whether the null hypothesis matches the statistic used.
4. Whether assumptions are sufficient and stated at the right level.
5. Whether degrees of freedom are correct.
6. Whether nuisance parameters, profiling, conditioning, standardization, and
   empirical-likelihood restrictions are handled consistently.
7. Whether the result should stay in the main text, move to an appendix, become
   a remark, or be shortened.
8. Whether the proof establishes exactly the stated conclusion.

At minimum, audit the parametric LRT proposition, the covariate-adjusted
semi-parametric LRT proposition, and the likelihood-ratio intervals proposition.
Also audit formal Bayesian model statements when they function as definitions.
```

## Expected Output Format

```text
Formal result: [label/location]
Purpose:
Statement audit:
Assumption audit:
Alignment with later use:
Proof status:
Recommended revision:
Revised statement, if needed:
```
