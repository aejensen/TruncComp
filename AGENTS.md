# Agent Instructions

## Locked IST-3 Application Fit

The Bayesian IST-3 fit used in the manuscript application is locked.

- Do not rerun, refresh, or replace the IST-3 Bayesian fit/cache when editing the manuscript application section unless the user explicitly asks for a rerun/refit/refresh.
- The locked cache is `manuscript/application-data/ist3-bayes-cache.rds`.
- The locked manuscript setting is the package bounded-score logit-normal model with `H = 2`, shared reporting-grid weights for `c(1, 5, 10)`, four chains, 1,000 warmup iterations, and 2,000 sampling iterations per chain.
- Text-only edits to the application section should preserve the cached fit, generated Bayesian tables, generated Bayesian figures, and diagnostics unless explicitly instructed otherwise.
