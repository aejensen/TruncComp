# Confidence Interval Simulation Study

This document outlines a proposed companion Monte Carlo study for confidence
intervals and confidence regions in `TruncComp2`. The existing
`simulation-study/` workflow is primarily a power and null-calibration study.
This companion study should answer a different question: whether the reported
intervals cover their target estimands at the advertised level, and how wide or
computationally costly those intervals are across representative data-generating
mechanisms.

The study should be treated as supplementary methodological validation. It
should not replace the current power study, and it should be smaller than the
power grid because profile and projected `delta` intervals can be expensive.

## Scientific Goals

The study should evaluate finite-sample behavior for:

- Component interval coverage for `mu_delta`, the observed-outcome mean
  difference among non-atom observations.
- Component interval coverage for `alpha_delta`, the odds ratio of being
  observed.
- Joint confidence-region coverage for `(mu_delta, log_or_delta)`.
- Derived combined-outcome contrast coverage for `delta` under the three
  implemented interval constructions: Welch, projected likelihood-ratio, and
  profile likelihood-ratio.
- Interval width, non-finite interval rate, fit failure rate, and runtime.

The main estimand hierarchy should match the package and manuscript:

- `mu_delta = E[Y | A = 1, R = 1] - E[Y | A = 1, R = 0]`.
- `log_or_delta = logit(P(A = 1 | R = 1)) - logit(P(A = 1 | R = 0))`.
- `alpha_delta = exp(log_or_delta)`.
- `delta = {pi1 * mu1 + (1 - pi1) * atom} - {pi0 * mu0 + (1 - pi0) * atom}`.

For non-normal survivor distributions, `mu0` and `mu1` are the true survivor
means under the generating distribution, not fitted normal-model means.

## Relationship To The Existing Power Study

The confidence-interval study should reuse as much of the current power-study
infrastructure as possible:

- Reuse `simulation-study/R/simulation-study.R` scenario definitions where the
  truth is already explicit.
- Reuse the same `simulate_truncated_data()`-style data representation:
  `Y`, `A`, `R`, and a scalar `atom`.
- Reuse local and Slurm execution patterns, with design cells split into
  chunked output files followed by an aggregation step.
- Keep output under a separate result directory, for example
  `simulation-study/results/trunccomp2-ci-study/`, so power and interval
  artifacts do not collide.

The CI study should not run the one-dimensional Wilcoxon and t-test power
comparators except for the Welch `delta` interval, because the scientific target
is interval calibration for the `TruncComp2` inferential objects.

## Proposed Design Grid

Use a reduced grid chosen to stress the interval calculations without repeating
the full power study. The current power-study defaults use 168 cells: six
scenarios, four effect levels, and the hybrid sample-size grid
`n = 50, 100, 150, 200, 300, 400, 500`, with 25,000 repetitions per cell.

Recommended scenarios:

- `S1` normal survivor shift. This is the parametric reference setting.
- `S2` survival-only effect. This stresses `alpha_delta` and joint-region
  behavior when `mu_delta = 0`.
- `S3` antagonistic cancellation. This stresses `delta` intervals when the
  combined-outcome contrast can be close to zero despite component effects. Use
  the revised power-study effect grid
  `pi1 = 0.50, 0.525, 0.55, 0.60`.
- `S5` gamma shape-change with mean shift. This stresses non-normal survivor
  distributions where the semi-parametric method should be less model-dependent.
- `S6` rare-outlier contamination. This stresses robustness and interval width.
  Use the revised power-study effect grid
  `delta_h = 0, 0.15, 0.25, 0.35`.

Recommended sample sizes:

- `n = 50`, small-sample behavior.
- `n = 100`, moderate finite-sample behavior.
- `n = 150`, the same focal sample size used in the power-effect figure.
- `n = 500`, large-sample calibration.

Recommended effect levels:

- `h = 0`, null configuration.
- `h = 2`, moderate non-null configuration.
- `h = 3`, strongest non-null configuration.

Recommended repetitions:

- Use `2,000` repetitions for development and routine checks.
- Use `10,000` repetitions for manuscript-quality tables if computationally
  feasible, because interval coverage is the primary estimand in this study.
- If runtime becomes prohibitive, reduce only the most expensive projected or
  profile `delta` intervals before reducing the cell-level repetition count.

This default grid has `5 * 4 * 3 = 60` cells before considering methods. At
`10,000` repetitions per cell, worst-case Monte Carlo standard error for a 95%
coverage estimate is about `sqrt(0.95 * 0.05 / 10000) = 0.0022`.

## Methods And Intervals

Fit both unadjusted `TruncComp2` methods in each replicate:

- `method = "lrt"`.
- `method = "splrt"`.

For each successful fit, collect:

- `mu_delta_ci` from the fitted object or `confint(fit, parameter = "mu_delta")`.
- `alpha_delta_ci` from the fitted object or
  `confint(fit, parameter = "alpha_delta")`.
- The joint surface from `confint(fit, parameter = "joint", plot = FALSE)`.
- `delta` Welch interval from
  `confint(fit, parameter = "delta", method = "welch")`.
- `delta` projected interval from
  `confint(fit, parameter = "delta", method = "projected",
  algorithm = "optimize")`.
- `delta` profile interval from
  `confint(fit, parameter = "delta", method = "profile",
  algorithm = "optimize")`.

Default interval settings:

- `conf.level = 0.95`.
- Include the Welch `delta` interval as the simple direct baseline.
- Use optimizer-based projected and profile `delta` intervals in the main
  simulation.
- Do not compute grid-based projected/profile `delta` intervals.
- Use `resolution = 35` for the joint confidence-region surface. The
  optimizer-based `delta` intervals do not use the grid resolution.

The study should record interval targets as `delta_welch`,
`delta_projected_optimize`, and `delta_profile_optimize`.

## Truth Calculations

Each design cell must store true values before simulation starts:

- `true_mu_delta = survivor_mean1 - survivor_mean0`.
- `true_log_or_delta = qlogis(pi1) - qlogis(pi0)`.
- `true_alpha_delta = exp(true_log_or_delta)`.
- `true_delta = pi1 * survivor_mean1 + (1 - pi1) * atom -
  {pi0 * survivor_mean0 + (1 - pi0) * atom}`.

Coverage indicators:

- `mu_delta_covered = lower_mu <= true_mu_delta <= upper_mu`.
- `alpha_delta_covered = lower_alpha <= true_alpha_delta <= upper_alpha`.
- `log_or_delta_covered = log(lower_alpha) <= true_log_or_delta <= log(upper_alpha)`,
  when the alpha interval is finite and strictly positive.
- `joint_covered = true pair lies in the evaluated joint confidence region`.
- `delta_welch_covered = lower_delta_welch <= true_delta <= upper_delta_welch`.
- `delta_projected_covered = lower_delta_projected <= true_delta <= upper_delta_projected`.
- `delta_profile_covered = lower_delta_profile <= true_delta <= upper_delta_profile`.

For the joint region, use the fitted surface with the usual chi-square threshold:

```r
threshold <- stats::qchisq(conf.level, df = 2)
joint_covered <- interpolated_or_grid_check(
  true_mu_delta,
  true_log_or_delta,
  surface$mu_delta,
  surface$log_or_delta,
  surface$surface,
  threshold
)
```

The first implementation can use a conservative grid check that reports:

- `TRUE` if the exact grid point or interpolated neighborhood contains the true
  pair below threshold.
- `FALSE` if the true pair lies inside the grid bounds and all nearby values are
  above threshold.
- `NA` if the true pair lies outside the evaluated grid bounds, because that is
  a grid construction failure rather than ordinary non-coverage.

Track the outside-grid rate separately.

## Per-Replicate Data Flow

For each replicate:

1. Draw a truncated two-arm dataset using the scenario-specific `pi0`, `pi1`,
   `f0`, `f1`, and `atom`.
2. Fit `trunc_comp(Y ~ R, atom = atom, data = data, method = "lrt")`.
3. Fit `trunc_comp(Y ~ R, atom = atom, data = data, method = "splrt")`.
4. For each successful fit, compute the requested component, joint, and `delta`
   intervals.
5. Store coverage indicators, interval widths, endpoint finiteness, runtime, and
   any error class/message.

Do not store full datasets for every replicate in the final results. If
debugging is needed, store a reproducible seed and enough metadata to regenerate
the failed replicate.

## Cell-Level Metrics

For each design cell, aggregate by method and interval type:

- `coverage_rate`.
- `coverage_mcse = sqrt(coverage_rate * (1 - coverage_rate) / successful_reps)`.
- `successful_reps`.
- `fit_failure_rate`.
- `interval_failure_rate`.
- `nonfinite_endpoint_rate`.
- `outside_grid_rate` for joint and grid-derived intervals.
- `mean_width`.
- `median_width`.
- `q25_width` and `q75_width`.
- `mean_runtime_sec`.
- `median_runtime_sec`.

For `alpha_delta`, record widths on both scales:

- Odds-ratio width: `upper_alpha - lower_alpha`.
- Log-odds-ratio width: `log(upper_alpha) - log(lower_alpha)`, when finite.

For `joint` regions, record practical size summaries rather than a one-dimensional
width:

- `joint_coverage_rate`.
- `joint_outside_grid_rate`.
- `joint_grid_area`, approximated by the number of accepted cells times grid-cell
  area.
- `joint_mu_range`.
- `joint_log_or_range`.

## Result Objects And Files

Use a separate CI-study version and result directory:

```text
simulation-study/results/trunccomp2-ci-study/
  config.rds
  pending-chunks.rds
  cells/
    S1-h0-n050-chunk001.rds
    S1-h0-n050-chunk002.rds
    ...
  ci-study-cell-metrics.csv
  ci-study-interval-metrics.csv
  ci-study.rds
  slurm/
    latest-submission.txt
    ...
```

Each cell `.rds` should contain:

- `version = "trunccomp2-ci-study-cell-v1"`.
- `chunk_summary`, one row with scenario, sample size, effect level, truth, seed,
  chunk range, repetitions, and completion metadata.
- `interval_metrics`, one row per method and interval target for the chunk.
- `replicate_diagnostics`, compact counts by error type.
- `interval_observations`, compact per-replicate interval diagnostics without
  storing simulated datasets.

The aggregated `ci-study.rds` should contain:

- `version = "trunccomp2-ci-study-v1"`.
- `config`.
- `design`.
- `scenarios`.
- `cell_metrics`.
- `interval_metrics`.
- `complete`, indicating whether every design cell has a result file.

## Proposed Code Structure

Mirror the current power-study organization without mixing outputs:

```text
simulation-study/R/ci-study.R
simulation-study/scripts/run-ci-study.R
simulation-study/scripts/run-ci-study-cell.R
simulation-study/scripts/collect-ci-study-results.R
simulation-study/scripts/submit-ci-study-slurm.R
simulation-study/slurm/run-ci-study-cell.sbatch
simulation-study/slurm/collect-ci-study-results.sbatch
```

The first implementation can avoid adding new Slurm scripts by extending the
existing scripts with a `--study=power|ci` option, but separate CI scripts are
simpler and safer because the metrics and output objects are different.

## Manuscript Outputs

The CI study should produce supplementary assets first:

- A coverage heatmap by scenario, sample size, method, and interval target.
- A compact coverage table at `n = 100` and `n = 150`.
- A width table for `delta` intervals comparing Welch, projected, and profile
  constructions.
- A failure/runtime table for computational transparency.

Recommended manuscript placement:

- Main text: one sentence noting that supplementary simulations assessed
  interval coverage and precision.
- Supplement: full design description, tables, and figures.

The study should only move into the main paper if it reveals a major qualitative
message, such as strong undercoverage of one interval in an important scenario
or a large conservatism gap between projected and profile `delta` intervals.

## Pilot Plan

Before launching the full CI study, run a pilot:

- Scenarios: `S1`, `S3`, `S5`.
- Sample sizes: `n = 50`, `n = 150`.
- Effect levels: `h = 0`, `h = 3`.
- Repetitions: `200`.
- Chunk repetitions: `100` or `200`.
- Joint resolution: `21`.

Pilot success criteria:

- Cell files are written and aggregate cleanly.
- Runtime per replicate is understood for both methods.
- Joint true-pair outside-grid rate is near zero.
- Interval failure rates are explainable.
- Coverage estimates are not used substantively, because `200` repetitions are
  for pipeline validation only.

After the pilot, tune:

- Chunk repetitions.
- Joint resolution or offset expansion for joint surfaces.
- Whether projected/profile optimizer intervals are fast enough for the default
  `chunk_reps = 500`.

## Quality Checks

Add unit or smoke tests for:

- Truth calculations in each scenario.
- Interval extraction from a known successful fit.
- Coverage indicator behavior at lower endpoint, upper endpoint, interior,
  exterior, and non-finite endpoints.
- Aggregation with one complete and one failed method.
- The collector's `complete` flag.

Add reproducibility checks:

- Replicate-level seeds must be deterministic from the cell seed and global
  replicate index, so results do not depend on chunk size.
- Re-running an existing chunk without `overwrite = TRUE` must not replace it.
- The aggregated CSV and RDS outputs must be reproducible from copied chunk
  files.

## Implemented Defaults

The first implementation uses all five recommended scenarios, optimizer-backed
projected/profile `delta` intervals as main metrics, Welch `delta` intervals as
the direct baseline, interpolation for joint-region truth evaluation, and no
adjusted or Bayesian fits.
