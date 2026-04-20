# Bayesian Simulation Study

This document outlines a proposed companion simulation study for the experimental
Bayesian pathway in `TruncComp2`, namely `trunc_comp_bayes()`. The existing
`simulation-study/` workflow evaluates frequentist null calibration and power for
the likelihood-ratio procedures. The confidence-interval companion study
evaluates frequentist interval coverage for the likelihood-ratio outputs. This
Bayesian study should answer a separate question: how well the Bayesian two-part
Dirichlet process mixture model recovers the scientific estimands, calibrates its
credible intervals and posterior probabilities, detects lack of fit through
posterior predictive checks, and behaves computationally.

The Bayesian implementation is currently experimental and no-covariate only.
This study should therefore be smaller than the power study and should initially
live in the supplement or repository validation materials.

## Scientific Goals

The study should evaluate:

- Frequentist coverage of Bayesian equal-tail credible intervals for
  `delta_atom`, `mu_delta`, `alpha_delta`, and `delta`.
- Bias and RMSE of posterior point summaries for the same estimands.
- Posterior contraction, measured by credible interval width as sample size
  increases.
- Calibration of directional posterior probabilities, for example
  `Pr(delta > 0 | data)`, `Pr(mu_delta > 0 | data)`, and
  `Pr(alpha_delta > 1 | data)`.
- Posterior predictive check behavior for atom counts and non-atom outcome
  distributions.
- MCMC reliability, including divergences, maximum `Rhat`, effective sample
  sizes, treedepth warnings, and runtime.
- Mixture truncation behavior, including selected mixture size, omitted-tail
  diagnostic summaries, and whether auto-selection changes inference.

The study should not be framed as a power study. Bayesian outputs are posterior
summaries, credible intervals, directional probabilities, predictive checks, and
diagnostics rather than p-values from a prespecified rejection rule.

## Estimands

The study should use the same observed-data estimands as the Bayesian model and
package summaries:

- `rho_0` and `rho_1`: arm-specific atom probabilities.
- `pi_0 = 1 - rho_0` and `pi_1 = 1 - rho_1`: arm-specific probabilities of
  being observed away from the atom.
- `mu_0_c` and `mu_1_c`: arm-specific means among non-atom outcomes.
- `delta_atom = rho_1 - rho_0`.
- `mu_delta = mu_1_c - mu_0_c`.
- `alpha_delta = {pi_1 / (1 - pi_1)} / {pi_0 / (1 - pi_0)}`.
- `delta = {atom * rho_1 + pi_1 * mu_1_c} -
  {atom * rho_0 + pi_0 * mu_0_c}`.

The truth should be computed analytically from each scenario whenever possible.
For mixture scenarios, store both the analytic mixture mean and a code-level
truth check to avoid transcription errors.

## Relationship To Existing Studies

This study should reuse infrastructure but keep outputs separate:

- Reuse scenario definitions from `simulation-study/R/simulation-study.R` when
  their truths are already available.
- Reuse `simulate_truncated_data()` or the same `Y`, `A`, `R`, `atom`
  representation used by the existing studies.
- Use a separate output directory, for example
  `simulation-study/results/trunccomp2-bayes-study/`.
- Keep Bayesian outputs out of the existing power-study CSVs because the metrics
  are not rejection rates.

The Bayesian study can be run after or independently of the frequentist studies.
Its result object should be independent enough that manuscript assets can be
rebuilt without requiring the power-study result object.

## Proposed Design Grid

Use a compact design because each replicate requires Stan sampling.

Recommended scenarios:

- `S1` normal survivor shift. This is the easiest setting and should favor the
  real-line Gaussian mixture.
- `S3` antagonistic cancellation. This checks whether posterior summaries recover
  component effects even when `delta` is near zero.
- `S5` gamma shape-change with mean shift. This checks the positive-support Gamma
  mixture and skewed survivor distributions.
- `S6` rare-outlier contamination. This checks robustness, tail behavior, and
  posterior predictive diagnostics.

Optional stress scenarios:

- A sparse-observation scenario with low `pi0` and `pi1`, included only after the
  main pipeline is stable.
- A strongly multimodal survivor distribution, included only if mixture behavior
  is a central manuscript claim.

Recommended sample sizes:

- `n = 50`, small-sample posterior behavior.
- `n = 100`, moderate finite-sample behavior.
- `n = 250`, the current power-study focal sample size.

Recommended effect levels:

- `h = 0`, null or no-difference configuration.
- `h = 3`, strong non-null configuration.

Recommended repetitions:

- `50` to `100` repetitions for development smoke tests.
- `200` repetitions for pilot manuscript diagnostics.
- `500` to `1000` repetitions only if runtime and diagnostics are acceptable.

The recommended first full design has `4 * 3 * 2 = 24` cells. At `500`
repetitions per cell, a nominal 95% coverage estimate has worst-case Monte Carlo
standard error about `sqrt(0.95 * 0.05 / 500) = 0.0097`.

## Bayesian Fit Settings

Default settings for pilot runs:

- `chains = 4`.
- `iter_warmup = 500`.
- `iter_sampling = 500`.
- `refresh = 0`.
- `control = list(adapt_delta = 0.95, max_treedepth = 12)`.
- `mixture_components = 10`.
- `auto_select_mixture_components = TRUE`.
- `mixture_components_max = 40`.

Default settings for final runs:

- `chains = 4`.
- `iter_warmup = 1000`.
- `iter_sampling = 1000`.
- Keep `adapt_delta = 0.95` unless pilot divergences suggest increasing to
  `0.99`.
- Keep auto-selection enabled unless it creates unacceptable runtime variability.

Support choices:

- Use `continuous_support = "real_line"` for normal and contamination scenarios.
- Use `continuous_support = "positive_real"` for gamma scenarios and any strictly
  positive non-atom outcome scenario.
- Do not force positive-support data through the real-line model in the primary
  study unless a deliberate misspecification comparison is planned.

Prior settings:

- Start with package defaults.
- Record all prior settings in the result object.
- Add sensitivity runs only after the main study is stable. Recommended
  sensitivity targets are concentration parameter priors and kernel scale priors.

## Metrics

For each estimand and fit, collect posterior summaries:

- Posterior median.
- Posterior mean.
- Equal-tail credible interval endpoints.
- Credible interval width.
- Posterior standard deviation.
- Directional posterior probability, where meaningful.

For each estimand and design cell, aggregate:

- Bias of posterior median.
- Bias of posterior mean.
- RMSE of posterior median.
- RMSE of posterior mean.
- Frequentist credible interval coverage.
- Mean credible interval width.
- Median credible interval width.
- Coverage Monte Carlo standard error.

For directional posterior probabilities, aggregate:

- Mean posterior probability.
- Median posterior probability.
- Proportion of replicates with posterior probability greater than `0.90`,
  `0.95`, and `0.975`.
- Under null cells, distribution of posterior probabilities to assess
  calibration and overconfidence.

For posterior predictive checks, aggregate:

- Atom PPC p-value median and interquartile range.
- Continuous PPC p-value median and interquartile range.
- Proportion of PPC p-values below `0.05` or above `0.95`.
- Discrepancy-specific summaries if multiple continuous discrepancies are stored.

For computation and diagnostics, aggregate:

- Fit success rate.
- Stan error rate.
- Divergent transition rate.
- Maximum `Rhat` across monitored parameters.
- Minimum bulk ESS.
- Minimum tail ESS.
- Treedepth saturation rate.
- Median runtime.
- Mean runtime.
- Selected mixture component level.
- Omitted-tail diagnostic median and upper quantiles.

## Coverage Definitions

Credible interval coverage should be evaluated as a frequentist repeated-sampling
property:

- `covered = lower <= true_value <= upper`.
- Missing or non-finite intervals should not be counted as covered.
- Store a separate `interval_failure` indicator when endpoints are missing or
  non-finite.

For `alpha_delta`, evaluate both:

- Coverage on the odds-ratio scale using `alpha_delta`.
- Coverage on the log-odds-ratio scale using `log(alpha_delta)` when posterior
  draws are finite and positive.

For `delta_atom`, remember that `delta_atom = rho_1 - rho_0`, where larger values
mean more atom mass in treatment. Interpret directional posterior probabilities
accordingly.

## Per-Replicate Workflow

For each replicate:

1. Set a deterministic replicate seed derived from the cell seed and replicate
   index.
2. Generate the dataset using the scenario's `pi0`, `pi1`, `f0`, `f1`, and
   `atom`.
3. Choose `continuous_support` from the scenario metadata.
4. Fit `trunc_comp_bayes(Y ~ R, atom = atom, data = data, ...)`.
5. If the fit fails, store the error message and diagnostics available from the
   failed object.
6. If the fit succeeds, extract posterior draws, summary table, credible
   intervals, posterior probabilities, diagnostics, mixture selection metadata,
   and posterior predictive p-values.
7. Compute coverage and error metrics against the stored cell truth.
8. Write compact replicate-level diagnostics or aggregate them in memory.

Do not store every posterior draw in the final aggregate by default. Posterior
draw storage should be optional and limited to a debugging subset because files
can become large quickly.

## Cell-Level Output

Each cell result should contain:

- `version = "trunccomp2-bayes-study-cell-v1"`.
- `cell_summary`, one row with scenario, sample size, effect level, truth, seed,
  repetition count, support choice, and fit settings.
- `estimand_metrics`, one row per estimand with coverage, bias, RMSE, interval
  width, and posterior probability summaries.
- `ppc_metrics`, one row per PPC discrepancy.
- `diagnostic_metrics`, one row per method/settings combination.
- `error_summary`, compact counts by failure type.

The cell file should not contain all generated datasets. It may contain:

- Seeds for failed replicates.
- A small sample of failed replicate metadata.
- Optional saved fits for a small debugging subset when explicitly requested.

## Aggregate Output

Use a separate result tree:

```text
simulation-study/results/trunccomp2-bayes-study/
  config.rds
  pending-cells.rds
  cells/
    S1-h0-n050.rds
    S1-h0-n100.rds
    ...
  bayes-study-cell-metrics.csv
  bayes-study-estimand-metrics.csv
  bayes-study-ppc-metrics.csv
  bayes-study-diagnostic-metrics.csv
  bayes-study.rds
  slurm/
    latest-submission.txt
    ...
```

The aggregate `bayes-study.rds` should contain:

- `version = "trunccomp2-bayes-study-v1"`.
- `config`.
- `design`.
- `scenarios`.
- `cell_metrics`.
- `estimand_metrics`.
- `ppc_metrics`.
- `diagnostic_metrics`.
- `error_summary`.
- `complete`, indicating whether every design cell has a result file.

## Proposed Code Structure

Mirror the current simulation-study structure, but keep Bayesian code separate:

```text
simulation-study/R/bayes-study.R
simulation-study/R/bayes-manuscript-assets.R
simulation-study/scripts/run-bayes-study.R
simulation-study/scripts/run-bayes-study-cell.R
simulation-study/scripts/collect-bayes-study-results.R
simulation-study/scripts/submit-bayes-study-slurm.R
simulation-study/slurm/run-bayes-study-cell.sbatch
simulation-study/slurm/collect-bayes-study-results.sbatch
```

Separate scripts are preferable to extending the power-study scripts because the
Bayesian metrics, runtime profile, and failure handling differ substantially.

## Slurm Strategy

Use one Slurm array task per design cell, as in the power study.

Suggested defaults:

- Request more memory than the power study because Stan fits are heavier.
- Use fewer simultaneous array tasks to avoid oversubscribing cores.
- Pass `cores = chains` to `trunc_comp_bayes()` only if the scheduler allocation
  actually provides those cores.
- Store Stan stdout/stderr per cell.
- Make the collector job dependent on successful or completed array jobs, then
  aggregate all available cell files even if some cells fail.

The submission wrapper should support:

- `--reps`.
- `--scenarios`.
- `--n`.
- `--h`.
- `--chains`.
- `--iter-warmup`.
- `--iter-sampling`.
- `--mixture-components`.
- `--mixture-components-max`.
- `--auto-select-mixture-components=true|false`.
- `--overwrite`.
- `--dry-run`.

## Manuscript Outputs

The first manuscript-facing assets should be supplementary:

- Coverage table for `delta_atom`, `mu_delta`, `alpha_delta`, and `delta`.
- Bias/RMSE table for posterior medians.
- Credible interval width plot by sample size and scenario.
- PPC p-value distribution plot by scenario.
- Diagnostic table showing success rate, divergence rate, maximum `Rhat`,
  minimum ESS, selected mixture size, and runtime.
- Optional prior-sensitivity table if sensitivity runs are performed.

Main-text use should be limited unless the Bayesian model becomes a central
contribution. A concise main-text sentence could state that supplementary
Bayesian simulations assessed posterior calibration, predictive adequacy, and
MCMC reliability under representative scenarios.

## Pilot Plan

Before launching a full Bayesian simulation, run a pilot:

- Scenarios: `S1`, `S5`, `S6`.
- Sample sizes: `n = 50`, `n = 250`.
- Effect levels: `h = 0`, `h = 3`.
- Repetitions: `20` to `50`.
- `chains = 2`.
- `iter_warmup = 250`.
- `iter_sampling = 250`.
- `mixture_components = 5`.
- `auto_select_mixture_components = FALSE` for the first runtime smoke test.

Pilot success criteria:

- All cells complete and aggregate.
- Runtime per replicate is understood.
- Fit failure and diagnostic rates are acceptable or explainable.
- Posterior summaries are finite for all primary estimands.
- PPC p-values can be computed reliably.
- The result object is small enough to store and move.

Second pilot:

- Turn on `auto_select_mixture_components = TRUE`.
- Use `chains = 4`, `iter_warmup = 500`, and `iter_sampling = 500`.
- Check selected mixture sizes, omitted-tail diagnostics, and runtime.

Only after both pilots should the study move to hundreds of repetitions.

## Quality Checks

Add smoke tests or small scripted checks for:

- Truth calculations for every included scenario.
- Support selection for real-line versus positive-real scenarios.
- Extraction of credible intervals from a known successful Bayesian fit.
- Coverage indicator behavior for endpoint, interior, exterior, missing, and
  non-finite intervals.
- Aggregation with successful, failed, and partially diagnostic fits.
- Deterministic cell and replicate seeds.
- Collector behavior when only a subset of cells is complete.

Do not add long Stan runs to routine unit tests. Use tiny smoke tests only, and
keep the full Bayesian simulation behind explicit scripts.

## Interpretation Rules

The Bayesian study should be interpreted carefully:

- Credible interval coverage is a repeated-sampling diagnostic, not the
  definition of Bayesian validity.
- PPC p-values are model-checking summaries, not classical p-values.
- Poor PPC behavior can indicate model mismatch even when posterior intervals
  appear calibrated for means.
- Good calibration in these scenarios does not validate covariate-adjusted
  Bayesian models, because the current implementation is no-covariate only.
- Mixture component counts are approximation diagnostics, not scientific
  estimands.

## Open Decisions

Before implementation, decide:

- Whether the Bayesian model should remain experimental in the manuscript or
  become a supported inferential pathway.
- Whether the main study should include `S3`, given that `delta` cancellation can
  make directional posterior probabilities for `delta` hard to interpret.
- Whether to include explicit model-misspecification scenarios, such as fitting
  real-line kernels to positive skewed data.
- Whether prior sensitivity is required for publication.
- Whether posterior draws should ever be saved, and if so for which debugging or
  manuscript subset.
- Whether to include Bayesian results in the same manuscript appendix as the
  confidence-interval study or in a separate Bayesian validation appendix.
