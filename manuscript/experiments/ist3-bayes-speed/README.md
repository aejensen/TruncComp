# IST-3 Bounded-Score Speed Experiment

This directory contains external experimental Stan/R workflows for the IST-3
bounded-score Bayesian fit. It does not modify or depend on modified package
Stan code in `packages/TruncComp2`.

Generated CmdStan executables, draws, benchmark CSVs, and diagnostics are written
under `manuscript/build/ist3-bayes-speed/`, which is ignored by the manuscript
build rules.

## Models

- `ist3_bounded_score_aggregated.stan`: H=2 aggregated bounded-score model with
  arm-level atom binomials, sparse arm-score likelihood cells, shared heaping
  weights, ordered component means, and bounded `phi` to avoid invalid beta-CDF
  proposals. For the no-heaping IST-3 experiment, run this with
  `--heaping-grids=1`.
- `ist3_bounded_score_single_beta.stan`: H=1 fallback with the same atom model,
  score support, heaping structure, sparse score-cell likelihood, and reported
  estimands. This removes the two-component mixture geometry that caused repeated
  H=2 divergences.

## Main Commands

```sh
Rscript manuscript/experiments/ist3-bayes-speed/run_ist3_bayes_speed.R --mode=validate --model=h2_ordered
Rscript manuscript/experiments/ist3-bayes-speed/run_ist3_bayes_speed.R --mode=screen --model=h2_ordered
Rscript manuscript/experiments/ist3-bayes-speed/run_ist3_bayes_speed.R --mode=fit --model=h2_ordered

Rscript manuscript/experiments/ist3-bayes-speed/run_ist3_bayes_speed.R \
  --mode=all \
  --model=h2_ordered \
  --heaping-grids=1 \
  --threads-per-chain=1 \
  --candidates=balanced_alpha_stable \
  --init-strategy=control_low_main \
  --final-warmup=300 \
  --final-sampling=300 \
  --seed=20260426

Rscript manuscript/experiments/ist3-bayes-speed/run_ist3_bayes_speed.R --mode=validate --model=h1_single
Rscript manuscript/experiments/ist3-bayes-speed/run_ist3_bayes_speed.R --mode=screen --model=h1_single
Rscript manuscript/experiments/ist3-bayes-speed/run_ist3_bayes_speed.R --mode=fit --model=h1_single --seed=20260420
```

## Current Accepted Runs

The current preferred H=2 no-heaping experimental fit uses the external
aggregated model with `heaping_grids = 1`, the `balanced_alpha_stable` prior,
`adapt_delta = 0.90`, and `control_low_main` deterministic initialization.
This targets the component-allocation mode that caused the earlier key-estimand
R-hat failures.

Manuscript-scale H=2 confirmation:

- 4 chains, 300 warmup + 300 sampling
- `adapt_delta = 0.90`
- elapsed wall time: 183.8 seconds
- divergences: 0
- max treedepth: 6
- minimum BFMI: 0.90
- max key-estimand R-hat: 1.0058
- minimum key-estimand bulk ESS: 975.2
- minimum key-estimand tail ESS: 732.8

A longer H=2 check with 800 warmup + 800 sampling also accepted in 468.3
seconds, with max key-estimand R-hat 1.0026 and minimum bulk ESS 3466.0.

The selected diagnostics and benchmark table are written under:

- `manuscript/build/ist3-bayes-speed/selected-fit-diagnostics-h2_ordered-no_heaping-control_low_main-gap1_1-stickv.txt`
- `manuscript/build/ist3-bayes-speed/selected-fit-diagnostics-w300s300-h2_ordered-no_heaping-control_low_main-gap1_1-stickv.txt`
- `manuscript/build/ist3-bayes-speed/benchmark-results-h2_ordered-no_heaping-control_low_main-gap1_1-stickv.csv`
- `manuscript/build/ist3-bayes-speed/screen-ranked-h2_ordered-no_heaping-control_low_main-gap1_1-stickv.csv`

The accepted H=1 fallback remains available as a simpler sensitivity check:

- `manuscript/build/ist3-bayes-speed/selected-fit-diagnostics-h1_single.txt`
- `manuscript/build/ist3-bayes-speed/benchmark-results-h1_single.csv`
- `manuscript/build/ist3-bayes-speed/screen-ranked-h1_single.csv`
