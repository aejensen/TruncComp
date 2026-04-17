# Simulation Study

This directory contains the Monte Carlo study used to build the manuscript
simulation figures and tables.

## Local execution

Run the existing local driver from the repository root:

```sh
Rscript simulation-study/scripts/run-simulation-study.R
```

Useful filters:

```sh
Rscript simulation-study/scripts/run-simulation-study.R \
  --reps=2000 \
  --scenarios=S2,S3 \
  --n=50,100,200 \
  --h=0,1,2 \
  --workers=4
```

## Slurm execution

The Slurm wrapper is an alternative entrypoint that keeps the same cell-level
outputs and aggregation logic, but submits one Slurm array task per pending
cell and then a dependent collector job.

Basic submission:

```sh
Rscript simulation-study/scripts/submit-simulation-study-slurm.R
```

Example with filters and scheduler settings:

```sh
Rscript simulation-study/scripts/submit-simulation-study-slurm.R \
  --reps=2000 \
  --scenarios=S2,S3 \
  --n=50,100,200 \
  --h=0,1,2 \
  --partition=standard \
  --time=04:00:00 \
  --mem=4G \
  --array-parallelism=80
```

Dry run without submitting:

```sh
Rscript simulation-study/scripts/submit-simulation-study-slurm.R \
  --scenarios=S2 \
  --n=50 \
  --h=0 \
  --dry-run
```

The Slurm wrapper writes:

- per-cell outputs under `simulation-study/results/.../cells/`
- the pending-cell manifest as `pending-cells.rds`
- Slurm stdout/stderr logs under `simulation-study/results/.../slurm/`
- the latest submission metadata as `simulation-study/results/.../slurm/latest-submission.txt`

Once the array finishes, the dependent collector job aggregates all available
cell files. If the full design is complete, it also publishes a single
manuscript-ready results object to `manuscript/simulation-study-results.rds`.

## Publishing manuscript-ready simulation results

If the raw per-cell `.rds` outputs were produced on another machine, copy them
into `simulation-study/results/trunccomp2-study/cells/` in this repository and
run:

```sh
Rscript simulation-study/scripts/collect-simulation-study-results.R
```

The collector will infer the study configuration from the copied cell files if
`config.rds` is not present, write the aggregated CSV/RDS outputs under
`simulation-study/results/trunccomp2-study/`, and, when the full design is
present, publish `manuscript/simulation-study-results.rds`.

The manuscript build then regenerates the simulation figures and tables from
that single published file:

```sh
cd manuscript
make assets
make pdf
```
