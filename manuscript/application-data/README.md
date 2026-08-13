# Application Data: IST-3

The manuscript application uses the third International Stroke Trial (IST-3)
dataset from Edinburgh DataShare:

- Record: <https://datashare.ed.ac.uk/handle/10283/1931>
- DOI: <https://doi.org/10.7488/ds/1350>
- Dataset citation: Sandercock, P; Wardlaw, J; Lindley, R; Cohen, G; Whiteley,
  W. (2016). The third International Stroke Trial (IST-3), 2000-2015 [dataset].
  University of Edinburgh & Edinburgh Clinical Trials Unit.

The embargo on the raw patient-level files expired on 25 January 2021. The
current DataShare record states that it is no longer necessary to apply to the
study investigators for access. Use nevertheless remains conditional on
acceptance of the IST-3 Data Use Agreement supplied with the record. This
includes its confidentiality, restricted-sharing, data-security, project-end
destruction, dataset-citation, and publication-acknowledgment requirements. The
patient-level files are ignored by Git and are not redistributed through this
repository. Do not describe these data as unrestricted open data.

## Raw files

The manuscript build can recreate the standardized analysis CSV from the raw
fixed-width file `ist3.dat`. If `manuscript/application-data/ist3.dat` is not
present, the build attempts to download it directly from Edinburgh DataShare:

```text
https://datashare.ed.ac.uk/bitstream/handle/10283/1931/ist3.dat?sequence=15&isAllowed=y
```

The fixed-width positions used in the build are taken from the IST-3 SAS syntax:

- `itt_treat`: column 27, allocated treatment.
- `dead6mo`: column 309, died before 6 months.
- `euroqol6`: columns 341-343, EuroQol health state at 6 months.
- `treatment`: columns 503-509, treatment text field used to verify coding.

The coding check requires `itt_treat == 0` to match `rt-PA` and `itt_treat == 1`
to match `Placebo`. The manuscript analysis recodes this as `R = 1` for rt-PA
and `R = 0` for placebo/control.

## Standardized file

Running the manuscript asset build creates:

```text
manuscript/application-data/ist3-standardized.csv
```

with the analysis columns:

- `Y`: `-1` for death by 6 months, otherwise 6-month EQ-VAS.
- `R`: binary randomized arm, `0` placebo/control and `1` rt-PA.
- `A`: indicator that `Y` is not the death atom.
- `atom`: `-1`.

Alive participants with missing 6-month EQ-VAS are excluded from the analysis
CSV. Deaths are retained through the `Y = -1` atom.

## Rebuilding manuscript assets

From the repository root:

```sh
Rscript manuscript/scripts/build-assets.R --target manuscript
```

The Bayesian analysis is cached at:

```text
manuscript/application-data/ist3-bayes-cache.rds
```

Ordinary builds validate and reuse this locked cache. They do not run a
Bayesian fit. Only when a new candidate fit has been explicitly requested,
generate it with:

```sh
TRUNCCOMP_REFRESH_IST3_BAYES=true Rscript manuscript/scripts/build-assets.R --target manuscript
```

The locked analysis uses the application-specific Stan source
`manuscript/stan/ist3_ordered_bounded_score_logit_normal_h2.stan`. It is an
ordered two-component bounded-score logit-normal model for EQ-VAS with
`score_min = 0`, `score_max = 100`, `score_step = 1`,
`heaping_grids = c(1, 5, 10)`, and shared reporting-grid weights. The fit used
four chains with 1,000 warmup and 2,000 sampling iterations per chain. This
application model is intentionally separate from the generic package
logit-normal mixtures under `inst/stan/`, whose component locations have
independent priors.

Default builds validate the exact settings, standardized data, embedded Stan
code, and locked cache fingerprints. A missing or mismatched cache stops the
build and never starts a fit. An explicitly requested candidate run uses the
application-specific ordered model and writes
`ist3-bayes-cache-candidate.rds`. It does not overwrite the locked publication
cache. Promoting a candidate requires an explicit review of its diagnostics and
posterior outputs followed by a deliberate update of the locked fingerprints.
Exploratory fits with different model or sampling settings must use a separate
script and output path rather than the publication cache path.
