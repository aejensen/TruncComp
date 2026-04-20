# Application Data: IST-3

The manuscript application uses the third International Stroke Trial (IST-3)
dataset from Edinburgh DataShare:

- Record: <https://datashare.ed.ac.uk/handle/10283/1931>
- DOI: <https://doi.org/10.7488/ds/1350>
- Dataset citation: Sandercock, P; Wardlaw, J; Lindley, R; Cohen, G; Whiteley,
  W. (2016). The third International Stroke Trial (IST-3), 2000-2015 [dataset].
  University of Edinburgh & Edinburgh Clinical Trials Unit.

The data are publicly downloadable after the embargo period, but reuse is
conditional on the IST-3 Data Use Agreement supplied with the DataShare record.
Do not describe these data as unrestricted open data.

## Raw files

The manuscript build can recreate the standardized analysis CSV from the raw
fixed-width file `ist3.dat`. If `manuscript/application-data/ist3.dat` is not
present, the build looks for the development copy at:

```text
/tmp/trunccomp-dataset-hunt/ist3/ist3.dat
```

If neither file is present, the build attempts to download `ist3.dat` directly
from Edinburgh DataShare:

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

## Rebuilding

From the repository root:

```sh
Rscript manuscript/scripts/build-assets.R --target manuscript
```

The Bayesian analysis is cached at:

```text
manuscript/application-data/ist3-bayes-cache.rds
```

Refresh the Bayesian cache with:

```sh
TRUNCCOMP_REFRESH_IST3_BAYES=true Rscript manuscript/scripts/build-assets.R --target manuscript
```

The default IST-3 Bayesian run uses the package bounded-score logit-normal
model for EQ-VAS: `continuous_support = "bounded_score"`,
`bounded_kernel = "logit_normal"`, `score_min = 0`, `score_max = 100`,
`score_step = 1`, `heaping_grids = c(1, 5, 10)`, and shared reporting-grid
weights. It uses two mixture components, four chains, and 1,000 warmup plus
2,000 sampling iterations per chain. Override the sampling length or truncation
level for an exploratory run, for example:

```sh
TRUNCCOMP_REFRESH_IST3_BAYES=true \
TRUNCCOMP_IST3_BAYES_COMPONENTS=4 \
TRUNCCOMP_IST3_BAYES_WARMUP=1000 \
TRUNCCOMP_IST3_BAYES_SAMPLING=2000 \
Rscript manuscript/scripts/build-assets.R --target manuscript
```

# Appendix Application Data: joineR liver

The appendix application uses the `liver` data distributed with the CRAN
package `joineR` version 1.2.8:

- Data-loading call: `data(liver, package = "joineR")`
- Package citation: Philipson, P; Sousa, I; Diggle, P; Williamson, P;
  Kolamunnage-Dona, R; Henderson, R; Hickey, G. (2018). `joineR`: Joint
  Modelling of Repeated Measurements and Time-to-Event Data.
- Data description: liver cirrhosis controlled trial of prednisone, with
  repeated prothrombin index measurements and survival/censoring information.

The joineR documentation states that `treatment` is coded `0 = placebo` and
`1 = prednisone`, and that `cens` is coded `1 = died` and `0 = censored`.

## Appendix endpoint

The appendix endpoint is a two-year landmark atom-plus-continuous endpoint:

- `Y = 0` for death before or at two years.
- `Y =` latest prothrombin index (%) at or before two years for subjects known
  alive at two years.
- `R = 0` for placebo and `R = 1` for prednisone.
- `A = 1(Y != 0)`.
- `atom = 0`.

Subjects censored before two years are excluded because their two-year death
atom status is not observed. In the packaged data this excludes 65 of 488
participants. All subjects known alive at two years have at least one
prothrombin measurement at or before the landmark.

Running the manuscript asset build creates:

```text
manuscript/application-data/liver-appendix-standardized.csv
```

The Bayesian appendix fit is cached at:

```text
manuscript/application-data/liver-appendix-bayes-cache.rds
```

Refresh the liver Bayesian cache with:

```sh
TRUNCCOMP_REFRESH_LIVER_BAYES=true Rscript manuscript/scripts/build-assets.R --target manuscript
```

The default liver Bayesian run uses four mixture components, positive-real
continuous support, four chains, and 300 warmup plus 300 sampling iterations per
chain. Override these for a longer run, for example:

```sh
TRUNCCOMP_REFRESH_LIVER_BAYES=true \
TRUNCCOMP_LIVER_BAYES_COMPONENTS=4 \
TRUNCCOMP_LIVER_BAYES_WARMUP=500 \
TRUNCCOMP_LIVER_BAYES_SAMPLING=500 \
Rscript manuscript/scripts/build-assets.R --target manuscript
```

# Appendix Application Data: licorice gargle randomized trial

The licorice gargle appendix application uses the `licorice_gargle` data
distributed with the CRAN package `medicaldata` version 0.2.0:

- Data-loading call: `data(licorice_gargle, package = "medicaldata")`
- Package citation: Higgins, P. (2021). `medicaldata`: Data Package for
  Medical Datasets.
- Data-source citation: Nowacki, A. S. (2017). Licorice Gargle Dataset,
  Teaching of Statistics in the Health Sciences Resources Portal.
- Source trial: Ruetzler et al. (2013). A randomized, double-blind comparison
  of licorice versus sugar-water gargle for prevention of postoperative sore
  throat and postextubation coughing.

The `medicaldata` documentation states that `treat` is coded `0 = Sugar 5g`
and `1 = Licorice 0.5g`. The source paper and package documentation describe
the pain scores as 11-point Likert scores, with 0 meaning no pain and 10 the
worst pain.

The CRAN package is publicly reproducible under the package license. The
original TSHS portal is an educational resource and states that publication
reuse permission is not granted by the portal itself, so the appendix records
that caveat explicitly.

## Licorice endpoint

The appendix endpoint uses the 30-minute PACU pain-on-swallowing score:

- `Y = 0` for no sore-throat pain during swallowing at 30 minutes.
- `Y =` the positive 30-minute pain-on-swallowing score for participants with
  symptoms.
- `R = 0` for sugar-water gargle and `R = 1` for licorice gargle.
- `A = 1(Y != 0)`.
- `atom = 0`.

Because larger pain scores are worse and the no-pain atom is favorable, clinical
benefit corresponds to lower `Delta`, lower `mu_delta`, an odds ratio
`alpha_delta < 1` for any pain versus no pain, and a higher no-pain atom
probability. Two packaged records have missing 30-minute pain-on-swallowing
scores, one in each randomized group, and are excluded from this appendix
analysis.

Running the manuscript asset build creates:

```text
manuscript/application-data/licorice-appendix-standardized.csv
```

The Bayesian appendix fit is cached at:

```text
manuscript/application-data/licorice-appendix-bayes-cache.rds
```

Refresh the licorice Bayesian cache with:

```sh
TRUNCCOMP_REFRESH_LICORICE_BAYES=true Rscript manuscript/scripts/build-assets.R --target manuscript
```

The default licorice Bayesian run uses two mixture components, positive-real
continuous support, four chains, and 300 warmup plus 300 sampling iterations per
chain. Override these for a longer run, for example:

```sh
TRUNCCOMP_REFRESH_LICORICE_BAYES=true \
TRUNCCOMP_LICORICE_BAYES_COMPONENTS=2 \
TRUNCCOMP_LICORICE_BAYES_WARMUP=500 \
TRUNCCOMP_LICORICE_BAYES_SAMPLING=500 \
Rscript manuscript/scripts/build-assets.R --target manuscript
```
