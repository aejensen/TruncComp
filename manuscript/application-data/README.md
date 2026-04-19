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

The default Bayesian run uses two mixture components, four chains, and 300
warmup plus 300 sampling iterations per chain. Override these for a longer run,
for example:

```sh
TRUNCCOMP_REFRESH_IST3_BAYES=true \
TRUNCCOMP_IST3_BAYES_COMPONENTS=4 \
TRUNCCOMP_IST3_BAYES_WARMUP=500 \
TRUNCCOMP_IST3_BAYES_SAMPLING=500 \
Rscript manuscript/scripts/build-assets.R --target manuscript
```
