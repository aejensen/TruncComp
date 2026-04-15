# Developing TruncComp2

`TruncComp2` is intended to evolve independently inside its own package
subtree.

Ordinary package work should only require changes under:

- `R/`
- `man/`
- `tests/`
- `tools/`
- package-local docs in this directory

No edits outside `packages/TruncComp2/` should be required for routine
development, testing, or release preparation.

## Package-Local Verification

From the package root, run:

```sh
Rscript tools/check-package.R
```

That script:

- installs helper packages needed for local verification if they are missing
- installs package dependencies for `TruncComp2`
- installs `TruncComp2` into a temporary library
- runs the package-local `testthat` suite
- runs `R CMD build`
- runs `R CMD check --no-manual --no-build-vignettes`

## Manual Commands

If you prefer to run the steps yourself from the package root:

```sh
R CMD build .
R CMD check --no-manual --no-build-vignettes TruncComp2_*.tar.gz
```

For the test suite after installation:

```sh
Rscript -e 'library(testthat); library(TruncComp2); testthat::test_dir("tests/testthat", reporter = "summary")'
```

## Fixture Maintenance

Frozen reference fixtures can be regenerated with:

```sh
Rscript tools/generate-el-fixture.R
```

This script is package-local and writes only into `tests/testthat/fixtures/`.

For the adjusted semi-parametric extension, a package-local validation script is
also available:

```sh
Rscript tools/validate-adjusted-splrt.R
```

That script runs fixed-seed mini-simulation checks for approximate null
calibration and basic power. It is intended for methodological validation, not
for routine CI.
