# TruncComp Monorepo

This repository hosts the `TruncComp` project material as a monorepo.

## Repository Layout

- `packages/TruncComp/`: the R package
- `manuscript/`: manuscript sources, simulations, and supporting project material

## Packages

### `TruncComp`

R package for two-sample comparison of truncated continuous outcomes with:

- a parametric likelihood-ratio test (`method = "LRT"`)
- a semi-parametric likelihood-ratio test (`method = "SPLRT"`)

Package README:
[packages/TruncComp/README.md](packages/TruncComp/README.md)

Install from GitHub:

```r
library(remotes)
install_github("aejensen/TruncComp", subdir = "packages/TruncComp")
```

## Development

The GitHub Actions workflow checks package directories from the monorepo layout.

For local package work:

```sh
cd packages/TruncComp
R CMD build .
R CMD check --no-manual --no-build-vignettes TruncComp_*.tar.gz
```
