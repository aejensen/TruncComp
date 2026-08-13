# Manuscript

This directory contains the source, generated assets, and submission variants
for *Two-component treatment comparisons for continuous outcomes truncated by a
substantive event*.

## Sources of truth

- `manuscript.tex` is the canonical manuscript source. Edit this file rather
  than any generated source under `self_contained/`.
- Every manuscript content change must be made first in this canonical
  reference TeX source. Run `make -C manuscript self-contained` afterward so
  the generated self-contained sources and PDFs match it.
- `bibliography.bib` is the manuscript bibliography.
- `Makefile` provides the supported build entry points.

## Directory structure

| Path | Purpose |
| --- | --- |
| `R/` | R functions for loading data, computing analyses, and generating manuscript values, tables, and figures. |
| `scripts/build-assets.R` | Entry point for generating the assets used by the LaTeX sources. |
| `scripts/build-self-contained.sh` | Creates and verifies the four flattened manuscript variants. |
| `build/generated-values.tex` | Generated LaTeX macros containing numerical results. |
| `build/tables/` | Generated LaTeX tables used by the manuscript. |
| `build/figures/` | Generated PDF figures used by the LaTeX sources. |
| `application-data/` | Application data, standardized analysis files, data documentation, and the locked IST-3 Bayesian cache. |
| `stan/` | Stan source for the locked IST-3 application model. |
| `simulation-study-results.rds` | Simulation results used to generate the manuscript summaries and figures. |
| `self_contained/` | Generated flattened TeX sources, compiled PDFs, the required proof style file, and copied figure dependencies. |
| `tests/` | Local checks for the locked application fit contract. |

Files under `build/` and `self_contained/` are generated. They should not be
edited by hand.

The asset-generation code loads the active `TruncComp` package directly from
the repository root. No nested package checkout is required. Set
`TRUNCCOMP_USE_LOCAL_PACKAGE=false` only when the manuscript should use an
installed `TruncComp` package instead.

## Building the manuscript

Run the supported targets from the repository root.

```sh
make -C manuscript pdf
make -C manuscript self-contained
make -C manuscript submission-pdf
```

`pdf` generates the manuscript assets and compiles `manuscript.pdf`.
`self-contained` uses the existing manuscript assets to produce and verify all
four flattened submission variants. `make -C manuscript all` builds the
ordinary manuscript PDF.

`submission-pdf` is the cache-free build path for a fresh clone. It compiles the
tracked appendix-proofs source and figure dependencies directly under
`self_contained/`. It does not generate assets, access the IST-3 data or cache,
or start a Bayesian fit. The full `pdf` target remains the validation path for
locally available application assets and stops if the locked cache is missing
or invalid.

The ordinary asset build validates and reuses the locked IST-3 Bayesian cache.
It does not refit that model.

Every PDF build removes LaTeX intermediate files. The cleanup can also be run
directly with

```sh
make -C manuscript clean-latex-temp
```

## Self-contained variants

The `self_contained/` directory contains the following generated TeX and PDF
pairs.

| Stem | Contents |
| --- | --- |
| `two_component_treatment_comparisons_inline_proofs` | Complete manuscript with each proof placed inline. |
| `two_component_treatment_comparisons_appendix_proofs` | Complete manuscript with proofs collected in the supplementary appendix. This is the preferred single-PDF submission version. |
| `two_component_treatment_comparisons_main` | Main manuscript for use with a separate supplementary document. |
| `two_component_treatment_comparisons_supplementary` | External proof stream and the remaining supplementary material. |

The flattened TeX files contain the generated values, tables, bibliography,
and, where applicable, the proof stream. Figures remain external dependencies
under `self_contained/figures/`. The build verifies both complete variants and
the separate main and supplementary pair before replacing the existing bundle.

## Statistics in Medicine submission

For a new submission, upload the following single file with the designation
**Main Document**.

```text
self_contained/two_component_treatment_comparisons_appendix_proofs.pdf
```

This combined PDF contains the complete article, figures, tables, references,
proofs, and supplementary material needed for peer review. Do not upload the
TeX sources or the separate main and supplementary PDFs at the initial stage
unless the submission system or editorial office specifically requests them.

Enter the following information separately in the submission system.

- Full title, author names, affiliations, and ORCID identifiers
- Short title `Two-component comparisons for truncated outcomes`
- Abstract and keywords from the manuscript
- Complete corresponding-author contact information for Andreas Kryger Jensen
- Funding information

Before creating the submission PDF, confirm that the manuscript contains the
final `Funding` and `Data and code availability` statements and that the code
availability statement contains the permanent repository URL or DOI.

For a revised submission, upload the rebuilt combined PDF as **Main Document**.
Also upload the following files with the designation **Supplemental Material
not for review**.

- `self_contained/two_component_treatment_comparisons_appendix_proofs.tex`
- `self_contained/proof-at-the-end.sty`
- Every referenced figure dependency under `self_contained/figures/`

Use the appendix-proofs source unless the editorial office asks for separate
main and supplementary documents. Delete superseded uploads in the submission
system.

## Locked IST-3 application fit

The publication fit at `application-data/ist3-bayes-cache.rds` is locked. It
uses the bounded-score logit-normal model with two ordered components, shared
reporting-grid weights for widths 1, 5, and 10, four chains, 1,000 warmup
iterations, and 2,000 sampling iterations per chain.

Ordinary manuscript builds validate and reuse this cache. They do not refit the
model. A missing or incompatible cache causes the build to stop. Do not rerun,
refresh, replace, or promote an alternative fit unless this has been explicitly
requested. See `application-data/README.md` for the data provenance and model
contract.
