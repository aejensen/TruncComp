# Application Data

Place the manually downloaded Dryad ARDS workbook here if you want the
manuscript build to use the clinical primary application instead of the
automatic RAND HIE fallback.

Supported local filenames:

- `dryad-ards-data.xlsx`
- `data.xlsx`
- `doi_10_5061_dryad_7d8k0__v20170929.zip`
- `dryad-ards-data.zip`
- `data.zip`
- `dryad-ards-standardized.csv`

The standardized CSV option should contain exactly two columns:

- `Y`: numeric outcome with atom coded as `0`
- `R`: binary group indicator coded `0/1`

If no supported local file is present, the manuscript build automatically uses
the public RAND Health Insurance Experiment fallback from Rdatasets.
