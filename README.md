# gasexchangeR

Shared R functions for adding/mapping species codes and fitting A-Ci and
temperature-response (Vcmax, Jmax, respiration) curves. Extracted for reuse
across multiple field study projects, it provides:

- `add_sp_code()`/`add_species()` for deriving/attaching 4-letter species
  codes,
- A-Ci curve fitting with either the [plantecophys](https://cran.r-project.org/package=plantecophys)
  or [PhotoGEA](https://cran.r-project.org/package=PhotoGEA) packages, and
- Arrhenius-type temperature-response fitting for the resulting
  photosynthetic parameters.

Reading/tidying raw LI-6800 exports and downstream sensitivity-analysis
functions are deliberately **not** part of this package — they stay local to
each project, since raw file formats/layouts and sensitivity-analysis needs
have so far diverged too much between projects to share safely.

## Installation

This package is not on CRAN. Install it directly from source with
[remotes](https://cran.r-project.org/package=remotes) or
[devtools](https://cran.r-project.org/package=devtools):

```r
# from a local clone
devtools::install("path/to/gasexchangeR")

# or straight from GitHub
remotes::install_github("CamilleSicangco/gasexchangeR")
```

A small worked example dataset, `aci_example`, is bundled with the package
(`data(aci_example)`; see `?aci_example`) — 3 individual *Eucalyptus
tereticornis* trees (`replicate`), each with a 5-temperature A-Ci curve
series (`curve_ID`, e.g. `"rep1_20C"` ... `"rep1_42C"`, named by replicate
and rounded leaf temperature) spanning ~20-42°C leaf temperature, subset
from Kumarathunge et al.'s published
[ACi-TGlob database](https://doi.org/10.6084/m9.figshare.7283567):

```r
data(aci_example)
head(aci_example)
#>   curve_ID replicate                 species        Ci          A    Tleaf
#> 1 rep1_20C         1 Eucalyptus tereticornis  49.39791  0.1614484 20.88185
#> 2 rep1_20C         1 Eucalyptus tereticornis  85.03205  3.5347838 21.32710
#> 3 rep1_20C         1 Eucalyptus tereticornis 121.85071  6.6689022 20.90660
#> 4 rep1_20C         1 Eucalyptus tereticornis 183.41322 10.9342214 19.41345
#> 5 rep1_20C         1 Eucalyptus tereticornis 270.21853 14.2373400 19.61558
#> 6 rep1_20C         1 Eucalyptus tereticornis 331.20626 15.7843464 19.66689
#>        Qin      Cond
#> 1 1799.048 0.7271726
#> 2 1799.366 0.7666477
#> 3 1798.858 0.7691114
#> 4 1800.062 0.7572762
#> 5 1799.365 0.7477825
#> 6 1799.925 0.7240140
```

It's used in the A-Ci and temperature-response examples below; see
`data-raw/aci_example.R` for how it was built.

## Workflow overview

### 1. Add/map species codes

```r
library(gasexchangeR)

df <- add_sp_code(df)                    # species -> 4-letter sp_code
df <- add_species(df, lookup = my_lookup) # sp_code -> species, via a study-specific lookup
```

Reading and tidying raw LI-6800 exports into a `data.frame` is left to each
project (raw file layouts and metadata-parsing conventions have diverged too
much between projects to share safely) — `add_sp_code()`/`add_species()`
pick up from there.

### 2. Fit A-Ci curves

```r
data(aci_example)

fits <- fit_ACis(aci_example, group_cols = "curve_ID", engine = "plantecophys")
r2 <- get_r2(fits, engine = "plantecophys")

curve_info <- unique(aci_example[, c("curve_ID", "replicate", "species", "Tleaf")])
curve_info$location <- "Richmond"
curve_info$method <- "whole-tree chamber"
params <- extract_curve_params(fits, curve_info, byname = "curve_ID", engine = "plantecophys")
```

`fit_ACis()` is a single, engine-aware entry point for fitting A-Ci curves
with either `PhotoGEA::fit_c3_aci()` (`engine = "PhotoGEA"`, taking a Licor
`exdf`) or `plantecophys::fitaci()` (`engine = "plantecophys"`, taking a
plain data frame — as above). `extract_curve_params()` builds a tidy data
frame of fitted Vcmax/Jmax-related quantities (including A at fixed Ci =
250/275/360 ppm and Rubisco-/RuBP-regeneration-limitation flags) merged
back onto the original curve-level data, for either engine.
`plot_c3_aci_fit_indiv()` and `plot_c3_aci_fit_indiv0()` give per-curve and
faceted diagnostic plots of PhotoGEA fits.

### 3. Fit temperature responses

```r
params$TsK <- params$Tleaf + 273.15

peak_fit <- fit_temperature_response(params, response = "Vcmax", return = "Peak")
arr_fit <- fit_arrhenius(params, response = "Jmax")
topt_fit <- fitA250vT(params, var = "A250")
```

`fit_temperature_response()` and `fit_arrhenius()` fit (peaked-)Arrhenius
models of Vcmax/Jmax vs. leaf temperature via `nlme`, with `replicate` as a
random effect. `fitA250vT()` fits a quadratic-Topt model of assimilation vs.
leaf temperature. Two respiration-vs-temperature fitters are provided, for
different models: `fitRvT_Ea()` fits an Arrhenius model (given a data frame
with `Rd` and `Tleaf` columns), and `fitRvT()` fits a quadratic-in-log mixed
model with a per-sample random effect (Heskel et al. 2016; given a data
frame with `R`, `sample`, and `Tleaf` columns). `arrhenius_formula()` and
`quadratic_topt_formula()` build the underlying model formulas directly if
you need a custom fit. `fit_regression()` is a general parameter-vs-MAT (or
similar) regression helper usable alongside these.

Climate/parameter-perturbation sensitivity analysis is deliberately **not**
part of this package — those functions stay local to each project.

## Documentation

Every exported function has a roxygen-generated help page (`man/`,
accessible via `?function_name` once the package is loaded/installed).

## Testing

Unit tests use [testthat](https://cran.r-project.org/package=testthat)
(edition 3):

```r
devtools::test()
```

## License

MIT — see [LICENSE.md](LICENSE.md).
