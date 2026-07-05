# unit-level-multi-type-model

This repository contains R code for fitting a unit-level multi-type model with one Gaussian response and one binary response under a spatial basis-function specification.

The current code is organized around two main tasks:

1. Region-level analysis using Illinois ACS PUMS data.
2. Empirical simulation based on repeated samples from the Illinois data.

The repository compares:
- a univariate Gaussian model,
- a univariate Bernoulli model,
- a joint multi-type model,
- and direct estimators used as benchmarks.

## Repository structure

- `packages.r`  
  Installs missing packages and loads all required libraries.

- `functions.r`  
  Contains the core functions used throughout the project, including:
  - `unis_gaus_fast()`
  - `unis_bios_fast()`
  - `MTSM_br_tau2_gaus_fast()`
  - `collapse_indicator_binomial()`
  - `gaus_post()`
  - `bios_post()`

- `region basis.r`  
  Runs the region-level basis-function analysis.

- `region basis results.r`  
  Loads the saved region-level fit and produces summary tables and figures.

- `empirical basis.r`  
  Runs the empirical simulation study.

- `empirical basis results.r`  
  Loads the saved simulation output and produces summary tables and figures.

- `HT_result.r`  
  Computes direct-estimator benchmark results for the simulation setting.

- `source data/IL23.csv`  
  Main input data file.

## Data and model setup

The code uses Illinois ACS PUMS data stored in `source data/IL23.csv`.

The responses are:
- Gaussian response: scaled log income
- Binary response: poverty indicator

The spatial effect is represented using basis functions built from the positive eigenvectors of the PUMA adjacency matrix.

The analysis is done at the PUMA level, with poststratification cells defined by:
- `PUMA`
- `SEX`
- `BACH`

## Requirements

The code is written in R.

All required packages are handled by:

```r
source("packages.r")
```

This script installs any missing packages and then loads them.

## Before running

Set the working directory to the root of this repository.

Make sure the following file is present:

```text
source data/IL23.csv
```

## Recommended workflow

Run the scripts from the repository root in the following order.

### 1. Load packages

```r
source("packages.r")
```

### 2. Load functions

```r
source("functions.r")
```

### 3. Run the region-level analysis

```r
source("region basis.r")
```

This script:
- reads the Illinois PUMS data,
- constructs poststratification cells,
- builds the spatial basis,
- fits the univariate Gaussian model,
- fits the univariate Bernoulli model,
- fits the joint multi-type model,
- saves the fitted objects to `data/region_basis_tau_in_gauss_results.RData`.

### 4. Summarize the region-level results

```r
source("region basis results.r")
```

This script:
- loads `data/region_basis_tau_in_gauss_results.RData`,
- computes region-level summaries,
- compares posterior variances across models,
- writes summary tables to `results/table/region/`,
- saves figures to `results/fig/region/`.

Typical outputs include:
- `results/table/region/basis_region_summary.csv`
- `results/table/region/basis_region_variance_summary.csv`
- `data/region_result_complete_basis.RData`
- figures in `results/fig/region/`

### 5. Run the empirical simulation study

```r
source("empirical basis.r")
```

This script:
- uses the Illinois data to define the area-level truth,
- builds the same spatial basis,
- repeatedly samples survey-like datasets,
- fits direct, univariate, and multi-type estimators for each replicate,
- stores the simulation output in `data/empirical_basis_tau_in_gauss_results.RData`.

Current default settings in the script are:
- `n_sim = 100`
- `sample_size = 1000`
- `nburn = 1000`
- `nsim = 1000`
- `nthin = 1`

### 6. Summarize the empirical simulation results

```r
source("empirical basis results.r")
```

This script:
- loads the saved simulation results,
- computes MSE, interval score, and coverage summaries,
- produces both replicate-level and area-level comparison tables,
- writes summary tables to `results/table/empirical/`,
- saves figures to `results/fig/empirical/`.

Typical outputs include:
- `results/table/empirical/basis_summary_table_area.csv`
- `results/table/empirical/basis_summary_table_rep.csv`
- `results/table/empirical/basis_absolute_table.csv`
- `results/table/empirical/basis_paired_table.csv`
- `results/table/empirical/basis_better_prop_table.csv`
- `results/table/empirical/basis_paired_replicate_metrics.csv`
- `results/table/empirical/basis_paired_area_metrics.csv`
- `results/table/empirical/basis_better_prop_area_table.csv`
- figures in `results/fig/empirical/`

### 7. Optional: direct-estimator benchmark

```r
source("HT_result.r")
```

This script computes direct-estimator benchmark metrics for the same simulation design.

## Output folders

The main scripts create these folders if they do not already exist:

```text
data/
results/
  fig/
    empirical/
    region/
  table/
    empirical/
    region/
```

Typical outputs are:
- saved `.RData` objects in `data/`
- `.csv` summary tables in `results/table/`
- `.png` figures in `results/fig/`

## Notes

- `region basis.r` will use `IL_poststrat_cells.csv` if it is already available.  
  Otherwise, it constructs the poststratification cells directly from the sample data.

- The basis-function specification keeps the eigenvectors associated with positive eigenvalues of the adjacency matrix.

- The same modeling functions in `functions.r` are used by both the region analysis and the empirical simulation.

## Minimal run example

If you want the main analysis pipeline only, run:

```r
source("packages.r")
source("functions.r")
source("region basis.r")
source("region basis results.r")
source("empirical basis.r")
source("empirical basis results.r")
```

## Reproducibility

For a clean run:

1. Clone the repository.
2. Set the working directory to the repository root.
3. Make sure `source data/IL23.csv` is available.
4. Run the scripts in the order listed above.

## Contact

This repository contains the code used for the unit-level multi-type modeling project.
