# Bayesian Small Area Estimation in R: Income and Poverty

R code comparing survey-weighted Gaussian, Binomial and joint spatial models for median income and poverty estimation across 88 Illinois PUMAs, using 2023 American Community Survey (ACS) data.

## Files

| File | Purpose |
|---|---|
| `main/packages.r` | Loads the required R packages. |
| `main/functions.r` | Defines `uni_gaus()` (Gaussian), `uni_bio()` (Binomial), `multi()` (joint), and shared transformation, prediction, Direct estimation and summary functions. |
| `main/empirical study.r` | Builds the fixed population, samples units, estimates lambda within each sample, fits three models and evaluates predictions against population truth. |
| `main/region study.r` | Processes the observed ACS sample and external population counts, fits three models, and produces regional estimates and maps. |
| `main/README.md` | Documents inputs, settings and execution. |
| `.gitignore` | Excludes generated results, local archives and manuscript files from Git. |

## Input data

Five source files are included in `main/data/`.

| File | Purpose |
|---|---|
| `IL23.csv` | 2023 Illinois ACS/PUMS person records for both studies. |
| `il_puma_2023_sf.rds` | PUMA geometry for the empirical spatial basis. |
| `il_puma_2020_cb_sf.rds` | 2020 cartographic PUMA geometry for the region study. |
| `acs2023_b15001_var_maps.rds` | B15001 variable mappings for sex and educational attainment. |
| `acs2023_il_puma_b15001_raw_cells.rds` | ACS estimates for the region study's population cells. |

The empirical population contains 92,817 positive-income records; cell counts come from those records. The region study uses 91,208 adults with positive income and external B15001 population counts totaling 9,848,062.

## Run the analyses

Install the packages listed in `main/packages.r`. Set the R working directory to the repository root, the parent of `main/`.

Run either study independently:

```r
source("./main/empirical study.r")
source("./main/region study.r")
```

Both scripts load their packages and functions automatically and can run independently.

Settings are defined at the top of each script:

| Study | Replicates | Sample size | Burn-in | Retained draws |
|---|---:|---:|---:|---:|
| Empirical | 100 | 1,000 per replicate | 1,000 | 1,000 |
| Region | One analysis | 91,208 | 10,000 | 30,000 |

Both use thinning 1. For a short execution check, use empirical `n_rep <- 2L` and `nburn <- 2L`, `nsim <- 4L` in either script.

Each run refits and predicts sequentially. The default region analysis is computationally intensive.

## Statistical targets

Both analyses predict all income units, including sampled units. Finite-population medians are computed on the original dollar scale after inverse transformation. Point estimates are posterior means of the median draws.

Box–Cox lambda is estimated from each empirical sample using a survey-weighted profile, and from the complete region sample using an unweighted profile.

## Outputs

| Location | Contents |
|---|---|
| `main/results/empirical/replicate_*.rds` | Per-replicate samples, fits, draws and estimates. |
| `main/results/empirical/empirical_results.rds` | Population truth, combined records and evaluation summaries. |
| `main/results/empirical/table2*.csv`, `table2.tex` | MSE, Ratio to Direct, interval score and coverage tables. |
| `main/results/region/model_fits.rds` | Posterior chains from all three models. |
| `main/results/region/area_draws.rds` | Regional prediction draws, prediction diagnostics, lambda and run settings. |
| `main/results/region/region_results.rds` | Regional estimates, variance comparisons and Direct diagnostics. |
| `main/results/region/tables/` | CSV tables of area estimates and model variance comparisons. |
| `main/results/region/figures/` | Maps and variance comparisons: `region_median_main.png` and `.pdf`. |

Generated results are ignored by Git. Rerunning a study overwrites its matching output files.

## Reproducibility

Empirical replicate `r` uses `set.seed(r)`; region uses seed `202309`. Fitting and prediction share the sequential random-number stream. Short-run repeatability has been checked; convergence and exact agreement with historical chains are separate validation tasks.
