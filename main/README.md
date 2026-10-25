# JSSAM code sample

R code for the empirical simulation (Table 2) and Illinois ACS application (Figure 2).

- `packages.r` loads the required packages.
- `functions.r` contains the model samplers `uni_gaus()`, `uni_bio()` and `multi()`, plus shared calculations.
- `empirical study.r` runs the simulation from the original data through evaluation.
- `region study.r` fits the regional models and produces tables and maps.

## Run

Install the packages listed in `packages.r`. Set the R working directory to the project root, the parent of `main/`, then run either study independently:

```r
source("./main/empirical study.r")
source("./main/region study.r")
```

Each script loads the packages and functions itself. Settings appear at the top of each script; paths are relative to the project root.

| Study | Replicates | Sample size | Burn-in | Retained draws |
|---|---:|---:|---:|---:|
| Empirical | 100 | 1,000 per replicate | 1,000 | 1,000 |
| Region | One analysis | 91,208 | 10,000 | 30,000 |

Both use thinning 1. For a quick check, use `n_rep <- 2L` in the empirical script and `nburn <- 2L`, `nsim <- 4L` in either script. These short runs are for checking execution only.

## Inputs

The five source files are included in `main/data/`; no downloads are required.

| File | Contents |
|---|---|
| `IL23.csv` | 2023 Illinois ACS/PUMS person records |
| `il_puma_2023_sf.rds` | PUMA geometry for the empirical study |
| `il_puma_2020_cb_sf.rds` | 2020 cartographic PUMA geometry for the region study |
| `acs2023_b15001_var_maps.rds` | B15001 sex and education variable mappings |
| `acs2023_il_puma_b15001_raw_cells.rds` | 2023 Illinois B15001 population cell estimates |

The empirical population consists of 92,817 positive-income records, with cell counts taken from those records. The region study uses 91,208 adults with positive income and B15001 population counts totaling 9,848,062. Both scripts show data processing and spatial basis construction explicitly.

## Analysis and outputs

Both studies fit Gaussian, Binomial and joint models. Income prediction covers all units, including sampled units. Finite-population median draws are computed on the original income scale and summarized by their posterior mean.

Results are written to `main/results/empirical/` and `main/results/region/`. Each run recomputes the analysis and overwrites files with matching names. Generated results are excluded from Git.

Each empirical replicate starts with `set.seed(r)`; the region study starts with seed `202309`. Fitting and prediction then consume random numbers in execution order, without separate prediction seeds. These settings reproduce new runs, while numerical results differ from earlier versions that reset seeds during prediction.

## Validation

Checks confirm unchanged sampler and prediction algorithms, matching single-draw predictions from the same initial random state, and reproducible complete short runs of both studies, including Table 2 and Figure 2 data.

These checks establish execution and reproducibility of short runs. Convergence and agreement of full-length runs with historical results have not been established by these checks.
