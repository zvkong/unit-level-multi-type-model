# Unit-Level Multi-Type Model

This repository contains R code for a Bayesian unit-level multi-type model that jointly analyzes a Gaussian outcome and a Binomial outcome with shared latent structure.

The project is motivated by small area estimation problems where related outcomes can borrow strength from each other and across domains. The repository includes model functions, an empirical example based on Illinois ACS/PUMS data, and a simulation study for comparing the proposed method with univariate alternatives and direct estimators.

## Repository Structure

- `packages.r`  
  Installs and loads the R packages required for the project.

- `functions.r`  
  Contains the main modeling and evaluation functions, including:
  - univariate Gaussian model
  - univariate Binomial model
  - multi-type model with Gaussian-specific random effects
  - multi-type model with Binomial-specific random effects
  - poststratification functions
  - interval score and utility functions

- `region.r`  
  Runs the empirical study using Illinois ACS/PUMS data and produces PUMA-level estimates, posterior variance comparisons, and spatial plots.

- `empirical simulation.r`  
  Runs the simulation study and compares direct estimators, univariate models, and multi-type models using metrics such as MSE and interval score.

- `README.md`  
  Project description and usage guide.

## Model Overview

This repository implements Bayesian unit-level models for two response types:

- a Gaussian response for transformed income
- a Binomial response for poverty status

The two responses are linked through a shared latent structure so that information can be borrowed across outcomes. The Binomial component is handled using Pólya–Gamma data augmentation for efficient Gibbs sampling.

The code includes:

- a univariate Gaussian model
- a univariate Binomial model
- a multi-type model with extra structure favoring the Gaussian block
- a multi-type model with extra structure favoring the Binomial block

After model fitting, posterior predictions are poststratified to the PUMA level.

## Data

The empirical script uses 2021 Illinois ACS Public Use Microdata Sample (PUMS) variables such as:

- `PINCP`
- `PWGTP`
- `PUMA`
- `SEX`
- `SCHL`
- `POVPIP`

In the empirical example:

- income is transformed using log income and then scaled to the unit interval
- poverty is defined as an indicator based on `POVPIP <= 100`
- education is simplified into a bachelor’s degree indicator

The simulation script currently expects a local input file:

- `IL21.csv`

This file is not included in the repository, so the simulation will not run out of the box unless you prepare the required data first.

## Requirements

The code is written in R.

Run the following command first:

```r
source("packages.r")
````

The main dependencies include:

* `Matrix`
* `MASS`
* `mvtnorm`
* `invgamma`
* `BayesLogit`
* `coda`
* `Metrics`
* `sf`
* `spdep`
* `spatialreg`
* `ggplot2`
* `RColorBrewer`
* `scales`
* `tidyverse`
* `dplyr`
* `tidycensus`
* `tigris`
* `sampling`
* `survey`
* `HDInterval`
* `readr`

## Typical Workflow

### 1. Load required packages

```r
source("packages.r")
```

### 2. Load model and utility functions

```r
source("functions.r")
```

### 3. Run the empirical analysis

```r
source("region.r")
```

### 4. Run the simulation study

```r
source("empirical simulation.r")
```

## Outputs

Depending on the script, the code can produce:

* posterior mean estimates
* posterior variance summaries
* PUMA-level poststratified estimates
* spatial comparison maps
* MSE comparisons
* interval score comparisons
* variance ratio plots
* trace plots for selected model parameters

## Notes

* The repository is organized as a research codebase rather than a packaged R library.
* Some inputs are created inside the scripts, while others must be prepared separately.
* The empirical example focuses on Illinois PUMAs.
* The simulation study compares multi-type and univariate models against direct estimators.

## Author

Zewei Kong
