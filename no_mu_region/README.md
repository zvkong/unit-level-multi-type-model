# Region No-Mu Diagnostic

This folder contains a non-block diagnostic version of the region study. It keeps the
same model, priors, sampler, predictions, poststratified summaries, parameter trace
plots, and MCMC diagnostic tables as the current non-block diagnostic script, but it
does not store `Mu.chain`, `Mu_1.chain`, or `Mu_2.chain`.

This is intended for longer parameter-level runs without the large memory cost of
storing individual-level fitted means.

Run the default diagnostic:

```bash
/opt/miniconda3/envs/R/bin/Rscript "no_mu_region/region_study_no_mu_diagnostic.r"
```

Default settings:

```r
REGION_DIAG_SEED = 20260706
REGION_DIAG_NBURN = 5000
REGION_DIAG_NSIM = 10000
REGION_DIAG_NTHIN = 1
REGION_DIAG_TAU2_INIT = 0
```

Default outputs:

```text
data/region_no_mu_tau0_burn5000_results.RData
results/fig/region/region_no_mu_tau0_burn5000_basis_main.png
results/fig/region/region_no_mu_tau0_burn5000_parameter_traceplots.png
results/table/region/region_no_mu_tau0_burn5000_mcmc_key_diagnostics.csv
results/table/region/region_no_mu_tau0_burn5000_mcmc_parameter_diagnostics.csv
results/table/region/region_no_mu_tau0_burn5000_region_summary.csv
```

This version intentionally does not produce `*_mcmc_mu_diagnostics.csv`, because
the Mu chains are not saved.
