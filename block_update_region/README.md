# Region Block-Update Diagnostic

This folder contains a sampler variant for the region study. It keeps the same
model, priors, and output structure as the current region workflow, but adds a
joint Metropolis-Hastings block move for `(tau_2, eta)` inside the multi-type
sampler.

Run the default diagnostic:

```bash
/opt/miniconda3/envs/R/bin/Rscript "block_update_region/region_study_block_update.r"
```

Default settings:

```r
REGION_DIAG_SEED = 20260706
REGION_DIAG_NBURN = 1000
REGION_DIAG_NSIM = 1000
REGION_DIAG_NTHIN = 1
REGION_DIAG_TAU2_INIT = -0.1
REGION_BLOCK_TAU_SD = 0.01
REGION_BLOCK_ETA_SD = 0.02
```

Outputs are ignored under `results/fig/region/`:

```text
region_block_update_basis_main.png
region_block_update_parameter_traceplots.png
```

Use `REGION_DIAG_LABEL` to keep multiple runs:

```bash
REGION_DIAG_LABEL=block_tau_neg1 \
REGION_DIAG_TAU2_INIT=-1 \
REGION_BLOCK_TAU_SD=0.01 \
REGION_BLOCK_ETA_SD=0.02 \
/opt/miniconda3/envs/R/bin/Rscript "block_update_region/region_study_block_update.r"
```
