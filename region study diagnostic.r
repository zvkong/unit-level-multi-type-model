## Region study diagnostic run.
##
## This script reruns the same region model with diagnostic-friendly MCMC
## settings. It does not change the model, priors, or sampler. By default it
## reuses the already saved region input objects so that ACS/tigris data do not
## need to be fetched again.

required_pkgs <- c(
  "dplyr", "tidyr", "tibble", "ggplot2", "RColorBrewer", "patchwork", "sf", "stringr",
  "coda"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
  library(sf)
  library(stringr)
  library(coda)
})

source("functions.r")

diag_seed <- as.integer(Sys.getenv("REGION_DIAG_SEED", unset = "20260706"))
nburn <- as.integer(Sys.getenv("REGION_DIAG_NBURN", unset = "3000"))
nsim <- as.integer(Sys.getenv("REGION_DIAG_NSIM", unset = "5000"))
nthin <- as.integer(Sys.getenv("REGION_DIAG_NTHIN", unset = "5"))
tau_2_init <- as.numeric(Sys.getenv("REGION_DIAG_TAU2_INIT", unset = "-0.1"))
diag_label <- Sys.getenv("REGION_DIAG_LABEL", unset = "region_diagnostic")

if (nsim %% nthin != 0) {
  stop("REGION_DIAG_NSIM must be divisible by REGION_DIAG_NTHIN.", call. = FALSE)
}

prep_file <- Sys.getenv(
  "REGION_DIAG_PREP_FILE",
  unset = file.path("data", "region_basis_tau_in_gauss_results.RData")
)

region_fig_dir <- file.path("results", "fig", "region")
region_table_dir <- file.path("results", "table", "region")
basis_main_file <- file.path(region_fig_dir, paste0(diag_label, "_basis_main.png"))
parameter_trace_file <- file.path(region_fig_dir, paste0(diag_label, "_parameter_traceplots.png"))
key_diagnostics_file <- file.path(region_table_dir, paste0(diag_label, "_mcmc_key_diagnostics.csv"))
parameter_diagnostics_file <- file.path(region_table_dir, paste0(diag_label, "_mcmc_parameter_diagnostics.csv"))
mu_diagnostics_file <- file.path(region_table_dir, paste0(diag_label, "_mcmc_mu_diagnostics.csv"))
region_summary_file <- file.path(region_table_dir, paste0(diag_label, "_region_summary.csv"))

dir.create(region_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(region_table_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading region inputs from: ", prep_file)
prep_env <- new.env(parent = emptyenv())
load(prep_file, envir = prep_env)
required_inputs <- c("pums", "pcells", "basis_mat", "puma_sf")
missing_inputs <- setdiff(required_inputs, ls(prep_env))
if (length(missing_inputs) > 0) {
  stop("Missing required input object(s): ", paste(missing_inputs, collapse = ", "), call. = FALSE)
}

pums <- get("pums", envir = prep_env)
pcells <- get("pcells", envir = prep_env)
basis_mat <- get("basis_mat", envir = prep_env)
puma_sf <- get("puma_sf", envir = prep_env)
if (exists("area_levels", envir = prep_env)) {
  area_levels <- get("area_levels", envir = prep_env)
}

## Drop the prep environment so saved model chains are not kept in memory.
rm(prep_env)
gc()

if (!exists("area_levels")) {
  area_levels <- rownames(basis_mat)
}

set.seed(diag_seed)

predX <- stats::model.matrix(~ SEX + BACH - 1, data = pcells)
predS <- basis_mat[match(pcells$PUMA, rownames(basis_mat)), , drop = FALSE]

modX <- stats::model.matrix(~ SEX + BACH - 1, data = pums)
modY <- pums$INCO
modZ <- pums$POV
modS <- basis_mat[match(pums$PUMA, rownames(basis_mat)), , drop = FALSE]

message(
  "Running diagnostic MCMC with seed=", diag_seed,
  ", nburn=", nburn,
  ", nsim=", nsim,
  ", nthin=", nthin,
  ", tau_2_init=", tau_2_init
)

unis_wage <- unis_gaus_fast(
  X = modX,
  Y = modY,
  S = modS,
  sig2b = 1000,
  wgt = pums$scaledWGT,
  predX = predX,
  predS = predS,
  nburn = nburn,
  nsim = nsim,
  nthin = nthin,
  a = 0.5,
  b = 0.5,
  a_eps = 0.1,
  b_eps = 0.1
)

unis_pov <- unis_bios_fast(
  X = modX,
  Y = modZ,
  S = modS,
  sig2b = 1000,
  wgt = pums$scaledWGT,
  n = NULL,
  predX = predX,
  predS = predS,
  nburn = nburn,
  nsim = nsim,
  nthin = nthin,
  a = 0.1,
  b = 0.1
)

mult_br <- MTSM_br_tau2_gaus_fast(
  X_1 = modX,
  X_2 = modX,
  Z_1 = modY,
  Z_2 = modZ,
  S = modS,
  sig2b = 1000,
  wgt = pums$scaledWGT,
  n = NULL,
  predX = predX,
  predS = predS,
  n_preds = NULL,
  nburn = nburn,
  nsim = nsim,
  nthin = nthin,
  sig2t = 10,
  sig2e = 10,
  tau_2_init = tau_2_init,
  a_eps = 0.1,
  b_eps = 0.1,
  aeta = 0.1,
  beta = 0.1,
  alambda = 2,
  blambda = 1
)

fit_ug <- unis_wage
fit_ub <- unis_pov
fit_mt <- mult_br

true_mean_dummy <- stats::setNames(rep(NA_real_, length(area_levels)), area_levels)

## Make the poststratified simulation step reproducible separately from MCMC.
set.seed(diag_seed + 100000L)
res_ug <- gaus_post(fit_ug$Preds, fit_ug$sig2.chain, true_mean_dummy, pcells$PUMA, pcells$popsize)
res_mg <- gaus_post(fit_mt$preds_gaus.chain, fit_mt$sig2.chain, true_mean_dummy, pcells$PUMA, pcells$popsize)
res_ub <- bios_post(fit_ub$Preds, true_mean_dummy, pcells$PUMA, pcells$popsize)
res_mb <- bios_post(fit_mt$preds_bios.chain, true_mean_dummy, pcells$PUMA, pcells$popsize)

puma_stats <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(
    n_sample = dplyr::n(),
    ht_inco = stats::weighted.mean(INCO, PWGTP),
    ht_pov = stats::weighted.mean(POV, PWGTP),
    .groups = "drop"
  ) |>
  dplyr::arrange(PUMA)

region_summary <- puma_stats |>
  dplyr::left_join(
    tibble::tibble(
      PUMA = area_levels,
      multi_gaussian_mean = res_mg$est,
      uni_gaussian_mean = res_ug$est,
      multi_gaussian_var = as.numeric(res_mg$sigma2),
      uni_gaussian_var = as.numeric(res_ug$sigma2),
      multi_binomial_mean = res_mb$est,
      uni_binomial_mean = res_ub$est,
      multi_binomial_var = as.numeric(res_mb$sigma2),
      uni_binomial_var = as.numeric(res_ub$sigma2)
    ),
    by = "PUMA"
  ) |>
  dplyr::mutate(
    gaussian_uni_over_multi_var_ratio = uni_gaussian_var / multi_gaussian_var,
    binomial_uni_over_multi_var_ratio = uni_binomial_var / multi_binomial_var
  )

chain_diagnostics <- function(x, parameter, group) {
  x <- as.numeric(x)
  n <- length(x)
  first <- x[seq_len(floor(n * 0.1))]
  last <- x[(n - floor(n * 0.1) + 1):n]
  ac <- stats::acf(x, lag.max = 10, plot = FALSE, na.action = na.pass)$acf
  ess <- suppressWarnings(as.numeric(coda::effectiveSize(coda::mcmc(x))))
  geweke_z <- suppressWarnings(as.numeric(coda::geweke.diag(coda::mcmc(x))$z))

  tibble::tibble(
    group = group,
    parameter = parameter,
    n_saved = n,
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    q025 = as.numeric(stats::quantile(x, 0.025, na.rm = TRUE)),
    q975 = as.numeric(stats::quantile(x, 0.975, na.rm = TRUE)),
    lag1_acf = as.numeric(ac[2]),
    lag10_acf = as.numeric(ac[11]),
    ess = ess,
    ess_ratio = ess / n,
    geweke_z = geweke_z,
    first_last_diff_sd = (mean(last, na.rm = TRUE) - mean(first, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE),
    flag = abs(geweke_z) > 2 | ess < 100 | as.numeric(ac[2]) > 0.9
  )
}

matrix_diagnostics <- function(chain_mat, group, prefix) {
  chain_mat <- as.matrix(chain_mat)
  parameter_names <- rownames(chain_mat)
  if (is.null(parameter_names)) {
    parameter_names <- paste0(prefix, "[", seq_len(nrow(chain_mat)), "]")
  }

  dplyr::bind_rows(lapply(seq_len(nrow(chain_mat)), function(i) {
    chain_diagnostics(chain_mat[i, ], parameter_names[i], group)
  }))
}

selected_mu_indices <- function(mu_mean, n_select = 5) {
  probs <- seq(0.1, 0.9, length.out = n_select)
  targets <- as.numeric(stats::quantile(mu_mean, probs = probs, na.rm = TRUE))
  unique(vapply(targets, function(target) {
    which.min(abs(mu_mean - target))
  }, integer(1)))
}

mu_selected_diagnostics <- function(chain_mat, selected_idx, group) {
  dplyr::bind_rows(lapply(seq_along(selected_idx), function(i) {
    chain_diagnostics(chain_mat[selected_idx[i], ], paste0("mu_q", i, "_row_", selected_idx[i]), group)
  }))
}

mu_difference_diagnostics <- function(multi_chain, uni_chain, selected_idx, group) {
  dplyr::bind_rows(lapply(seq_along(selected_idx), function(i) {
    chain_diagnostics(
      multi_chain[selected_idx[i], ] - uni_chain[selected_idx[i], ],
      paste0("diff_q", i, "_row_", selected_idx[i]),
      group
    )
  }))
}

parameter_diagnostics <- dplyr::bind_rows(
  matrix_diagnostics(fit_ug$Beta.chain, "Univariate Gaussian beta", "beta"),
  matrix_diagnostics(fit_ub$Beta.chain, "Univariate Bernoulli beta", "beta"),
  matrix_diagnostics(fit_mt$Beta_1.chain, "Multi-type Gaussian beta_1", "beta_1"),
  matrix_diagnostics(fit_mt$Beta_2.chain, "Multi-type Bernoulli beta_2", "beta_2"),
  matrix_diagnostics(fit_ug$U.chain, "Univariate Gaussian spatial U", "U"),
  matrix_diagnostics(fit_ub$U.chain, "Univariate Bernoulli spatial U", "U"),
  matrix_diagnostics(fit_mt$eta.chain, "Multi-type shared eta", "eta"),
  matrix_diagnostics(fit_mt$lambda.chain, "Multi-type residual lambda", "lambda"),
  chain_diagnostics(fit_ug$sig2.chain, "sigma2", "Univariate Gaussian variance"),
  chain_diagnostics(fit_ug$Sigma2_u.chain, "Sigma2_u", "Univariate Gaussian variance"),
  chain_diagnostics(fit_ub$Sigma2_u.chain, "Sigma2_u", "Univariate Bernoulli variance"),
  chain_diagnostics(fit_mt$sig2.chain, "sigma2", "Multi-type Gaussian variance"),
  chain_diagnostics(fit_mt$Sigma2_eta.chain, "Sigma2_eta", "Multi-type shared variance"),
  chain_diagnostics(fit_mt$Sigma2_lambda.chain, "Sigma2_lambda", "Multi-type residual variance"),
  chain_diagnostics(fit_mt$tau_2.chain, "tau_2", "Multi-type shared loading")
) |>
  dplyr::mutate(
    diag_label = diag_label,
    seed = diag_seed,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    tau_2_init = tau_2_init,
    .before = 1
  )

eta_lambda_summary <- parameter_diagnostics |>
  dplyr::filter(group %in% c("Multi-type shared eta", "Multi-type residual lambda")) |>
  dplyr::group_by(group) |>
  dplyr::summarise(
    diag_label = dplyr::first(diag_label),
    seed = dplyr::first(seed),
    nburn = dplyr::first(nburn),
    nsim = dplyr::first(nsim),
    nthin = dplyr::first(nthin),
    tau_2_init = dplyr::first(tau_2_init),
    parameter = "component_summary",
    n_saved = dplyr::first(n_saved),
    mean = NA_real_,
    sd = NA_real_,
    q025 = NA_real_,
    q975 = NA_real_,
    lag1_acf = max(lag1_acf, na.rm = TRUE),
    lag10_acf = max(lag10_acf, na.rm = TRUE),
    ess = min(ess, na.rm = TRUE),
    ess_ratio = min(ess_ratio, na.rm = TRUE),
    geweke_z = max(abs(geweke_z), na.rm = TRUE),
    first_last_diff_sd = max(abs(first_last_diff_sd), na.rm = TRUE),
    flag = any(flag, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::select(dplyr::all_of(names(parameter_diagnostics)))

key_diagnostics <- parameter_diagnostics |>
  dplyr::filter(
    stringr::str_detect(group, "beta|variance|loading")
  ) |>
  dplyr::bind_rows(eta_lambda_summary) |>
  dplyr::arrange(dplyr::desc(flag), group, parameter)

gaussian_mu_idx <- selected_mu_indices((rowMeans(fit_ug$Mu.chain, na.rm = TRUE) + rowMeans(fit_mt$Mu_1.chain, na.rm = TRUE)) / 2)
binomial_mu_idx <- selected_mu_indices((rowMeans(fit_ub$Mu.chain, na.rm = TRUE) + rowMeans(fit_mt$Mu_2.chain, na.rm = TRUE)) / 2)

mu_diagnostics <- dplyr::bind_rows(
  mu_selected_diagnostics(fit_ug$Mu.chain, gaussian_mu_idx, "Univariate Gaussian Mu"),
  mu_selected_diagnostics(fit_mt$Mu_1.chain, gaussian_mu_idx, "Multi-type Gaussian Mu_1"),
  mu_difference_diagnostics(fit_mt$Mu_1.chain, fit_ug$Mu.chain, gaussian_mu_idx, "Gaussian Mu_1 - univariate Mu"),
  mu_selected_diagnostics(fit_ub$Mu.chain, binomial_mu_idx, "Univariate Bernoulli Mu"),
  mu_selected_diagnostics(fit_mt$Mu_2.chain, binomial_mu_idx, "Multi-type Bernoulli Mu_2"),
  mu_difference_diagnostics(fit_mt$Mu_2.chain, fit_ub$Mu.chain, binomial_mu_idx, "Bernoulli Mu_2 - univariate Mu")
) |>
  dplyr::mutate(
    diag_label = diag_label,
    seed = diag_seed,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    tau_2_init = tau_2_init,
    .before = 1
  )

utils::write.csv(parameter_diagnostics, parameter_diagnostics_file, row.names = FALSE)
utils::write.csv(key_diagnostics, key_diagnostics_file, row.names = FALSE)
utils::write.csv(mu_diagnostics, mu_diagnostics_file, row.names = FALSE)
utils::write.csv(region_summary, region_summary_file, row.names = FALSE)

plot_parameter_trace <- function(trace_df, title, y_lab = "Value") {
  ggplot(trace_df, aes(x = iteration, y = value, color = parameter)) +
    geom_line(linewidth = 0.25, alpha = 0.85) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = title, x = "Iteration", y = y_lab, color = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "bottom",
      legend.key.width = unit(0.8, "cm"),
      legend.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
}

trace_matrix_df <- function(chain_mat, prefix) {
  chain_mat <- as.matrix(chain_mat)
  param_names <- rownames(chain_mat)
  if (is.null(param_names)) {
    param_names <- paste0(prefix, "[", seq_len(nrow(chain_mat)), "]")
  }

  as.data.frame(t(chain_mat)) |>
    tibble::as_tibble(.name_repair = "minimal") |>
    stats::setNames(param_names) |>
    dplyr::mutate(iteration = dplyr::row_number()) |>
    tidyr::pivot_longer(
      cols = -iteration,
      names_to = "parameter",
      values_to = "value"
    ) |>
    dplyr::mutate(parameter = factor(parameter, levels = param_names))
}

trace_vector_df <- function(chain_vec, parameter) {
  tibble::tibble(
    iteration = seq_along(chain_vec),
    parameter = parameter,
    value = as.numeric(chain_vec)
  )
}

variance_trace_df <- dplyr::bind_rows(
  trace_vector_df(fit_ug$sig2.chain, "ug sigma2"),
  trace_vector_df(fit_ug$Sigma2_u.chain, "ug Sigma2_u"),
  trace_vector_df(fit_ub$Sigma2_u.chain, "ub Sigma2_u"),
  trace_vector_df(fit_mt$sig2.chain, "mt sigma2"),
  trace_vector_df(fit_mt$Sigma2_eta.chain, "mt Sigma2_eta"),
  trace_vector_df(fit_mt$Sigma2_lambda.chain, "mt Sigma2_lambda")
) |>
  dplyr::mutate(parameter = factor(parameter, levels = unique(parameter)))

parameter_trace_plot <- (
  plot_parameter_trace(trace_matrix_df(fit_ug$Beta.chain, "beta"), "Univariate Gaussian: beta") |
    plot_parameter_trace(trace_matrix_df(fit_ub$Beta.chain, "beta"), "Univariate Bernoulli: beta")
) / (
  plot_parameter_trace(trace_matrix_df(fit_mt$Beta_1.chain, "beta_1"), "Multi-type Gaussian: beta_1") |
    plot_parameter_trace(trace_matrix_df(fit_mt$Beta_2.chain, "beta_2"), "Multi-type Bernoulli: beta_2")
) / (
  plot_parameter_trace(variance_trace_df, "Variance parameters: sigma2 / Sigma2") |
    plot_parameter_trace(trace_vector_df(fit_mt$tau_2.chain, "tau_2"), "Multi-type: tau_2")
) +
  patchwork::plot_annotation(
    title = "Region Diagnostic Parameter Trace Plots",
    subtitle = paste0(
      "seed=", diag_seed,
      ", nburn=", nburn,
      ", nsim=", nsim,
      ", nthin=", nthin,
      ", tau_2_init=", tau_2_init
    )
  )

ggsave(
  parameter_trace_file,
  parameter_trace_plot,
  width = 12,
  height = 10,
  dpi = 300
)

fill_pal <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
ct_est <- 0.99

map_data_full <- puma_sf |>
  dplyr::left_join(region_summary, by = "PUMA")

inco_est_long <- map_data_full |>
  dplyr::select(
    PUMA,
    `Horvitz-Thompson` = ht_inco,
    `Multi-type` = multi_gaussian_mean,
    `Univariate (Gaussian)` = uni_gaussian_mean,
    geometry
  ) |>
  tidyr::pivot_longer(
    cols = c(`Horvitz-Thompson`, `Multi-type`, `Univariate (Gaussian)`),
    names_to = "source",
    values_to = "estimate"
  ) |>
  dplyr::mutate(
    source = factor(source, levels = c("Horvitz-Thompson", "Multi-type", "Univariate (Gaussian)")),
    estimate_plot = ifelse(estimate > stats::quantile(estimate, ct_est, na.rm = TRUE), NA_real_, estimate)
  )

pov_est_long <- map_data_full |>
  dplyr::select(
    PUMA,
    `Horvitz-Thompson` = ht_pov,
    `Multi-type` = multi_binomial_mean,
    `Univariate (Bernoulli)` = uni_binomial_mean,
    geometry
  ) |>
  tidyr::pivot_longer(
    cols = c(`Horvitz-Thompson`, `Multi-type`, `Univariate (Bernoulli)`),
    names_to = "source",
    values_to = "estimate"
  ) |>
  dplyr::mutate(
    source = factor(source, levels = c("Horvitz-Thompson", "Multi-type", "Univariate (Bernoulli)")),
    estimate_plot = ifelse(estimate > stats::quantile(estimate, ct_est, na.rm = TRUE), NA_real_, estimate)
  )

plot_gaus <- ggplot(inco_est_long) +
  geom_sf(aes(fill = estimate_plot), colour = NA) +
  facet_wrap(~source, nrow = 1) +
  scale_fill_gradientn(colours = fill_pal, name = "Transformed income", na.value = "grey90") +
  labs(title = "Transformed Income by PUMA", subtitle = "INCO = min-max scaled log income") +
  theme_minimal(base_size = 13)

p_sig_gaus <- ggplot(map_data_full, aes(x = multi_gaussian_var, y = uni_gaussian_var)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(
    x = expression(sigma^2[Multi-type]),
    y = expression(sigma^2[Univariate]),
    title = "Gaussian Variance Comparison"
  ) +
  theme_minimal(base_size = 13)

plot_bios <- ggplot(pov_est_long) +
  geom_sf(aes(fill = estimate_plot), colour = NA) +
  facet_wrap(~source, nrow = 1) +
  scale_fill_gradientn(colours = fill_pal, name = "Poverty rate", na.value = "grey90") +
  labs(title = "Poverty Rate by PUMA") +
  theme_minimal(base_size = 13)

p_sig_bios <- ggplot(map_data_full, aes(x = multi_binomial_var, y = uni_binomial_var)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(
    x = expression(sigma^2[Multi-type]),
    y = expression(sigma^2[Univariate]),
    title = "Bernoulli Variance Comparison"
  ) +
  theme_minimal(base_size = 13)

basis_main_plot <- ((plot_gaus | p_sig_gaus) + patchwork::plot_layout(widths = c(3, 1))) /
  ((plot_bios | p_sig_bios) + patchwork::plot_layout(widths = c(3, 1))) +
  patchwork::plot_layout(heights = c(1, 1))

ggsave(
  basis_main_file,
  basis_main_plot,
  width = 13,
  height = 8.5,
  dpi = 300
)

message("Saved diagnostic basis main plot to: ", basis_main_file)
message("Saved diagnostic parameter trace plot to: ", parameter_trace_file)
message("Saved MCMC key diagnostics to: ", key_diagnostics_file)
message("Saved full MCMC parameter diagnostics to: ", parameter_diagnostics_file)
message("Saved representative Mu diagnostics to: ", mu_diagnostics_file)
message(
  "Gaussian multi variance smaller: ",
  round(mean(region_summary$multi_gaussian_var < region_summary$uni_gaussian_var, na.rm = TRUE), 4)
)
message(
  "Bernoulli multi variance smaller: ",
  round(mean(region_summary$multi_binomial_var < region_summary$uni_binomial_var, na.rm = TRUE), 4)
)
message(
  "Median uni/multi variance ratio: Gaussian=",
  round(median(region_summary$gaussian_uni_over_multi_var_ratio, na.rm = TRUE), 4),
  ", Bernoulli=",
  round(median(region_summary$binomial_uni_over_multi_var_ratio, na.rm = TRUE), 4)
)
save(list = intersect(c("fit_ug", "fit_ub", "fit_mt", "res_ug", "res_ub", "res_mg", "res_mb", "region_summary", "parameter_diagnostics", "key_diagnostics", "mu_diagnostics", "diag_seed", "nburn", "nsim", "nthin", "tau_2_init", "prep_file", "basis_main_file", "parameter_trace_file", "key_diagnostics_file", "parameter_diagnostics_file", "mu_diagnostics_file", "region_summary_file"), ls()), file = file.path("data", paste0(diag_label, "_results.RData")))
