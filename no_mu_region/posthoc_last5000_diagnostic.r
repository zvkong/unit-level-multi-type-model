## Post-hoc diagnostic using only the last 5000 saved samples from a no-Mu
## region run. This does not rerun MCMC.

required_pkgs <- c(
  "dplyr", "tidyr", "tibble", "ggplot2", "RColorBrewer", "patchwork", "sf",
  "stringr", "coda"
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

source(file.path("no_mu_region", "functions_no_mu.r"))

source_result_file <- Sys.getenv(
  "REGION_POSTHOC_SOURCE",
  unset = file.path("data", "region_no_mu_tau0_burn5000_results.RData")
)
diag_label <- Sys.getenv(
  "REGION_POSTHOC_LABEL",
  unset = "region_no_mu_tau0_burn5000_last5000"
)
keep_start <- as.integer(Sys.getenv("REGION_POSTHOC_KEEP_START", unset = "5001"))
keep_end <- as.integer(Sys.getenv("REGION_POSTHOC_KEEP_END", unset = "10000"))

region_fig_dir <- file.path("results", "fig", "region")
region_table_dir <- file.path("results", "table", "region")
basis_main_file <- file.path(region_fig_dir, paste0(diag_label, "_basis_main.png"))
parameter_trace_file <- file.path(region_fig_dir, paste0(diag_label, "_parameter_traceplots.png"))
key_diagnostics_file <- file.path(region_table_dir, paste0(diag_label, "_mcmc_key_diagnostics.csv"))
parameter_diagnostics_file <- file.path(region_table_dir, paste0(diag_label, "_mcmc_parameter_diagnostics.csv"))
region_summary_file <- file.path(region_table_dir, paste0(diag_label, "_region_summary.csv"))

dir.create(region_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(region_table_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading saved run from: ", source_result_file)
run_env <- new.env(parent = emptyenv())
load(source_result_file, envir = run_env)

for (nm in c("fit_ug", "fit_ub", "fit_mt")) {
  if (!exists(nm, envir = run_env)) {
    stop("Missing object in source result file: ", nm, call. = FALSE)
  }
}

fit_ug <- get("fit_ug", envir = run_env)
fit_ub <- get("fit_ub", envir = run_env)
fit_mt <- get("fit_mt", envir = run_env)
source_diag_label <- if (exists("diag_label", envir = run_env)) get("diag_label", envir = run_env) else NA_character_
diag_seed <- if (exists("diag_seed", envir = run_env)) get("diag_seed", envir = run_env) else 20260706L
source_nburn <- if (exists("nburn", envir = run_env)) get("nburn", envir = run_env) else NA_integer_
source_nsim <- if (exists("nsim", envir = run_env)) get("nsim", envir = run_env) else NA_integer_
source_nthin <- if (exists("nthin", envir = run_env)) get("nthin", envir = run_env) else NA_integer_
tau_2_init <- if (exists("tau_2_init", envir = run_env)) get("tau_2_init", envir = run_env) else NA_real_
prep_file <- if (exists("prep_file", envir = run_env)) {
  get("prep_file", envir = run_env)
} else {
  file.path("data", "region_basis_tau_in_gauss_results.RData")
}
rm(run_env)
gc()

n_saved_total <- length(fit_mt$tau_2.chain)
if (keep_start < 1 || keep_end > n_saved_total || keep_start > keep_end) {
  stop(
    "Invalid keep range: ", keep_start, ":", keep_end,
    " for ", n_saved_total, " saved samples.",
    call. = FALSE
  )
}
keep_idx <- keep_start:keep_end

subset_matrix_cols <- function(x, idx) {
  if (is.null(x)) return(NULL)
  x <- as.matrix(x)
  x[, idx, drop = FALSE]
}

subset_vector <- function(x, idx) {
  if (is.null(x)) return(NULL)
  as.numeric(x)[idx]
}

fit_ug$Beta.chain <- subset_matrix_cols(fit_ug$Beta.chain, keep_idx)
fit_ug$U.chain <- subset_matrix_cols(fit_ug$U.chain, keep_idx)
fit_ug$Preds <- subset_matrix_cols(fit_ug$Preds, keep_idx)
fit_ug$sig2.chain <- subset_vector(fit_ug$sig2.chain, keep_idx)
fit_ug$Sigma2_u.chain <- subset_vector(fit_ug$Sigma2_u.chain, keep_idx)

fit_ub$Beta.chain <- subset_matrix_cols(fit_ub$Beta.chain, keep_idx)
fit_ub$U.chain <- subset_matrix_cols(fit_ub$U.chain, keep_idx)
fit_ub$Preds <- subset_matrix_cols(fit_ub$Preds, keep_idx)
fit_ub$logit_bios.chain <- subset_matrix_cols(fit_ub$logit_bios.chain, keep_idx)
fit_ub$Sigma2_u.chain <- subset_vector(fit_ub$Sigma2_u.chain, keep_idx)

fit_mt$Beta_1.chain <- subset_matrix_cols(fit_mt$Beta_1.chain, keep_idx)
fit_mt$Beta_2.chain <- subset_matrix_cols(fit_mt$Beta_2.chain, keep_idx)
fit_mt$eta.chain <- subset_matrix_cols(fit_mt$eta.chain, keep_idx)
fit_mt$lambda.chain <- subset_matrix_cols(fit_mt$lambda.chain, keep_idx)
fit_mt$preds_gaus.chain <- subset_matrix_cols(fit_mt$preds_gaus.chain, keep_idx)
fit_mt$preds_bios.chain <- subset_matrix_cols(fit_mt$preds_bios.chain, keep_idx)
fit_mt$logit_bios.chain <- subset_matrix_cols(fit_mt$logit_bios.chain, keep_idx)
fit_mt$sig2.chain <- subset_vector(fit_mt$sig2.chain, keep_idx)
fit_mt$Sigma2_eta.chain <- subset_vector(fit_mt$Sigma2_eta.chain, keep_idx)
fit_mt$Sigma2_lambda.chain <- subset_vector(fit_mt$Sigma2_lambda.chain, keep_idx)
fit_mt$tau_2.chain <- subset_vector(fit_mt$tau_2.chain, keep_idx)

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
area_levels <- if (exists("area_levels", envir = prep_env)) {
  get("area_levels", envir = prep_env)
} else {
  rownames(basis_mat)
}
rm(prep_env)
gc()

true_mean_dummy <- stats::setNames(rep(NA_real_, length(area_levels)), area_levels)

set.seed(as.integer(diag_seed) + 200000L)
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
    source_diag_label = source_diag_label,
    source_nburn = source_nburn,
    source_nsim = source_nsim,
    source_nthin = source_nthin,
    kept_saved_start = keep_start,
    kept_saved_end = keep_end,
    effective_nburn = source_nburn + keep_start - 1,
    tau_2_init = tau_2_init,
    .before = 1
  )

eta_lambda_summary <- parameter_diagnostics |>
  dplyr::filter(group %in% c("Multi-type shared eta", "Multi-type residual lambda")) |>
  dplyr::group_by(group) |>
  dplyr::summarise(
    dplyr::across(
      c(
        diag_label, source_diag_label, source_nburn, source_nsim, source_nthin,
        kept_saved_start, kept_saved_end, effective_nburn, tau_2_init
      ),
      dplyr::first
    ),
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
  dplyr::filter(stringr::str_detect(group, "beta|variance|loading")) |>
  dplyr::bind_rows(eta_lambda_summary) |>
  dplyr::arrange(dplyr::desc(flag), group, parameter)

utils::write.csv(parameter_diagnostics, parameter_diagnostics_file, row.names = FALSE)
utils::write.csv(key_diagnostics, key_diagnostics_file, row.names = FALSE)
utils::write.csv(region_summary, region_summary_file, row.names = FALSE)

plot_parameter_trace <- function(trace_df, title, y_lab = "Value") {
  ggplot(trace_df, aes(x = iteration, y = value, color = parameter)) +
    geom_line(linewidth = 0.25, alpha = 0.85) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = title, x = "Saved sample index after post-hoc burn-in", y = y_lab, color = NULL) +
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
    title = "Region Post-hoc Parameter Trace Plots",
    subtitle = paste0(
      "source=", source_diag_label,
      ", kept saved samples=", keep_start, ":", keep_end,
      ", effective nburn=", source_nburn + keep_start - 1,
      ", tau_2_init=", tau_2_init
    )
  )

ggsave(parameter_trace_file, parameter_trace_plot, width = 12, height = 10, dpi = 300)

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

ggsave(basis_main_file, basis_main_plot, width = 13, height = 8.5, dpi = 300)

save(
  list = c(
    "fit_ug", "fit_ub", "fit_mt", "res_ug", "res_ub", "res_mg", "res_mb",
    "region_summary", "parameter_diagnostics", "key_diagnostics",
    "source_result_file", "source_diag_label", "source_nburn", "source_nsim",
    "source_nthin", "keep_start", "keep_end", "tau_2_init", "prep_file",
    "basis_main_file", "parameter_trace_file", "key_diagnostics_file",
    "parameter_diagnostics_file", "region_summary_file"
  ),
  file = file.path("data", paste0(diag_label, "_results.RData"))
)

message("Saved post-hoc basis main plot to: ", basis_main_file)
message("Saved post-hoc parameter trace plot to: ", parameter_trace_file)
message("Saved post-hoc MCMC key diagnostics to: ", key_diagnostics_file)
message("Saved post-hoc full MCMC parameter diagnostics to: ", parameter_diagnostics_file)
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
message("Key Geweke diagnostics:")
print(
  key_diagnostics |>
    dplyr::select(group, parameter, n_saved, mean, lag1_acf, ess, geweke_z, first_last_diff_sd, flag),
  n = Inf
)
