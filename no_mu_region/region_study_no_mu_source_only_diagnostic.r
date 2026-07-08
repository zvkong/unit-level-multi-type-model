## Region study diagnostic run without storing Mu chains, source-data only.
##
## This script does not load the prepared data/*.RData region-input file. It
## rebuilds pums, pcells, and a region basis directly from source data/IL23.csv.
## Because no tigris/ACS objects are reused or fetched, the basis is a
## source-only PUMA contrast basis rather than a shapefile-adjacency basis.

required_pkgs <- c(
  "readr", "dplyr", "tidyr", "tibble", "ggplot2", "RColorBrewer", "patchwork", "stringr",
  "coda"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
  library(stringr)
  library(coda)
})

source(file.path("no_mu_region", "functions_no_mu.r"))

diag_seed <- as.integer(Sys.getenv("REGION_DIAG_SEED", unset = "20260706"))
nburn <- as.integer(Sys.getenv("REGION_DIAG_NBURN", unset = "5000"))
nsim <- as.integer(Sys.getenv("REGION_DIAG_NSIM", unset = "10000"))
nthin <- as.integer(Sys.getenv("REGION_DIAG_NTHIN", unset = "1"))
tau_2_init <- as.numeric(Sys.getenv("REGION_DIAG_TAU2_INIT", unset = "0"))
diag_label <- Sys.getenv("REGION_DIAG_LABEL", unset = "region_no_mu_source_only_tau0_burn5000")

if (nsim %% nthin != 0) {
  stop("REGION_DIAG_NSIM must be divisible by REGION_DIAG_NTHIN.", call. = FALSE)
}

data_file <- Sys.getenv(
  "REGION_SOURCE_DATA_FILE",
  unset = file.path("source data", "IL23.csv")
)
basis_source <- "source-only PUMA contrast basis"

region_fig_dir <- file.path("results", "fig", "region")
region_table_dir <- file.path("results", "table", "region")
summary_plot_file <- file.path(region_fig_dir, paste0(diag_label, "_summary.png"))
parameter_trace_file <- file.path(region_fig_dir, paste0(diag_label, "_parameter_traceplots.png"))
key_diagnostics_file <- file.path(region_table_dir, paste0(diag_label, "_mcmc_key_diagnostics.csv"))
parameter_diagnostics_file <- file.path(region_table_dir, paste0(diag_label, "_mcmc_parameter_diagnostics.csv"))
region_summary_file <- file.path(region_table_dir, paste0(diag_label, "_region_summary.csv"))

dir.create(region_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(region_table_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading source data from: ", data_file)
if (!file.exists(data_file)) {
  stop("Source data file not found: ", data_file, call. = FALSE)
}

source_required_cols <- c("PUMA", "PWGTP", "AGEP", "SEX", "POVPIP", "PINCP", "SCHL")
pums_raw <- readr::read_csv(data_file, show_col_types = FALSE)
missing_source_cols <- setdiff(source_required_cols, names(pums_raw))
if (length(missing_source_cols) > 0) {
  stop("Missing required source column(s): ", paste(missing_source_cols, collapse = ", "), call. = FALSE)
}

pums <- pums_raw |>
  dplyr::select(dplyr::all_of(source_required_cols)) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    PUMA = sprintf("%05d", as.integer(PUMA)),
    PWGTP = as.numeric(PWGTP),
    AGEP = as.numeric(AGEP),
    POVPIP = as.numeric(POVPIP),
    PINCP = as.numeric(PINCP),
    SCHL = as.numeric(SCHL),
    SEX = as.character(SEX),
    LOGINCO = dplyr::if_else(PINCP > 0, log(PINCP), NA_real_),
    POV = dplyr::if_else(POVPIP <= 100, 1, 0),
    BACH = dplyr::if_else(SCHL >= 21, "1", "0")
  ) |>
  dplyr::filter(AGEP >= 18, is.finite(LOGINCO), is.finite(PWGTP), PWGTP > 0) |>
  dplyr::mutate(
    SEX = factor(SEX, levels = sort(unique(SEX))),
    BACH = factor(BACH, levels = c("0", "1")),
    INCO = (LOGINCO - min(LOGINCO)) / (max(LOGINCO) - min(LOGINCO)),
    scaledWGT = PWGTP * dplyr::n() / sum(PWGTP)
  ) |>
  dplyr::select(PUMA, PWGTP, SEX, BACH, INCO, POV, scaledWGT) |>
  dplyr::arrange(PUMA, SEX, BACH)

if (dplyr::n_distinct(pums$PUMA) < 2) {
  stop("At least two PUMA regions are required to build a source-only region basis.", call. = FALSE)
}

area_levels <- sort(unique(pums$PUMA))

build_puma_contrast_basis <- function(area_levels) {
  n_area <- length(area_levels)
  centered_identity <- scale(diag(n_area), center = TRUE, scale = FALSE)
  basis_qr <- qr(centered_identity)
  basis_mat <- qr.Q(basis_qr, complete = FALSE)[, seq_len(basis_qr$rank), drop = FALSE]
  rownames(basis_mat) <- area_levels
  colnames(basis_mat) <- paste0("bf", seq_len(ncol(basis_mat)))
  basis_mat
}

basis_mat <- build_puma_contrast_basis(area_levels)

cell_pop <- pums |>
  dplyr::group_by(PUMA, SEX, BACH) |>
  dplyr::summarise(popsize = sum(PWGTP), .groups = "drop")

pcells <- tidyr::expand_grid(
  PUMA = area_levels,
  SEX = levels(pums$SEX),
  BACH = levels(pums$BACH)
) |>
  dplyr::mutate(
    SEX = factor(SEX, levels = levels(pums$SEX)),
    BACH = factor(BACH, levels = levels(pums$BACH))
  ) |>
  dplyr::left_join(cell_pop, by = c("PUMA", "SEX", "BACH")) |>
  dplyr::mutate(popsize = dplyr::coalesce(popsize, 0)) |>
  dplyr::arrange(factor(PUMA, levels = area_levels), SEX, BACH)

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
  b_eps = 0.1,
  save_mu = FALSE
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
  b = 0.1,
  save_mu = FALSE
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
  blambda = 1,
  save_mu = FALSE
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
    data_file = data_file,
    basis_source = basis_source,
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
  matrix_diagnostics(fit_ug$U.chain, "Univariate Gaussian region U", "U"),
  matrix_diagnostics(fit_ub$U.chain, "Univariate Bernoulli region U", "U"),
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
    data_file = data_file,
    basis_source = basis_source,
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
    data_file = dplyr::first(data_file),
    basis_source = dplyr::first(basis_source),
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

utils::write.csv(parameter_diagnostics, parameter_diagnostics_file, row.names = FALSE)
utils::write.csv(key_diagnostics, key_diagnostics_file, row.names = FALSE)
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

ordered_region_summary <- region_summary |>
  dplyr::arrange(factor(PUMA, levels = area_levels)) |>
  dplyr::mutate(puma_order = dplyr::row_number())

estimate_long <- dplyr::bind_rows(
  ordered_region_summary |>
    dplyr::transmute(
      puma_order,
      PUMA,
      outcome = "Transformed income",
      `Horvitz-Thompson` = ht_inco,
      `Multi-type` = multi_gaussian_mean,
      `Univariate` = uni_gaussian_mean
    ),
  ordered_region_summary |>
    dplyr::transmute(
      puma_order,
      PUMA,
      outcome = "Poverty rate",
      `Horvitz-Thompson` = ht_pov,
      `Multi-type` = multi_binomial_mean,
      `Univariate` = uni_binomial_mean
    )
) |>
  tidyr::pivot_longer(
    cols = c(`Horvitz-Thompson`, `Multi-type`, `Univariate`),
    names_to = "source",
    values_to = "estimate"
  ) |>
  dplyr::mutate(
    outcome = factor(outcome, levels = c("Transformed income", "Poverty rate")),
    source = factor(source, levels = c("Horvitz-Thompson", "Multi-type", "Univariate"))
  )

variance_long <- dplyr::bind_rows(
  ordered_region_summary |>
    dplyr::transmute(
      outcome = "Gaussian variance",
      multi_var = multi_gaussian_var,
      uni_var = uni_gaussian_var
    ),
  ordered_region_summary |>
    dplyr::transmute(
      outcome = "Bernoulli variance",
      multi_var = multi_binomial_var,
      uni_var = uni_binomial_var
    )
)

estimate_plot <- ggplot(estimate_long, aes(x = puma_order, y = estimate, color = source)) +
  geom_line(linewidth = 0.35, alpha = 0.85) +
  geom_point(size = 0.7, alpha = 0.8) +
  facet_grid(outcome ~ ., scales = "free_y") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "PUMA order", y = "Estimate", color = NULL, title = "Region Estimates by PUMA") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

variance_plot <- ggplot(variance_long, aes(x = multi_var, y = uni_var)) +
  geom_point(size = 1.1, alpha = 0.75) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_wrap(~outcome, scales = "free") +
  labs(
    x = expression(sigma^2[Multi-type]),
    y = expression(sigma^2[Univariate]),
    title = "Posterior Variance Comparison"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

summary_plot <- (estimate_plot / variance_plot) +
  patchwork::plot_layout(heights = c(2, 1)) +
  patchwork::plot_annotation(
    title = "Source-Only Region Diagnostic Summary",
    subtitle = paste0(
      "data=", data_file,
      "; basis=", basis_source,
      "; seed=", diag_seed,
      ", nburn=", nburn,
      ", nsim=", nsim,
      ", nthin=", nthin
    )
  )

ggsave(
  summary_plot_file,
  summary_plot,
  width = 13,
  height = 8.5,
  dpi = 300
)

message("Saved source-only diagnostic summary plot to: ", summary_plot_file)
message("Saved diagnostic parameter trace plot to: ", parameter_trace_file)
message("Saved MCMC key diagnostics to: ", key_diagnostics_file)
message("Saved full MCMC parameter diagnostics to: ", parameter_diagnostics_file)
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
save(
  list = intersect(
    c(
      "fit_ug", "fit_ub", "fit_mt", "res_ug", "res_ub", "res_mg", "res_mb",
      "pums", "pcells", "basis_mat", "area_levels", "region_summary",
      "parameter_diagnostics", "key_diagnostics", "diag_seed", "nburn", "nsim",
      "nthin", "tau_2_init", "data_file", "basis_source", "summary_plot_file",
      "parameter_trace_file", "key_diagnostics_file", "parameter_diagnostics_file",
      "region_summary_file"
    ),
    ls()
  ),
  file = file.path("data", paste0(diag_label, "_results.RData"))
)
