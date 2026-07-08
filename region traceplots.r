required_pkgs <- c("ggplot2", "dplyr", "tidyr", "patchwork", "RColorBrewer")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required package: ", pkg, call. = FALSE)
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

region_result_file <- Sys.getenv(
  "REGION_TRACE_RESULT_FILE",
  unset = file.path("data", "region_basis_tau_in_gauss_results_previous.RData")
)
region_traceplot_dir <- file.path("results", "fig", "region")
region_traceplot_label <- tools::file_path_sans_ext(basename(region_result_file))
region_traceplot_file <- file.path(region_traceplot_dir, paste0(region_traceplot_label, "_mu_traceplots.png"))

dir.create(region_traceplot_dir, recursive = TRUE, showWarnings = FALSE)
load(region_result_file)

if (!exists("fit_ug") && exists("unis_wage")) {
  fit_ug <- unis_wage
}
if (!exists("fit_ub") && exists("unis_pov")) {
  fit_ub <- unis_pov
}
if (!exists("fit_mt") && exists("mult_br")) {
  fit_mt <- mult_br
}

required_objects <- c("fit_ug", "fit_ub", "fit_mt")
missing_objects <- setdiff(required_objects, ls())
if (length(missing_objects) > 0) {
  stop(
    "Missing model object(s) in ", region_result_file, ": ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

required_chains <- list(
  fit_ug_Mu = fit_ug$Mu.chain,
  fit_ub_Mu = fit_ub$Mu.chain,
  fit_mt_Mu_1 = fit_mt$Mu_1.chain,
  fit_mt_Mu_2 = fit_mt$Mu_2.chain
)

missing_chains <- names(required_chains)[vapply(required_chains, is.null, logical(1))]
if (length(missing_chains) > 0) {
  stop(
    "Missing Mu chain(s): ",
    paste(missing_chains, collapse = ", "),
    call. = FALSE
  )
}

select_mu_indices <- function(mu_mean, n_select = 5) {
  probs <- seq(0.1, 0.9, length.out = n_select)
  targets <- as.numeric(stats::quantile(mu_mean, probs = probs, na.rm = TRUE))
  unique(vapply(targets, function(target) {
    which.min(abs(mu_mean - target))
  }, integer(1)))
}

mu_trace_df <- function(chain_mat, selected_idx, label_prefix) {
  selected_idx <- as.integer(selected_idx)
  selected_mat <- as.matrix(chain_mat[selected_idx, , drop = FALSE])
  rownames(selected_mat) <- paste0(label_prefix, " q", seq_along(selected_idx))

  as.data.frame(t(selected_mat)) |>
    tibble::as_tibble(.name_repair = "minimal") |>
    stats::setNames(rownames(selected_mat)) |>
    dplyr::mutate(iteration = dplyr::row_number()) |>
    tidyr::pivot_longer(
      cols = -iteration,
      names_to = "parameter",
      values_to = "value"
    ) |>
    dplyr::mutate(parameter = factor(parameter, levels = rownames(selected_mat)))
}

mu_difference_trace_df <- function(multi_chain, uni_chain, selected_idx, label_prefix) {
  selected_idx <- as.integer(selected_idx)
  diff_mat <- as.matrix(multi_chain[selected_idx, , drop = FALSE] - uni_chain[selected_idx, , drop = FALSE])
  rownames(diff_mat) <- paste0(label_prefix, " q", seq_along(selected_idx))

  as.data.frame(t(diff_mat)) |>
    tibble::as_tibble(.name_repair = "minimal") |>
    stats::setNames(rownames(diff_mat)) |>
    dplyr::mutate(iteration = dplyr::row_number()) |>
    tidyr::pivot_longer(
      cols = -iteration,
      names_to = "parameter",
      values_to = "value"
    ) |>
    dplyr::mutate(parameter = factor(parameter, levels = rownames(diff_mat)))
}

plot_mu_trace <- function(trace_df, title, y_lab = "Mu") {
  ggplot(trace_df, aes(x = iteration, y = value, color = parameter)) +
    geom_line(linewidth = 0.25, alpha = 0.85) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = title, x = "Iteration", y = y_lab, color = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "bottom",
      legend.key.width = unit(0.8, "cm"),
      legend.text = element_text(size = 7),
      panel.grid.minor = element_blank()
    )
}

gaussian_mean_score <- (rowMeans(fit_ug$Mu.chain, na.rm = TRUE) + rowMeans(fit_mt$Mu_1.chain, na.rm = TRUE)) / 2
bernoulli_mean_score <- (rowMeans(fit_ub$Mu.chain, na.rm = TRUE) + rowMeans(fit_mt$Mu_2.chain, na.rm = TRUE)) / 2

gaussian_idx <- select_mu_indices(gaussian_mean_score)
bernoulli_idx <- select_mu_indices(bernoulli_mean_score)

p_ug_mu <- mu_trace_df(fit_ug$Mu.chain, gaussian_idx, "mu") |>
  plot_mu_trace("Univariate Gaussian: Mu")

p_mg_mu <- mu_trace_df(fit_mt$Mu_1.chain, gaussian_idx, "mu") |>
  plot_mu_trace("Multi-type Gaussian: Mu_1")

p_gaussian_diff <- mu_difference_trace_df(fit_mt$Mu_1.chain, fit_ug$Mu.chain, gaussian_idx, "diff") |>
  plot_mu_trace("Gaussian: Multi-type Mu_1 - Univariate Mu", y_lab = "Difference")

p_ub_mu <- mu_trace_df(fit_ub$Mu.chain, bernoulli_idx, "mu") |>
  plot_mu_trace("Univariate Bernoulli: Mu")

p_mb_mu <- mu_trace_df(fit_mt$Mu_2.chain, bernoulli_idx, "mu") |>
  plot_mu_trace("Multi-type Bernoulli: Mu_2")

p_bernoulli_diff <- mu_difference_trace_df(fit_mt$Mu_2.chain, fit_ub$Mu.chain, bernoulli_idx, "diff") |>
  plot_mu_trace("Bernoulli: Multi-type Mu_2 - Univariate Mu", y_lab = "Difference")

traceplot_grid <- (
  p_ug_mu | p_mg_mu
) / (
  p_gaussian_diff | p_ub_mu
) / (
  p_mb_mu | p_bernoulli_diff
) +
  patchwork::plot_annotation(
    title = "Region Basis MCMC Mu Trace Plots",
    subtitle = "Five representative units selected by posterior mean quantiles"
  )

ggsave(
  region_traceplot_file,
  traceplot_grid,
  width = 12,
  height = 10,
  dpi = 300
)

message("Saved region Mu trace plots to: ", region_traceplot_file)
message("Gaussian selected Mu rows: ", paste(gaussian_idx, collapse = ", "))
message("Bernoulli selected Mu rows: ", paste(bernoulli_idx, collapse = ", "))
