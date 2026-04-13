source("basis version/packages.r")
source("basis version/functions.r")

load("data/region_basis_tau_in_gauss_results.RData")
if (!exists("fit_ug") && exists("unis_wage")) {
  fit_ug <- unis_wage
} else if (exists("fit_ug") && !exists("unis_wage")) {
  unis_wage <- fit_ug
}

if (!exists("fit_ub") && exists("unis_pov")) {
  fit_ub <- unis_pov
} else if (exists("fit_ub") && !exists("unis_pov")) {
  unis_pov <- fit_ub
}

if (!exists("fit_mt") && exists("mult_br")) {
  fit_mt <- mult_br
} else if (exists("fit_mt") && !exists("mult_br")) {
  mult_br <- fit_mt
}
puma_levels <- sort(unique(pcells$PUMA))
true_mean_dummy <- stats::setNames(rep(NA_real_, length(puma_levels)), puma_levels)

puma_stats <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(
    n_sample = dplyr::n(),
    ht_inco = stats::weighted.mean(INCO, PWGTP),
    ht_pov = stats::weighted.mean(POV, PWGTP),
    .groups = "drop"
  ) |>
  dplyr::arrange(PUMA)

res_ug <- gaus_post(fit_ug$Preds, fit_ug$sig2.chain, true_mean_dummy, pcells$PUMA, pcells$popsize)
res_mg <- gaus_post(fit_mt$preds_gaus.chain, fit_mt$sig2.chain, true_mean_dummy, pcells$PUMA, pcells$popsize)
res_ub <- bios_post(fit_ub$Preds, true_mean_dummy, pcells$PUMA, pcells$popsize)
res_mb <- bios_post(fit_mt$preds_bios.chain, true_mean_dummy, pcells$PUMA, pcells$popsize)

region_summary <- puma_stats |>
  dplyr::left_join(
    tibble::tibble(
      PUMA = puma_levels,
      multi_gaussian_mean = res_mg$est, uni_gaussian_mean = res_ug$est,
      multi_gaussian_var = as.numeric(res_mg$sigma2), uni_gaussian_var = as.numeric(res_ug$sigma2),
      multi_gaussian_lb = res_mg$lb, multi_gaussian_ub = res_mg$ub,
      uni_gaussian_lb = res_ug$lb, uni_gaussian_ub = res_ug$ub,
      
      multi_binomial_mean = res_mb$est, uni_binomial_mean = res_ub$est,
      multi_binomial_var = as.numeric(res_mb$sigma2), uni_binomial_var = as.numeric(res_ub$sigma2),
      multi_binomial_lb = res_mb$lb, multi_binomial_ub = res_mb$ub,
      uni_binomial_lb = res_ub$lb, uni_binomial_ub = res_ub$ub
    ), by = "PUMA"
  ) |>
  dplyr::mutate(
    gaussian_uni_over_multi_var_ratio = uni_gaussian_var / multi_gaussian_var,
    binomial_uni_over_multi_var_ratio = uni_binomial_var / multi_binomial_var
  )

variance_summary <- tibble::tibble(
  Response = c("Gaussian", "Bernoulli"),
  mean_ht = c(mean(region_summary$ht_inco, na.rm = TRUE), mean(region_summary$ht_pov, na.rm = TRUE)),
  mean_uni = c(mean(region_summary$uni_gaussian_mean, na.rm = TRUE), mean(region_summary$uni_binomial_mean, na.rm = TRUE)),
  mean_multi = c(mean(region_summary$multi_gaussian_mean, na.rm = TRUE), mean(region_summary$multi_binomial_mean, na.rm = TRUE)),
  mean_uni_var = c(mean(region_summary$uni_gaussian_var, na.rm = TRUE), mean(region_summary$uni_binomial_var, na.rm = TRUE)),
  mean_multi_var = c(mean(region_summary$multi_gaussian_var, na.rm = TRUE), mean(region_summary$multi_binomial_var, na.rm = TRUE)),
  prop_multi_var_smaller = c(
    mean(region_summary$multi_gaussian_var < region_summary$uni_gaussian_var, na.rm = TRUE),
    mean(region_summary$multi_binomial_var < region_summary$uni_binomial_var, na.rm = TRUE)
  ),
  median_uni_over_multi_var_ratio = c(
    median(region_summary$gaussian_uni_over_multi_var_ratio, na.rm = TRUE),
    median(region_summary$binomial_uni_over_multi_var_ratio, na.rm = TRUE)
  )
)

fill_pal <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
ct_est <- 0.99

map_data_full <- puma_sf |>
  dplyr::left_join(region_summary, by = "PUMA")

inco_est_long <- map_data_full |>
  dplyr::select(PUMA, `Horvitz–Thompson` = ht_inco, `Multi-type` = multi_gaussian_mean, `Univariate (Gaussian)` = uni_gaussian_mean, geometry) |>
  tidyr::pivot_longer(cols = c(`Horvitz–Thompson`, `Multi-type`, `Univariate (Gaussian)`), names_to = "source", values_to = "estimate") |>
  dplyr::mutate(
    source = factor(source, levels = c("Horvitz–Thompson", "Multi-type", "Univariate (Gaussian)")),
    estimate_plot = ifelse(estimate > stats::quantile(estimate, ct_est, na.rm = TRUE), NA_real_, estimate)
  )

plot_gaus <- ggplot(inco_est_long) +
  geom_sf(aes(fill = estimate_plot), colour = NA) +
  facet_wrap(~source, nrow = 1) +
  scale_fill_gradientn(colours = fill_pal, name = "Transformed income", na.value = "grey90") +
  labs(title = "Transformed Income by PUMA", subtitle = "INCO = scaled log income") +
  theme_minimal()

p_sig_gaus <- ggplot(map_data_full, aes(x = multi_gaussian_var, y = uni_gaussian_var)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = expression(sigma^2[Multi-type]), y = expression(sigma^2[Univariate]), title = "Gaussian Variance Comparison") +
  theme_minimal()

p_gaussian_base <- (plot_gaus | p_sig_gaus) + patchwork::plot_layout(widths = c(3, 1))

pov_est_long <- map_data_full |>
  dplyr::select(PUMA, `Horvitz–Thompson` = ht_pov, `Multi-type` = multi_binomial_mean, `Univariate (Bernoulli)` = uni_binomial_mean, geometry) |>
  tidyr::pivot_longer(cols = c(`Horvitz–Thompson`, `Multi-type`, `Univariate (Bernoulli)`), names_to = "source", values_to = "estimate") |>
  dplyr::mutate(
    source = factor(source, levels = c("Horvitz–Thompson", "Multi-type", "Univariate (Bernoulli)")),
    estimate_plot = ifelse(estimate > stats::quantile(estimate, ct_est, na.rm = TRUE), NA_real_, estimate)
  )

plot_bios <- ggplot(pov_est_long) +
  geom_sf(aes(fill = estimate_plot), colour = NA) +
  facet_wrap(~source, nrow = 1) +
  scale_fill_gradientn(colours = fill_pal, name = "Poverty rate", na.value = "grey90") +
  labs(title = "Poverty Rate by PUMA") +
  theme_minimal()

p_sig_bios <- ggplot(map_data_full, aes(x = multi_binomial_var, y = uni_binomial_var)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = expression(sigma^2[Multi-type]), y = expression(sigma^2[Univariate]), title = "Bernoulli Variance Comparison") +
  theme_minimal()

p_binomial_base <- (plot_bios | p_sig_bios) + patchwork::plot_layout(widths = c(3, 1))

plot_advanced_metrics <- function(region_df, var_uni_col, var_multi_col, lb_uni, ub_uni, est_uni, lb_multi, ub_multi, est_multi, response_name) {
  
  puma_order <- region_df |> dplyr::arrange(n_sample) |> dplyr::pull(PUMA)
  
  cat_data <- dplyr::bind_rows(
    tibble::tibble(PUMA = region_df$PUMA, Model = "Univariate", est = region_df[[est_uni]], lb = region_df[[lb_uni]], ub = region_df[[ub_uni]]),
    tibble::tibble(PUMA = region_df$PUMA, Model = "Multi-type", est = region_df[[est_multi]], lb = region_df[[lb_multi]], ub = region_df[[ub_multi]])
  ) |> dplyr::mutate(PUMA = factor(PUMA, levels = puma_order))
  
  p_cat <- ggplot(cat_data, aes(x = PUMA, y = est, color = Model, shape = Model)) +
    geom_point(position = position_dodge(width = 0.6), size = 1.5) +
    geom_errorbar(aes(ymin = lb, ymax = ub), position = position_dodge(width = 0.6), width = 0, alpha = 0.7, linewidth = 0.8) +
    scale_color_manual(values = c("Multi-type" = "#1f78b4", "Univariate" = "#e31a1c")) +
    labs(title = paste("95% Credible Intervals by Region", response_name), subtitle = "Regions ordered by increasing sample size", x = "Regions (PUMA)", y = "Estimate") +
    theme_minimal(base_size = 14) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = "bottom")
  
  p_rat <- ggplot(region_df, aes(x = n_sample, y = .data[[var_multi_col]] / .data[[var_uni_col]])) +
    geom_hline(yintercept = 1, color = "black", linewidth = 1) +
    geom_point(alpha = 0.7, color = "#2c3e50", size = 2.5) +
    geom_smooth(method = "loess", se = FALSE, color = "#e74c3c", linetype = "dashed", linewidth = 1.2) +
    scale_x_log10() +
    labs(title = paste("Variance Ratio vs Sample Size", response_name), subtitle = "Values below line indicate variance reduction", x = "Sample Size (log scale)", y = "Variance Ratio (Multi/Uni)") +
    theme_minimal(base_size = 14)
  
  p_45 <- ggplot(region_df, aes(x = .data[[var_multi_col]], y = .data[[var_uni_col]], color = log(n_sample))) +
    geom_abline(slope = 1, intercept = 0, color = "#e74c3c", linetype = "dashed", linewidth = 1) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_viridis_c(name = "Log(Size)", option = "plasma") +
    scale_x_log10() + scale_y_log10() + coord_fixed(ratio = 1) +
    labs(title = paste("Posterior Variance Comparison", response_name), x = expression(Log~Var(Multi-type)), y = expression(Log~Var(Univariate))) +
    theme_minimal(base_size = 14)
  
  list(cat = p_cat, comb = (p_rat + p_45 + patchwork::plot_layout(widths = c(1, 1))))
}

adv_gaus <- plot_advanced_metrics(region_summary, "uni_gaussian_var", "multi_gaussian_var", "uni_gaussian_lb", "uni_gaussian_ub", "uni_gaussian_mean", "multi_gaussian_lb", "multi_gaussian_ub", "multi_gaussian_mean", "(Gaussian)")
adv_bios <- plot_advanced_metrics(region_summary, "uni_binomial_var", "multi_binomial_var", "uni_binomial_lb", "uni_binomial_ub", "uni_binomial_mean", "multi_binomial_lb", "multi_binomial_ub", "multi_binomial_mean", "(Binomial)")

dir.create("data", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

save.image("data/region_result_complete_basis.RData")
readr::write_csv(region_summary, "basis_region_summary.csv")
readr::write_csv(variance_summary, "basis_region_variance_summary.csv")

ggsave(file.path("figs", "p_gaussian_tau_in_gauss_basis.png"), p_gaussian_base, width = 10.5, height = 3, dpi = 300)
ggsave(file.path("figs", "p_binomial_tau_in_gauss_basis.png"), p_binomial_base, width = 10.5, height = 3, dpi = 300)
## save name as basis
ggsave(file.path("figs", "p_gaussian_base_basis.png"), adv_gaus$comb, width = 10.5, height = 6, dpi = 300)
ggsave(file.path("figs", "p_binomial_base_basis.png"), adv_bios$comb, width = 10.5, height = 6, dpi = 300)
print(region_summary)
print(variance_summary)

print(p_gaussian_base)
print(p_binomial_base)

print(adv_gaus$cat)
print(adv_gaus$comb)

print(adv_bios$cat)
print(adv_bios$comb)

plot_gaus <- ggplot(inco_est_long) +
  geom_sf(aes(fill = estimate_plot), colour = NA) +
  facet_wrap(~source, nrow = 1) +
  scale_fill_gradientn(
    colours = fill_pal,
    name = "Transformed income",
    na.value = "grey90"
  ) +
  labs(
    title = "Transformed Income by PUMA",
    subtitle = "INCO = min–max scaled log income"
  ) +
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
  scale_fill_gradientn(
    colours = fill_pal,
    name = "Poverty rate",
    na.value = "grey90"
  ) +
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

top_row <- (plot_gaus | p_sig_gaus) +
  patchwork::plot_layout(widths = c(3, 1))

bottom_row <- (plot_bios | p_sig_bios) +
  patchwork::plot_layout(widths = c(3, 1))

final_region_plot <- top_row / bottom_row +
  patchwork::plot_layout(heights = c(1, 1))

ggsave(
  file.path("figs", "region_basis_main.png"),
  final_region_plot,
  width = 13,
  height = 8.5,
  dpi = 300
)

print(final_region_plot)

sum(region_summary$uni_binomial_var / region_summary$multi_binomial_var >1)/88
sum(region_summary$uni_gaussian_var / region_summary$multi_gaussian_var >1)/88
