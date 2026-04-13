source("basis version/packages.r")
source("basis version/functions.r")
interval_score <- function(lower, upper, truth, alpha = 0.05) {
  width <- upper - lower
  penalty_low  <- (2 / alpha) * pmax(lower - truth, 0)
  penalty_high <- (2 / alpha) * pmax(truth - upper, 0)
  width + penalty_low + penalty_high
}

load(file.path("data", "empirical_basis_tau_in_gauss.RData"))
## traceplot of fit
n_area <- nrow(truth)
n_rep <- dim(ubios_pre)[2]

cr_ub <- cr_mb <- cr_ug <- cr_mg <- numeric(n_area)
is_ub <- is_mb <- is_ug <- is_mg <- numeric(n_area)
mse_ub <- mse_mb <- mse_ug <- mse_mg <- mse_db <- mse_dg <- numeric(n_area)

for (j in seq_len(n_area)) {
  cr_ub[j] <- mean(
    ubios_qual[j, 1, ] <= truth$POV[j] & truth$POV[j] <= ubios_qual[j, 2, ],
    na.rm = TRUE
  )
  cr_mb[j] <- mean(
    mbios_qual[j, 1, ] <= truth$POV[j] & truth$POV[j] <= mbios_qual[j, 2, ],
    na.rm = TRUE
  )
  cr_ug[j] <- mean(
    ugaus_qual[j, 1, ] <= truth$INCO[j] & truth$INCO[j] <= ugaus_qual[j, 2, ],
    na.rm = TRUE
  )
  cr_mg[j] <- mean(
    mgaus_qual[j, 1, ] <= truth$INCO[j] & truth$INCO[j] <= mgaus_qual[j, 2, ],
    na.rm = TRUE
  )

  is_ub[j] <- mean(interval_score(ubios_qual[j, 1, ], ubios_qual[j, 2, ], truth$POV[j]), na.rm = TRUE)
  is_mb[j] <- mean(interval_score(mbios_qual[j, 1, ], mbios_qual[j, 2, ], truth$POV[j]), na.rm = TRUE)
  is_ug[j] <- mean(interval_score(ugaus_qual[j, 1, ], ugaus_qual[j, 2, ], truth$INCO[j]), na.rm = TRUE)
  is_mg[j] <- mean(interval_score(mgaus_qual[j, 1, ], mgaus_qual[j, 2, ], truth$INCO[j]), na.rm = TRUE)

  mse_ub[j] <- mean((ubios_pre[j, ] - truth$POV[j])^2, na.rm = TRUE)
  mse_mb[j] <- mean((mbios_pre[j, ] - truth$POV[j])^2, na.rm = TRUE)
  mse_ug[j] <- mean((ugaus_pre[j, ] - truth$INCO[j])^2, na.rm = TRUE)
  mse_mg[j] <- mean((mgaus_pre[j, ] - truth$INCO[j])^2, na.rm = TRUE)
  mse_db[j] <- mean((dbio_pre[j, ] - truth$POV[j])^2, na.rm = TRUE)
  mse_dg[j] <- mean((dgaus_pre[j, ] - truth$INCO[j])^2, na.rm = TRUE)
}

summary_table_area <- tibble(
  `Response / Model` = c(
    "Bernoulli, HT",
    "Bernoulli, Univariate",
    "Bernoulli, Multi-type",
    "Gaussian, HT",
    "Gaussian, Univariate",
    "Gaussian, Multi-type"
  ),
  MSE = c(
    mean(mse_db, na.rm = TRUE),
    mean(mse_ub, na.rm = TRUE),
    mean(mse_mb, na.rm = TRUE),
    mean(mse_dg, na.rm = TRUE),
    mean(mse_ug, na.rm = TRUE),
    mean(mse_mg, na.rm = TRUE)
  ),
  Ratio_to_direct = c(
    1,
    mean(mse_ub / mse_db, na.rm = TRUE),
    mean(mse_mb / mse_db, na.rm = TRUE),
    1,
    mean(mse_ug / mse_dg, na.rm = TRUE),
    mean(mse_mg / mse_dg, na.rm = TRUE)
  ),
  Mean_IS = c(
    NA_real_,
    mean(is_ub, na.rm = TRUE),
    mean(is_mb, na.rm = TRUE),
    NA_real_,
    mean(is_ug, na.rm = TRUE),
    mean(is_mg, na.rm = TRUE)
  ),
  Mean_CR = c(
    NA_real_,
    mean(cr_ub, na.rm = TRUE),
    mean(cr_mb, na.rm = TRUE),
    NA_real_,
    mean(cr_ug, na.rm = TRUE),
    mean(cr_mg, na.rm = TRUE)
  ),
  Coverage_Gap = c(
    NA_real_,
    mean(abs(cr_ub - 0.95), na.rm = TRUE),
    mean(abs(cr_mb - 0.95), na.rm = TRUE),
    NA_real_,
    mean(abs(cr_ug - 0.95), na.rm = TRUE),
    mean(abs(cr_mg - 0.95), na.rm = TRUE)
  )
)

mse_db_rep <- colMeans((dbio_pre - truth$POV)^2, na.rm = TRUE)
mse_ub_rep <- colMeans((ubios_pre - truth$POV)^2, na.rm = TRUE)
mse_mb_rep <- colMeans((mbios_pre - truth$POV)^2, na.rm = TRUE)

mse_dg_rep <- colMeans((dgaus_pre - truth$INCO)^2, na.rm = TRUE)
mse_ug_rep <- colMeans((ugaus_pre - truth$INCO)^2, na.rm = TRUE)
mse_mg_rep <- colMeans((mgaus_pre - truth$INCO)^2, na.rm = TRUE)

is_ub_rep <- sapply(seq_len(n_rep), function(k) {
  mean(interval_score(ubios_qual[, 1, k], ubios_qual[, 2, k], truth$POV), na.rm = TRUE)
})
is_mb_rep <- sapply(seq_len(n_rep), function(k) {
  mean(interval_score(mbios_qual[, 1, k], mbios_qual[, 2, k], truth$POV), na.rm = TRUE)
})
is_ug_rep <- sapply(seq_len(n_rep), function(k) {
  mean(interval_score(ugaus_qual[, 1, k], ugaus_qual[, 2, k], truth$INCO), na.rm = TRUE)
})
is_mg_rep <- sapply(seq_len(n_rep), function(k) {
  mean(interval_score(mgaus_qual[, 1, k], mgaus_qual[, 2, k], truth$INCO), na.rm = TRUE)
})

cr_ub_rep <- sapply(seq_len(n_rep), function(k) {
  mean(ubios_qual[, 1, k] <= truth$POV & truth$POV <= ubios_qual[, 2, k], na.rm = TRUE)
})
cr_mb_rep <- sapply(seq_len(n_rep), function(k) {
  mean(mbios_qual[, 1, k] <= truth$POV & truth$POV <= mbios_qual[, 2, k], na.rm = TRUE)
})
cr_ug_rep <- sapply(seq_len(n_rep), function(k) {
  mean(ugaus_qual[, 1, k] <= truth$INCO & truth$INCO <= ugaus_qual[, 2, k], na.rm = TRUE)
})
cr_mg_rep <- sapply(seq_len(n_rep), function(k) {
  mean(mgaus_qual[, 1, k] <= truth$INCO & truth$INCO <= mgaus_qual[, 2, k], na.rm = TRUE)
})

summary_table_rep <- tibble(
  `Response / Model` = c(
    "Bernoulli, HT",
    "Bernoulli, Univariate",
    "Bernoulli, Multi-type",
    "Gaussian, HT",
    "Gaussian, Univariate",
    "Gaussian, Multi-type"
  ),
  MSE = c(
    mean(mse_db_rep, na.rm = TRUE),
    mean(mse_ub_rep, na.rm = TRUE),
    mean(mse_mb_rep, na.rm = TRUE),
    mean(mse_dg_rep, na.rm = TRUE),
    mean(mse_ug_rep, na.rm = TRUE),
    mean(mse_mg_rep, na.rm = TRUE)
  ),
  Ratio_to_direct = c(
    1,
    mean(mse_ub_rep / mse_db_rep, na.rm = TRUE),
    mean(mse_mb_rep / mse_db_rep, na.rm = TRUE),
    1,
    mean(mse_ug_rep / mse_dg_rep, na.rm = TRUE),
    mean(mse_mg_rep / mse_dg_rep, na.rm = TRUE)
  ),
  Mean_IS = c(
    NA_real_,
    mean(is_ub_rep, na.rm = TRUE),
    mean(is_mb_rep, na.rm = TRUE),
    NA_real_,
    mean(is_ug_rep, na.rm = TRUE),
    mean(is_mg_rep, na.rm = TRUE)
  ),
  Mean_CR = c(
    NA_real_,
    mean(cr_ub_rep, na.rm = TRUE),
    mean(cr_mb_rep, na.rm = TRUE),
    NA_real_,
    mean(cr_ug_rep, na.rm = TRUE),
    mean(cr_mg_rep, na.rm = TRUE)
  ),
  Coverage_Gap = c(
    NA_real_,
    mean(abs(cr_ub_rep - 0.95), na.rm = TRUE),
    mean(abs(cr_mb_rep - 0.95), na.rm = TRUE),
    NA_real_,
    mean(abs(cr_ug_rep - 0.95), na.rm = TRUE),
    mean(abs(cr_mg_rep - 0.95), na.rm = TRUE)
  )
)

paired_rep <- tibble(
  rep = seq_len(n_rep),
  mse_db = mse_db_rep,
  mse_ub = mse_ub_rep,
  mse_mb = mse_mb_rep,
  mse_dg = mse_dg_rep,
  mse_ug = mse_ug_rep,
  mse_mg = mse_mg_rep,
  is_ub = is_ub_rep,
  is_mb = is_mb_rep,
  is_ug = is_ug_rep,
  is_mg = is_mg_rep,
  cr_ub = cr_ub_rep,
  cr_mb = cr_mb_rep,
  cr_ug = cr_ug_rep,
  cr_mg = cr_mg_rep
) |>
  mutate(
    mse_ratio_mb_ub = mse_mb / mse_ub,
    mse_ratio_mg_ug = mse_mg / mse_ug,
    mse_ratio_ub_db = mse_ub / mse_db,
    mse_ratio_ug_dg = mse_ug / mse_dg,
    is_ratio_mb_ub = is_mb / is_ub,
    is_ratio_mg_ug = is_mg / is_ug,
    cr_gap_ub = abs(cr_ub - 0.95),
    cr_gap_mb = abs(cr_mb - 0.95),
    cr_gap_ug = abs(cr_ug - 0.95),
    cr_gap_mg = abs(cr_mg - 0.95),
    multi_better_mse_bin = mse_ratio_mb_ub < 1,
    multi_better_mse_gau = mse_ratio_mg_ug < 1,
    multi_better_is_bin = is_ratio_mb_ub < 1,
    multi_better_is_gau = is_ratio_mg_ug < 1,
    multi_better_cr_bin = cr_gap_mb < cr_gap_ub,
    multi_better_cr_gau = cr_gap_mg < cr_gap_ug
  )

basis_absolute_table <- tibble(
  Response = c("Bernoulli", "Bernoulli", "Bernoulli",
               "Gaussian", "Gaussian", "Gaussian"),
  Model = c("HT", "Univariate", "Multi-type",
            "HT", "Univariate", "Multi-type"),
  MSE = c(
    mean(mse_db_rep, na.rm = TRUE),
    mean(mse_ub_rep, na.rm = TRUE),
    mean(mse_mb_rep, na.rm = TRUE),
    mean(mse_dg_rep, na.rm = TRUE),
    mean(mse_ug_rep, na.rm = TRUE),
    mean(mse_mg_rep, na.rm = TRUE)
  ),
  Ratio_to_HT = c(
    1,
    mean(mse_ub_rep / mse_db_rep, na.rm = TRUE),
    mean(mse_mb_rep / mse_db_rep, na.rm = TRUE),
    1,
    mean(mse_ug_rep / mse_dg_rep, na.rm = TRUE),
    mean(mse_mg_rep / mse_dg_rep, na.rm = TRUE)
  ),
  Interval_Score = c(
    NA_real_,
    mean(is_ub_rep, na.rm = TRUE),
    mean(is_mb_rep, na.rm = TRUE),
    NA_real_,
    mean(is_ug_rep, na.rm = TRUE),
    mean(is_mg_rep, na.rm = TRUE)
  ),
  Coverage_Rate = c(
    NA_real_,
    mean(cr_ub_rep, na.rm = TRUE),
    mean(cr_mb_rep, na.rm = TRUE),
    NA_real_,
    mean(cr_ug_rep, na.rm = TRUE),
    mean(cr_mg_rep, na.rm = TRUE)
  ),
  Coverage_Gap = c(
    NA_real_,
    mean(abs(cr_ub_rep - 0.95), na.rm = TRUE),
    mean(abs(cr_mb_rep - 0.95), na.rm = TRUE),
    NA_real_,
    mean(abs(cr_ug_rep - 0.95), na.rm = TRUE),
    mean(abs(cr_mg_rep - 0.95), na.rm = TRUE)
  )
)

basis_paired_table <- tibble(
  Response = c("Bernoulli", "Gaussian"),
  Median_Multi_Uni_MSE_Ratio = c(
    median(paired_rep$mse_ratio_mb_ub, na.rm = TRUE),
    median(paired_rep$mse_ratio_mg_ug, na.rm = TRUE)
  ),
  Mean_Multi_Uni_MSE_Ratio = c(
    mean(paired_rep$mse_ratio_mb_ub, na.rm = TRUE),
    mean(paired_rep$mse_ratio_mg_ug, na.rm = TRUE)
  ),
  Prop_MSE_Multi_Better = c(
    mean(paired_rep$multi_better_mse_bin, na.rm = TRUE),
    mean(paired_rep$multi_better_mse_gau, na.rm = TRUE)
  ),
  Median_Multi_Uni_IS_Ratio = c(
    median(paired_rep$is_ratio_mb_ub, na.rm = TRUE),
    median(paired_rep$is_ratio_mg_ug, na.rm = TRUE)
  ),
  Mean_Multi_Uni_IS_Ratio = c(
    mean(paired_rep$is_ratio_mb_ub, na.rm = TRUE),
    mean(paired_rep$is_ratio_mg_ug, na.rm = TRUE)
  ),
  Prop_IS_Multi_Better = c(
    mean(paired_rep$multi_better_is_bin, na.rm = TRUE),
    mean(paired_rep$multi_better_is_gau, na.rm = TRUE)
  ),
  Prop_CoverageGap_Multi_Better = c(
    mean(paired_rep$multi_better_cr_bin, na.rm = TRUE),
    mean(paired_rep$multi_better_cr_gau, na.rm = TRUE)
  )
)

better_prop_table <- tibble(
  Response = c("Bernoulli", "Gaussian"),
  prop_MSE_multi_better = c(
    mean(paired_rep$multi_better_mse_bin, na.rm = TRUE),
    mean(paired_rep$multi_better_mse_gau, na.rm = TRUE)
  ),
  prop_IS_multi_better = c(
    mean(paired_rep$multi_better_is_bin, na.rm = TRUE),
    mean(paired_rep$multi_better_is_gau, na.rm = TRUE)
  ),
  prop_CR_multi_better = c(
    mean(paired_rep$multi_better_cr_bin, na.rm = TRUE),
    mean(paired_rep$multi_better_cr_gau, na.rm = TRUE)
  )
)

ratio_plot_df <- bind_rows(
  tibble(rep = seq_len(n_rep), Response = "Bernoulli", Metric = "MSE", Ratio = paired_rep$mse_ratio_mb_ub),
  tibble(rep = seq_len(n_rep), Response = "Gaussian",  Metric = "MSE", Ratio = paired_rep$mse_ratio_mg_ug),
  tibble(rep = seq_len(n_rep), Response = "Bernoulli", Metric = "IS",  Ratio = paired_rep$is_ratio_mb_ub),
  tibble(rep = seq_len(n_rep), Response = "Gaussian",  Metric = "IS",  Ratio = paired_rep$is_ratio_mg_ug)
)

coverage_plot_df <- bind_rows(
  tibble(rep = seq_len(n_rep), Response = "Bernoulli", Model = "Univariate", Gap = paired_rep$cr_gap_ub),
  tibble(rep = seq_len(n_rep), Response = "Bernoulli", Model = "Multi-type", Gap = paired_rep$cr_gap_mb),
  tibble(rep = seq_len(n_rep), Response = "Gaussian",  Model = "Univariate", Gap = paired_rep$cr_gap_ug),
  tibble(rep = seq_len(n_rep), Response = "Gaussian",  Model = "Multi-type", Gap = paired_rep$cr_gap_mg)
)

better_plot_df <- better_prop_table |>
  pivot_longer(
    cols = -Response,
    names_to = "Metric",
    values_to = "Proportion"
  ) |>
  mutate(
    Metric = recode(
      Metric,
      prop_MSE_multi_better = "MSE: multi < uni",
      prop_IS_multi_better = "IS: multi < uni",
      prop_CR_multi_better = "Coverage gap: multi better"
    )
  )

p_ratio <- ggplot(ratio_plot_df, aes(x = Response, y = Ratio, fill = Response)) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 0.9) +
  geom_violin(trim = FALSE, alpha = 0.55, colour = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.08, height = 0, size = 1.1, alpha = 0.45, colour = "grey35") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Bernoulli" = "#F4A261")) +
  labs(x = NULL, y = "Multi / Univariate ratio") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

p_covgap <- ggplot(coverage_plot_df, aes(x = Model, y = Gap, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.55, colour = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.08, height = 0, size = 1.1, alpha = 0.45, colour = "grey35") +
  facet_wrap(~ Response, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Univariate" = "#BDBDBD", "Multi-type" = "#4C78A8")) +
  labs(x = NULL, y = "|Coverage - 0.95|") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

p_better <- ggplot(better_plot_df, aes(x = Response, y = Proportion, fill = Response)) +
  geom_col(width = 0.65) +
  facet_wrap(~ Metric, nrow = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Bernoulli" = "#F4A261")) +
  labs(x = NULL, y = "Proportion of replicates") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

p_direct_compare <- bind_rows(
  tibble(rep = seq_len(n_rep), Response = "Bernoulli", Comparison = "Univariate / HT", Ratio = paired_rep$mse_ratio_ub_db),
  tibble(rep = seq_len(n_rep), Response = "Bernoulli", Comparison = "Multi-type / HT", Ratio = paired_rep$mse_mb / paired_rep$mse_db),
  tibble(rep = seq_len(n_rep), Response = "Gaussian",  Comparison = "Univariate / HT", Ratio = paired_rep$mse_ratio_ug_dg),
  tibble(rep = seq_len(n_rep), Response = "Gaussian",  Comparison = "Multi-type / HT", Ratio = paired_rep$mse_mg / paired_rep$mse_dg)
) |>
  ggplot(aes(x = Comparison, y = Ratio, fill = Response)) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 0.9) +
  geom_violin(trim = FALSE, alpha = 0.55, colour = NA) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    alpha = 0.9,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
    size = 1.1,
    alpha = 0.45,
    colour = "grey35"
  ) +
  facet_wrap(~ Response, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Bernoulli" = "#F4A261")) +
  labs(x = NULL, y = "Ratio to HT") +
  theme_bw(base_size = 14)

final_plot <- (p_ratio / p_covgap) / p_better +
  patchwork::plot_layout(heights = c(1, 1, 0.9))

final_plot2 <- p_direct_compare

print(summary_table_area)
print(summary_table_rep)
print(basis_absolute_table)
print(basis_paired_table)
print(better_prop_table)
print(final_plot)
print(final_plot2)

readr::write_csv(summary_table_area, "basis_summary_table_area.csv")
readr::write_csv(summary_table_rep, "basis_summary_table_rep.csv")
readr::write_csv(basis_absolute_table, "basis_absolute_table.csv")
readr::write_csv(basis_paired_table, "basis_paired_table.csv")
readr::write_csv(better_prop_table, "basis_better_prop_table.csv")
readr::write_csv(paired_rep, "basis_paired_replicate_metrics.csv")

ggsave(file.path("figs", "empirical_basis_main_plot.png"), final_plot, width = 14, height = 12, dpi = 300)
ggsave(file.path("figs", "empirical_basis_direct_compare_plot.png"), final_plot2, width = 12, height = 5, dpi = 300)


area_id <- if ("PUMA" %in% names(truth)) as.character(truth$PUMA) else sprintf("Area_%03d", seq_len(n_area))

paired_area <- tibble(
  area = area_id,
  mse_db = mse_db,
  mse_ub = mse_ub,
  mse_mb = mse_mb,
  mse_dg = mse_dg,
  mse_ug = mse_ug,
  mse_mg = mse_mg,
  is_ub = is_ub,
  is_mb = is_mb,
  is_ug = is_ug,
  is_mg = is_mg,
  cr_ub = cr_ub,
  cr_mb = cr_mb,
  cr_ug = cr_ug,
  cr_mg = cr_mg
) |>
  mutate(
    mse_ratio_mb_ub = mse_mb / mse_ub,
    mse_ratio_mg_ug = mse_mg / mse_ug,
    mse_ratio_ub_db = mse_ub / mse_db,
    mse_ratio_ug_dg = mse_ug / mse_dg,
    mse_ratio_mb_db = mse_mb / mse_db,
    mse_ratio_mg_dg = mse_mg / mse_dg,
    is_ratio_mb_ub = is_mb / is_ub,
    is_ratio_mg_ug = is_mg / is_ug,
    cr_gap_ub = abs(cr_ub - 0.95),
    cr_gap_mb = abs(cr_mb - 0.95),
    cr_gap_ug = abs(cr_ug - 0.95),
    cr_gap_mg = abs(cr_mg - 0.95),
    multi_better_mse_bin = mse_ratio_mb_ub < 1,
    multi_better_mse_gau = mse_ratio_mg_ug < 1,
    multi_better_is_bin = is_ratio_mb_ub < 1,
    multi_better_is_gau = is_ratio_mg_ug < 1,
    multi_better_cr_bin = cr_gap_mb < cr_gap_ub,
    multi_better_cr_gau = cr_gap_mg < cr_gap_ug
  )

better_prop_area_table <- tibble(
  Response = c("Bernoulli", "Gaussian"),
  prop_MSE_multi_better = c(
    mean(paired_area$multi_better_mse_bin, na.rm = TRUE),
    mean(paired_area$multi_better_mse_gau, na.rm = TRUE)
  ),
  prop_IS_multi_better = c(
    mean(paired_area$multi_better_is_bin, na.rm = TRUE),
    mean(paired_area$multi_better_is_gau, na.rm = TRUE)
  ),
  prop_CR_multi_better = c(
    mean(paired_area$multi_better_cr_bin, na.rm = TRUE),
    mean(paired_area$multi_better_cr_gau, na.rm = TRUE)
  )
)

ratio_plot_area_df <- bind_rows(
  tibble(area = area_id, Response = "Bernoulli", Metric = "MSE", Ratio = paired_area$mse_ratio_mb_ub),
  tibble(area = area_id, Response = "Gaussian",  Metric = "MSE", Ratio = paired_area$mse_ratio_mg_ug),
  tibble(area = area_id, Response = "Bernoulli", Metric = "IS",  Ratio = paired_area$is_ratio_mb_ub),
  tibble(area = area_id, Response = "Gaussian",  Metric = "IS",  Ratio = paired_area$is_ratio_mg_ug)
)

coverage_plot_area_df <- bind_rows(
  tibble(area = area_id, Response = "Bernoulli", Model = "Univariate", Gap = paired_area$cr_gap_ub),
  tibble(area = area_id, Response = "Bernoulli", Model = "Multi-type", Gap = paired_area$cr_gap_mb),
  tibble(area = area_id, Response = "Gaussian",  Model = "Univariate", Gap = paired_area$cr_gap_ug),
  tibble(area = area_id, Response = "Gaussian",  Model = "Multi-type", Gap = paired_area$cr_gap_mg)
)

better_plot_area_df <- better_prop_area_table |>
  pivot_longer(
    cols = -Response,
    names_to = "Metric",
    values_to = "Proportion"
  ) |>
  mutate(
    Metric = recode(
      Metric,
      prop_MSE_multi_better = "MSE: multi < uni",
      prop_IS_multi_better = "IS: multi < uni",
      prop_CR_multi_better = "Coverage gap: multi better"
    )
  )

p_ratio_area <- ggplot(ratio_plot_area_df, aes(x = Response, y = Ratio, fill = Response)) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 0.9) +
  geom_violin(trim = FALSE, alpha = 0.55, colour = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.08, height = 0, size = 1.1, alpha = 0.45, colour = "grey35") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Bernoulli" = "#F4A261")) +
  labs(x = NULL, y = "Multi / Univariate ratio") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

p_covgap_area <- ggplot(coverage_plot_area_df, aes(x = Model, y = Gap, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.55, colour = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.08, height = 0, size = 1.1, alpha = 0.45, colour = "grey35") +
  facet_wrap(~ Response, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Univariate" = "#BDBDBD", "Multi-type" = "#4C78A8")) +
  labs(x = NULL, y = "|Coverage - 0.95|") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

p_better_area <- ggplot(better_plot_area_df, aes(x = Response, y = Proportion, fill = Response)) +
  geom_col(width = 0.65) +
  facet_wrap(~ Metric, nrow = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Bernoulli" = "#F4A261")) +
  labs(x = NULL, y = "Proportion of areas") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

p_direct_compare_area <- bind_rows(
  tibble(area = area_id, Response = "Bernoulli", Comparison = "Univariate / HT", Ratio = paired_area$mse_ratio_ub_db),
  tibble(area = area_id, Response = "Bernoulli", Comparison = "Multi-type / HT", Ratio = paired_area$mse_ratio_mb_db),
  tibble(area = area_id, Response = "Gaussian",  Comparison = "Univariate / HT", Ratio = paired_area$mse_ratio_ug_dg),
  tibble(area = area_id, Response = "Gaussian",  Comparison = "Multi-type / HT", Ratio = paired_area$mse_ratio_mg_dg)
) |>
  ggplot(aes(x = Comparison, y = Ratio, fill = Response)) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 0.9) +
  geom_violin(trim = FALSE, alpha = 0.55, colour = NA) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    alpha = 0.9,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
    size = 1.1,
    alpha = 0.45,
    colour = "grey35"
  ) +
  facet_wrap(~ Response, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Bernoulli" = "#F4A261")) +
  labs(x = NULL, y = "Ratio to HT") +
  theme_bw(base_size = 14)

final_plot_area <- (p_ratio_area / p_covgap_area) / p_better_area +
  patchwork::plot_layout(heights = c(1, 1, 0.9))

final_plot2_area <- p_direct_compare_area

print(paired_area)
print(better_prop_area_table)
print(final_plot_area)
print(final_plot2_area)

readr::write_csv(paired_area, "basis_paired_area_metrics.csv")
readr::write_csv(better_prop_area_table, "basis_better_prop_area_table.csv")

ggsave(file.path("figs", "empirical_basis_main_plot_area.png"), final_plot_area, width = 14, height = 12, dpi = 300)
ggsave(file.path("figs", "empirical_basis_direct_compare_plot_area.png"), final_plot2_area, width = 12, height = 5, dpi = 300)