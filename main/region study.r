source("./main/packages.r")
source("./main/functions.r")

## Run settings
nburn <- 10000L
nsim <- 30000L
nthin <- 1L
seed <- 202309L
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
output_dir <- "./main/results/region"

## ACS 2023 Illinois person records: adult, positive-income analysis sample
variables <- c("PUMA", "PWGTP", "SEX", "POVPIP", "PINCP", "SCHL", "AGEP")
pums <- readr::read_csv("./main/data/IL23.csv", show_col_types = FALSE) |>
  dplyr::select(dplyr::all_of(variables)) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    PUMA = sprintf("%05d", as.integer(PUMA)),
    PWGTP = as.numeric(PWGTP), PINCP_NUM = as.numeric(PINCP),
    POVPIP = as.numeric(POVPIP), SEX = factor(SEX),
    POV = dplyr::if_else(POVPIP <= 100, 1, 0),
    BACH = factor(dplyr::if_else(SCHL >= 21, 1, 0))
  ) |>
  dplyr::filter(is.finite(PINCP_NUM), PINCP_NUM > 0) |>
  dplyr::arrange(PUMA, SEX, BACH) |>
  dplyr::filter(as.numeric(AGEP) >= 18)
areas <- unique(pums$PUMA)

## Unweighted region Box-Cox profile, conditional on sex and education
# The region lambda uses the complete analysis sample, before the weighted model fit.
lambda_grid <- seq(-2, 2, by = 0.01)
lambda_design <- model.matrix(~ ., as.data.frame(pums[c("SEX", "BACH")]))
lambda_profile <- rep(-Inf, length(lambda_grid))
log_income_sum <- sum(log(pums$PINCP_NUM))
for (g in seq_along(lambda_grid)) {
  transformed <- boxcox_transform(pums$PINCP_NUM, lambda_grid[g])
  working_fit <- lm.fit(lambda_design, transformed)
  rss <- sum(working_fit$residuals^2)
  if (is.finite(rss) && rss > 0) {
    lambda_profile[g] <- -nrow(pums) / 2 * log(rss / nrow(pums)) +
      (lambda_grid[g] - 1) * log_income_sum
  }
}
lambda <- lambda_grid[which.max(lambda_profile)]
pums$INCO <- boxcox_transform(pums$PINCP_NUM, lambda)

## 2020 cartographic PUMA boundaries and the adjacency basis
puma_sf <- readRDS("./main/data/il_puma_2020_cb_sf.rds")
puma_sf$PUMA <- sprintf("%05d", as.integer(puma_sf$PUMACE20))
puma_sf <- puma_sf[match(areas, puma_sf$PUMA), , drop = FALSE]
stopifnot(identical(puma_sf$PUMA, areas))
neighbors <- spdep::poly2nb(puma_sf, queen = TRUE)
adjacency <- spdep::nb2mat(neighbors, style = "B", zero.policy = TRUE)
adjacency <- (adjacency + t(adjacency)) / 2
eig <- eigen(adjacency, symmetric = TRUE)
# Region retains eigenvalues >1e-8; the empirical study retains eigenvalues >0.
basis <- eig$vectors[, eig$values > 1e-8, drop = FALSE]
rownames(basis) <- areas

## ACS 2023 one-year B15001: sex by education for the population aged 18+
# Sex totals are B15001_002/_043; degree cells include bachelor's and graduate degrees.
maps <- readRDS("./main/data/acs2023_b15001_var_maps.rds")
raw_cells <- readRDS("./main/data/acs2023_il_puma_b15001_raw_cells.rds")
sex_totals <- raw_cells |>
  dplyr::inner_join(maps$sex_total_vars, by = c("variable" = "var")) |>
  dplyr::group_by(GEOID, NAME, SEX) |>
  dplyr::summarise(total = sum(estimate), .groups = "drop")
bachelor_totals <- raw_cells |>
  dplyr::inner_join(maps$bach_vars, by = c("variable" = "var")) |>
  dplyr::group_by(GEOID, NAME, SEX) |>
  dplyr::summarise(bach = sum(estimate), .groups = "drop")
cells <- dplyr::left_join(sex_totals, bachelor_totals, by = c("GEOID", "NAME", "SEX")) |>
  dplyr::mutate(bach = dplyr::coalesce(bach, 0))
stopifnot(all(is.finite(cells$total) & is.finite(cells$bach)),
  all(cells$total >= cells$bach & cells$bach >= 0))
pcells <- cells |>
  dplyr::transmute(
    PUMA = sprintf("%05d", as.integer(substr(GEOID, nchar(GEOID) - 4, nchar(GEOID)))),
    SEX = factor(dplyr::if_else(SEX == "Male", "1", "2"), levels = levels(pums$SEX)),
    bach_yes = bach, bach_no = total - bach
  ) |>
  tidyr::pivot_longer(c(bach_yes, bach_no), names_to = "BACH", values_to = "popsize") |>
  dplyr::mutate(BACH = factor(dplyr::if_else(BACH == "bach_yes", "1", "0"),
    levels = levels(pums$BACH))) |>
  dplyr::filter(PUMA %in% areas) |>
  dplyr::arrange(PUMA, SEX, BACH)
stopifnot(nrow(pcells) == length(areas) * 4L, !anyNA(pcells),
  all(pcells$popsize >= 0 & pcells$popsize == floor(pcells$popsize)))

## Model and prediction matrices; survey weights sum to the sample size
X <- model.matrix(~ SEX + BACH - 1, pums)
S <- basis[match(pums$PUMA, rownames(basis)), , drop = FALSE]
predX <- model.matrix(~ SEX + BACH - 1, pcells)
predS <- basis[match(pcells$PUMA, rownames(basis)), , drop = FALSE]
wgt <- pums$PWGTP * nrow(pums) / sum(pums$PWGTP)
cell_area <- match(pcells$PUMA, areas)
cells_by_area <- split(seq_len(nrow(pcells)), as.character(pcells$PUMA))
# Population size: Sum the cell counts within each PUMA.
population_by_area <- sapply(cells_by_area[areas], function(ids) sum(pcells$popsize[ids]))
area_weights <- Matrix::sparseMatrix(
  i = cell_area, j = seq_len(nrow(pcells)),
  x = 1 / population_by_area[cell_area], dims = c(length(areas), nrow(pcells))
)
rownames(area_weights) <- areas
draw_ids <- seq_len(nsim / nthin)
message(nrow(pums), " adults; ", length(areas), " PUMAs; ", ncol(basis),
  " basis columns; lambda = ", lambda, "; population = ", sum(pcells$popsize))
set.seed(seed)

## Gaussian model and all-unit income prediction
fit_ug <- uni_gaus(
  X = X, Y = pums$INCO, S = S, sig2b = 2000, wgt = wgt, n = NULL,
  predX = predX, predS = predS, nburn = nburn, nsim = nsim, nthin = nthin,
  a = 0.5, b = 100, a_eps = 0.1, b_eps = 20, on_draw = NULL
)
income_ug <- pure_r_region_medians(fit_ug$Preds, as.numeric(fit_ug$sig2.chain),
  lambda, pcells$popsize, cell_area)
dimnames(income_ug$medians) <- list(areas, as.character(draw_ids))

## Binomial model and finite-population poverty prediction
fit_ub <- uni_bio(
  X = X, Y = pums$POV, S = S, sig2b = 2000, wgt = wgt, n = NULL,
  predX = predX, predS = predS, nburn = nburn, nsim = nsim, nthin = nthin,
  a = 0.1, b = 0.1, on_draw = NULL
)
poverty_ub <- predict_poverty_region(fit_ub$Preds, pcells$popsize, area_weights)

## Multi-type model: retain paired Gaussian/binomial posterior columns
fit_mt <- multi(
  X_1 = X, X_2 = X, Z_1 = pums$INCO, Z_2 = pums$POV, S = S,
  sig2b = 2000, wgt = wgt, n = NULL, predX = predX, predS = predS, n_preds = NULL,
  nburn = nburn, nsim = nsim, nthin = nthin,
  sig2t = 2000, sig2e = 10, tau_2_init = 0,
  a_eps = 0.1, b_eps = 20, aeta = 0.1, beta = 0.1, alambda = 2, blambda = 1,
  on_draw = NULL
)
income_mg <- pure_r_region_medians(fit_mt$preds_gaus.chain, as.numeric(fit_mt$sig2.chain),
  lambda, pcells$popsize, cell_area)
dimnames(income_mg$medians) <- list(areas, as.character(draw_ids))
poverty_mb <- predict_poverty_region(fit_mt$preds_bios.chain, pcells$popsize, area_weights)
area_draws <- list(ug = income_ug$medians, mg = income_mg$medians,
  ub = poverty_ub, mb = poverty_mb)

## Design-based Direct estimators
pums$W <- pums$PWGTP
direct <- md_direct(pums, areas)
direct_diagnostics <- attr(direct, "diagnostics")
direct$interval_type <- "Approximate design-based Woodruff 95% interval"
design <- survey::svydesign(ids = ~1, weights = ~W, data = pums)
direct_pov <- survey::svyby(~POV, ~PUMA, design, survey::svymean,
  vartype = "se", keep.names = FALSE)
ix <- match(areas, as.character(direct_pov$PUMA))
estimate <- direct_pov$POV[ix]
se <- as.numeric(survey::SE(direct_pov)[ix])
direct_poverty <- data.frame(PUMA = areas, response = "poverty", model = "Direct",
  estimate = estimate, lower = pmax(0, estimate - qnorm(.975) * se),
  upper = pmin(1, estimate + qnorm(.975) * se), variance = se^2,
  interval_type = "Approximate design-based normal 95% interval")

## Posterior means, equal-tailed intervals and variances
ug_summary <- summarize_draws(area_draws$ug)
ug_summary$response <- "income"
ug_summary$model <- "Univariate"
mg_summary <- summarize_draws(area_draws$mg)
mg_summary$response <- "income"
mg_summary$model <- "Multi-type"
ub_summary <- summarize_draws(area_draws$ub)
ub_summary$response <- "poverty"
ub_summary$model <- "Univariate"
mb_summary <- summarize_draws(area_draws$mb)
mb_summary$response <- "poverty"
mb_summary$model <- "Multi-type"
ug_summary$interval_type <- mg_summary$interval_type <- ub_summary$interval_type <-
  mb_summary$interval_type <- "Posterior equal-tailed 95% interval"
estimates <- rbind(direct, direct_poverty, ug_summary, mg_summary, ub_summary, mb_summary)
estimates$n_sample <- as.integer(table(factor(pums$PUMA, areas)))[match(estimates$PUMA, areas)]
estimates$population <- as.numeric(rowsum(pcells$popsize, cell_area))[match(estimates$PUMA, areas)]
estimates <- estimates[order(match(estimates$response, c("income", "poverty")),
  match(estimates$model, c("Direct", "Univariate", "Multi-type")), estimates$PUMA), ]
comparisons <- vector("list", 2L)
variance_summaries <- vector("list", 2L)
responses <- c("income", "poverty")
for (i in seq_along(responses)) {
  response <- responses[i]
  u <- estimates[estimates$response == response & estimates$model == "Univariate", ]
  m <- estimates[estimates$response == response & estimates$model == "Multi-type", ]
  stopifnot(identical(u$PUMA, m$PUMA))
  comparisons[[i]] <- data.frame(PUMA = u$PUMA, response = response,
    Univariate = u$variance, `Multi-type` = m$variance,
    multi_variance_lower = m$variance < u$variance,
    multi_to_uni_ratio = m$variance / u$variance, check.names = FALSE)
  variance_summaries[[i]] <- data.frame(response = response, n_PUMA = nrow(u),
    n_multi_lower = sum(m$variance < u$variance),
    proportion_multi_lower = mean(m$variance < u$variance),
    mean_multi_variance = mean(m$variance), mean_uni_variance = mean(u$variance),
    ratio_of_mean_variances = mean(m$variance) / mean(u$variance))
}
variance_comparison <- do.call(rbind, comparisons)
variance_summary <- do.call(rbind, variance_summaries)

## Export estimates and Figure 2 data
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(estimates, file.path(table_dir, "region_area_estimates.csv"), row.names = FALSE)
write.csv(variance_comparison, file.path(table_dir, "region_variance_comparison.csv"), row.names = FALSE)
write.csv(variance_summary, file.path(table_dir, "region_variance_summary.csv"), row.names = FALSE)
saveRDS(list(ug = fit_ug, ub = fit_ub, mt = fit_mt), file.path(output_dir, "model_fits.rds"))
saveRDS(list(area_draws = area_draws,
  audit = list(ug = income_ug$audit, mg = income_mg$audit), lambda = lambda,
  settings = list(nburn = nburn, nsim = nsim, nthin = nthin, seed = seed)),
  file.path(output_dir, "area_draws.rds"))
saveRDS(list(estimates = estimates, variance_comparison = variance_comparison,
  variance_summary = variance_summary, direct_diagnostics = direct_diagnostics),
  file.path(output_dir, "region_results.rds"))

## Figure 2: Direct / Multi-type / Univariate maps and variance comparisons
geometry <- sf::st_as_sf(puma_sf)
map_data <- estimates |>
  dplyr::select(PUMA, response, model, estimate) |>
  dplyr::mutate(model = factor(model, c("Direct", "Multi-type", "Univariate")))
fill_pal <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
# map_plot: Map one response for the three estimators by PUMA.
map_plot <- function(response_name) {
  d <- geometry |> dplyr::select(PUMA) |>
    dplyr::left_join(dplyr::filter(map_data, response == response_name), by = "PUMA")
  ggplot(d) + geom_sf(aes(fill = estimate), color = "white", linewidth = .07) +
    facet_wrap(~model, nrow = 1) + coord_sf(datum = NA) +
    scale_fill_gradientn(colours = fill_pal,
      name = if (response_name == "income") "Income ($)" else "Poverty rate",
      labels = if (response_name == "income") scales::label_comma() else scales::label_number(),
      guide = guide_colorbar(title.position = "top", barheight = grid::unit(1.3, "cm"))) +
    labs(title = if (response_name == "income") "Median Income by PUMA" else "Poverty Rate by PUMA") +
    theme_void(base_size = 11) + theme(
      plot.title = element_text(face = "bold", size = 15),
      strip.text = element_text(size = 10), legend.title = element_text(size = 10),
      legend.text = element_text(size = 9), plot.margin = margin(3, 2, 3, 0))
}
# variance_plot: Compare multi-type and univariate posterior variances.
variance_plot <- function(response_name) {
  mult <- if (response_name == "income") 1e-6 else 1e6
  d <- variance_comparison |> dplyr::filter(response == response_name) |>
    dplyr::mutate(multi = .data[["Multi-type"]] * mult, uni = Univariate * mult,
      comparison = ifelse(multi < uni, "Smaller", "Larger"))
  limits <- c(0, max(d$multi, d$uni) * 1.05)
  units <- if (response_name == "income") " (10^6 dollars^2)" else " (x 1,000,000)"
  ggplot(d, aes(multi, uni, color = comparison)) +
    geom_abline(slope = 1, intercept = 0, color = "grey40", linetype = "dashed") +
    geom_point(size = 1.6, alpha = .85) +
    scale_color_manual(values = c(Smaller = "#2166AC", Larger = "#B2182B"), guide = "none") +
    coord_equal(xlim = limits, ylim = limits) +
    labs(title = paste(if (response_name == "income") "Income" else "Bernoulli", "Variance Comparison"),
      x = paste0("Multi-type", units), y = paste0("Univariate", units)) +
    theme_minimal(base_size = 11) + theme(
      plot.title = element_text(face = "bold", size = 13),
      axis.title = element_text(size = 9), panel.grid.minor = element_blank())
}
top <- patchwork::wrap_plots(map_plot("income"), variance_plot("income"),
  nrow = 1, widths = c(2.2, 1.4))
bottom <- patchwork::wrap_plots(map_plot("poverty"), variance_plot("poverty"),
  nrow = 1, widths = c(2.2, 1.4))
figure <- patchwork::wrap_plots(top, bottom, ncol = 1)
stem <- file.path(figure_dir, "region_median_main")
ggsave(paste0(stem, ".png"), figure, device = grDevices::png,
  width = 13.5, height = 8.5, dpi = 300, bg = "white")
ggsave(paste0(stem, ".pdf"), figure, device = grDevices::pdf,
  width = 13.5, height = 8.5, bg = "white")
if (interactive()) print(figure)
print(variance_summary)
