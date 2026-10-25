source("./main/packages.r")
source("./main/functions.r")

## Run settings
n_rep <- 100L
sample_size <- 1000L
nburn <- 1000L
nsim <- 1000L
nthin <- 1L
chunk_units <- 100000L
lambda_grid <- seq(-2, 2, by = 0.01)
pri <- list(beta_variance = 2000, tau_variance = 2000,
  gaussian_area_shape = 0.5, gaussian_area_rate = 100,
  gaussian_error_shape = 0.1, gaussian_error_rate = 20,
  binary_area_shape = 0.1, binary_area_rate = 0.1,
  shared_shape = 0.1, shared_rate = 0.1,
  binary_specific_shape = 2, binary_specific_rate = 1)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
stopifnot(n_rep >= 1L, sample_size >= 1L, nburn >= 0L,
  nsim >= 2L, nthin >= 1L, nsim %% nthin == 0L)
output_dir <- "./main/results/empirical"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Build the fixed finite population from the original ACS records
pums <- readr::read_csv("./main/data/IL23.csv", show_col_types = FALSE) |>
  dplyr::select(PUMA, PWGTP, SEX, POVPIP, PINCP, SCHL) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    PUMA = sprintf("%05d", as.integer(PUMA)), PWGTP = as.numeric(PWGTP),
    PINCP_NUM = as.numeric(PINCP), POVPIP = as.numeric(POVPIP), SEX = factor(SEX),
    POV = dplyr::if_else(POVPIP <= 100, 1, 0),
    BACH = factor(dplyr::if_else(SCHL >= 21, 1, 0))
  ) |>
  dplyr::filter(is.finite(PINCP_NUM), PINCP_NUM > 0) |>
  dplyr::arrange(PUMA, SEX, BACH)
stopifnot(all(is.finite(pums$PWGTP) & pums$PWGTP > 0))
population_poverty <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(poverty = mean(POV), .groups = "drop") |>
  dplyr::arrange(PUMA)
areas <- population_poverty$PUMA
# Population truth: Calculate each PUMA original-income median.
truth <- data.frame(PUMA = areas, income = as.numeric(vapply(areas, function(area)
  stats::median(pums$PINCP_NUM[as.character(pums$PUMA) == area]), numeric(1))),
  poverty = population_poverty$poverty, stringsAsFactors = FALSE)

## Spatial basis and population prediction cells
geometry <- readRDS("./main/data/il_puma_2023_sf.rds")
geometry$PUMA <- sprintf("%05d", as.integer(geometry$PUMACE20))
geometry <- geometry[match(areas, geometry$PUMA), , drop = FALSE]
stopifnot(identical(geometry$PUMA, areas), !anyNA(geometry$PUMA))
adjacency <- spdep::nb2mat(spdep::poly2nb(geometry, queen = TRUE),
  style = "B", zero.policy = TRUE)
adjacency <- (adjacency + t(adjacency)) / 2
eig <- eigen(adjacency, symmetric = TRUE)
basis <- eig$vectors[, eig$values > 0, drop = FALSE]
rownames(basis) <- areas
stopifnot(ncol(basis) > 0L)
pcells <- pums |>
  dplyr::count(PUMA, SEX, BACH, name = "popsize") |>
  dplyr::arrange(PUMA, SEX, BACH)
stopifnot(!anyNA(pcells), all(pcells$popsize > 0),
  identical(unique(pcells$PUMA), areas))
predX <- model.matrix(~ SEX + BACH - 1, pcells)
predS <- basis[match(pcells$PUMA, rownames(basis)), , drop = FALSE]
cell_area <- match(as.character(pcells$PUMA), areas)
totals <- as.numeric(rowsum(pcells$popsize, cell_area, reorder = FALSE))
area_weights <- Matrix::sparseMatrix(i = cell_area, j = seq_len(nrow(pcells)),
  x = pcells$popsize / totals[cell_area], dims = c(length(areas), nrow(pcells)))
population_counts <- rowsum(pcells$popsize, as.character(pcells$PUMA))
population_sizes <- setNames(as.numeric(population_counts[, 1L]), rownames(population_counts))
# PPS oversamples poor units; inverse inclusion weights are normalized below for fitting.
prob <- sampling::inclusionprobabilities(
  pums$PWGTP * (1 + 5 * (pums$POV == 1)), sample_size)
record_list <- vector("list", n_rep)

## Sample, fit all three models and predict each replicate
for (r in seq_len(n_rep)) {
  set.seed(r)
  sample_rows <- which(as.logical(sampling::UPsystematic(prob)))
  samp <- pums[sample_rows, , drop = FALSE]
  samp$W <- 1 / prob[sample_rows]
  # Estimate lambda from this sample's survey-weighted working-regression profile.
  lambda_info <- sample_boxcox_lambda(samp$PINCP_NUM,
    samp[, c("SEX", "BACH"), drop = FALSE], samp$W, grid = lambda_grid)
  lambda <- lambda_info$lambda
  check_prediction_lambda(lambda)
  samp$INCO <- boxcox_transform(samp$PINCP_NUM, lambda)
  X <- model.matrix(~ SEX + BACH - 1, samp)
  S <- basis[match(samp$PUMA, rownames(basis)), , drop = FALSE]
  wgt <- samp$W * nrow(samp) / sum(samp$W)
  direct_income <- md_direct(samp, areas)
  direct_poverty <- emp_poverty_direct_ht(samp, areas, population_sizes)
  message("Empirical replicate ", r, "/", n_rep, "; sample lambda = ", lambda)

  ## Gaussian model and all-unit income medians
  fit_ug <- uni_gaus(X = X, Y = samp$INCO, S = S,
    sig2b = pri$beta_variance, wgt = wgt, predX = predX, predS = predS,
    nburn = nburn, nsim = nsim, nthin = nthin,
    a = pri$gaussian_area_shape, b = pri$gaussian_area_rate,
    a_eps = pri$gaussian_error_shape, b_eps = pri$gaussian_error_rate,
    on_draw = NULL)
  income_ug <- pure_r_empirical_medians(fit_ug$Preds,
    as.numeric(fit_ug$sig2.chain), lambda, pcells$popsize, cell_area,
    chunk_units = chunk_units)
  rownames(income_ug$median) <- areas

  ## Binomial model and all-unit poverty proportions
  fit_ub <- uni_bio(X = X, Y = samp$POV, S = S,
    sig2b = pri$beta_variance, wgt = wgt, predX = predX, predS = predS,
    nburn = nburn, nsim = nsim, nthin = nthin,
    a = pri$binary_area_shape, b = pri$binary_area_rate, on_draw = NULL)
  poverty_cells_ub <- predict_cell_poverty(fit_ub$Preds, pcells$popsize)
  poverty_ub <- poststratify_cells(poverty_cells_ub, area_weights, areas)

  ## Joint model with paired income and poverty posterior columns
  fit_mt <- multi(X_1 = X, X_2 = X,
    Z_1 = samp$INCO, Z_2 = samp$POV, S = S,
    sig2b = pri$beta_variance, wgt = wgt, predX = predX, predS = predS,
    nburn = nburn, nsim = nsim, nthin = nthin,
    sig2t = pri$tau_variance, sig2e = 10, tau_2_init = 0,
    a_eps = pri$gaussian_error_shape, b_eps = pri$gaussian_error_rate,
    aeta = pri$shared_shape, beta = pri$shared_rate,
    alambda = pri$binary_specific_shape, blambda = pri$binary_specific_rate,
    on_draw = NULL)
  # Income and poverty predictions use the same joint posterior column.
  income_mt <- pure_r_empirical_medians(fit_mt$preds_gaus.chain,
    as.numeric(fit_mt$sig2.chain), lambda, pcells$popsize, cell_area,
    chunk_units = chunk_units, probability = fit_mt$preds_bios.chain)
  rownames(income_mt$median) <- areas
  poverty_mt <- poststratify_cells(income_mt$poverty, area_weights, areas)

  ## Original-scale posterior summaries and Direct comparisons
  income_uni <- summarize_draws(income_ug$median)
  income_uni$response <- "income"; income_uni$model <- "Univariate"
  poverty_uni <- summarize_draws(poverty_ub)
  poverty_uni$response <- "poverty"; poverty_uni$model <- "Univariate"
  income_multi <- summarize_draws(income_mt$median)
  income_multi$response <- "income"; income_multi$model <- "Multi-type"
  poverty_multi <- summarize_draws(poverty_mt)
  poverty_multi$response <- "poverty"; poverty_multi$model <- "Multi-type"
  columns <- c("PUMA", "response", "model", "estimate", "lower", "upper", "variance")
  records <- rbind(direct_income[columns], direct_poverty[columns], income_uni[columns],
    poverty_uni[columns], income_multi[columns], poverty_multi[columns])
  records$replicate <- r
  rownames(records) <- NULL
  record_list[[r]] <- records

  saveRDS(list(sample_rows = sample_rows, sample = samp, lambda = lambda_info,
    models = list(ug = fit_ug, ub = fit_ub, mt = fit_mt),
    income_draws = list(ug = income_ug$median, mt = income_mt$median),
    poverty_draws = list(ub = poverty_ub, mb = poverty_mt),
    audit = list(ug = income_ug$audit, mt = income_mt$audit), records = records),
    file.path(output_dir, sprintf("replicate_%03d.rds", r)))
}

## Evaluate over replicates and PUMAs
records <- do.call(rbind, record_list)
index <- match(as.character(records$PUMA), as.character(truth$PUMA))
stopifnot(!anyNA(index), !anyDuplicated(records[c("PUMA", "response", "model", "replicate")]),
  nrow(records) == length(areas) * n_rep * 6L)
target <- ifelse(records$response == "income", truth$income[index], truth$poverty[index])
scored <- cbind(records, score_intervals(records$estimate, records$lower, records$upper, target))
area_summary <- aggregate(scored[c("squared_error", "interval_score", "covered")],
  scored[c("PUMA", "response", "model")], mean)
names(area_summary)[match(c("squared_error", "interval_score", "covered"), names(area_summary))] <-
  c("MSE", "IS", "CR")
area_summary$Ratio_to_direct <- NA_real_
responses <- c("poverty", "income")
models <- c("Direct", "Univariate", "Multi-type")
for (response in responses) {
  direct <- area_summary[area_summary$response == response & area_summary$model == "Direct", ]
  stopifnot(all(direct$MSE > 0))
  take <- area_summary$response == response
  area_summary$Ratio_to_direct[take] <- area_summary$MSE[take] /
    direct$MSE[match(area_summary$PUMA[take], direct$PUMA)]
}
table_rows <- vector("list", 6L)
i <- 0L
for (response in responses) {
  for (model in models) {
    i <- i + 1L
    one <- area_summary[area_summary$response == response & area_summary$model == model, ]
    table_rows[[i]] <- data.frame(response = response, model = model,
      MSE = mean(one$MSE), Ratio_to_direct = mean(one$Ratio_to_direct),
      IS = mean(one$IS), CR = mean(one$CR), n_PUMA = nrow(one), stringsAsFactors = FALSE)
  }
}
summary_table <- do.call(rbind, table_rows)
rownames(summary_table) <- NULL
display_table <- summary_table
display_table$MSE <- ifelse(display_table$response == "income",
  display_table$MSE / 1e6, display_table$MSE * 1e3)
display_table$MSE_units <- ifelse(display_table$response == "income",
  "million dollars squared", "1e-3 proportion squared")
display_table$IS_units <- ifelse(display_table$response == "income", "dollars", "proportion")

## Save the paper Table 2
write.csv(summary_table, file.path(output_dir, "table2_original_scale.csv"), row.names = FALSE)
write.csv(display_table, file.path(output_dir, "table2_display_units.csv"), row.names = FALSE)
latex_rows <- character(nrow(display_table))
for (i in seq_len(nrow(display_table))) {
  one <- display_table[i, ]
  latex_rows[i] <- sprintf("%s & %s & %.3f & %.3f & %s & %.3f \\\\",
    if (one$response == "income") "Median income" else "Bernoulli", one$model,
    one$MSE, one$Ratio_to_direct,
    if (one$response == "income") formatC(one$IS, format = "f", digits = 1, big.mark = ",")
    else sprintf("%.3f", one$IS), one$CR)
}
writeLines(c("% Requires booktabs.", "\\begin{table}[htbp!]", "\\centering\\small",
  paste0("\\caption{Empirical simulation of all-unit finite-population median income and poverty proportions. ",
    "Results use ", n_rep, " samples of ", sample_size, " people and ", length(areas), " PUMAs. ",
    "Each income draw predicts all units, including sampled units, and takes the median after inverse Box--Cox transformation. ",
    "Point estimates are posterior means of median draws; intervals use their 2.5th and 97.5th percentiles. ",
    "Direct income uses the weighted median and Woodruff interval; Direct poverty is the HT total divided by known population size. ",
    "Metrics average equally over PUMAs. Ratio to Direct averages PUMA-specific MSE ratios. ",
    "Income MSE is in million dollars squared; poverty MSE is in $10^{-3}$.}"),
  "\\label{tab:empirical_all_unit_median}", "\\begin{tabular}{llrrrr}", "\\toprule",
  "Response & Model & MSE & Ratio to Direct & IS & CR \\\\", "\\midrule",
  latex_rows[1:3], "\\addlinespace", latex_rows[4:6],
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), file.path(output_dir, "table2.tex"))
saveRDS(list(summary = summary_table, display = display_table, truth = truth,
  records = records, area_summary = area_summary), file.path(output_dir, "empirical_results.rds"))
print(display_table, row.names = FALSE, digits = 12)
