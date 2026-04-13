source("basis version/packages.r")
suppressPackageStartupMessages(library(survey))

data_file <- "source data/IL23.csv"

dir.create("data", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

interval_score <- function(lower, upper, truth, alpha = 0.05) {
  width <- upper - lower
  penalty_low  <- (2 / alpha) * pmax(lower - truth, 0)
  penalty_high <- (2 / alpha) * pmax(truth - upper, 0)
  width + penalty_low + penalty_high
}

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

pums_raw <- readr::read_csv(data_file, show_col_types = FALSE)

pums <- pums_raw |>
  dplyr::select(PUMA, PWGTP, SEX, POVPIP, PINCP, SCHL) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    PUMA = sprintf("%05d", as.integer(PUMA)),
    PWGTP = as.numeric(PWGTP),
    POVPIP = as.numeric(POVPIP),
    SEX = factor(SEX),
    LOGINCO = log(as.numeric(PINCP)),
    POV = dplyr::if_else(POVPIP <= 100, 1, 0),
    BACH = factor(dplyr::if_else(SCHL >= 21, 1, 0))
  ) |>
  dplyr::filter(is.finite(LOGINCO)) |>
  dplyr::mutate(
    INCO = (LOGINCO - min(LOGINCO)) / (max(LOGINCO) - min(LOGINCO))
  ) |>
  arrange(PUMA, SEX, BACH)

truth <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(
    INCO = mean(INCO),
    POV = mean(POV),
    .groups = "drop"
  ) |>
  dplyr::arrange(PUMA)

area_levels <- truth$PUMA

n_sim <- 100
sample_size <- 1000
n_area <- nrow(truth)
zcrit <- qnorm(0.975)

dgaus_pre  <- matrix(NA_real_, nrow = n_area, ncol = n_sim)
dbio_pre   <- matrix(NA_real_, nrow = n_area, ncol = n_sim)
dgaus_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))
dbio_qual  <- array(NA_real_, dim = c(n_area, 2, n_sim))
cor_set    <- numeric(n_sim)

for (k in seq_len(n_sim)) {
  set.seed(k)

  prob <- sampling::inclusionprobabilities(
    pums$PWGTP * (1 + 5 * (pums$POV == 1)),
    sample_size
  )

  Ind <- sampling::UPsystematic(prob)
  samp <- pums[as.logical(Ind), , drop = FALSE]
  samp$P <- prob[as.logical(Ind)]
  samp$W <- 1 / samp$P

  cor_set[k] <- stats::cor(samp$POV, samp$INCO)

  des <- survey::svydesign(
    ids = ~1,
    weights = ~W,
    data = samp
  )

  direct_df <- survey::svyby(
    ~INCO + POV,
    ~PUMA,
    design = des,
    FUN = survey::svymean,
    vartype = "se",
    na.rm = TRUE,
    keep.names = TRUE,
    drop.empty.groups = FALSE
  ) |>
    as.data.frame() |>
    dplyr::mutate(
      PUMA = sprintf("%05d", as.integer(as.character(PUMA))),
      dgaus_est = INCO,
      dbio_est  = POV,
      dgaus_lb  = pmax(0, pmin(1, INCO - zcrit * se.INCO)),
      dgaus_ub  = pmax(0, pmin(1, INCO + zcrit * se.INCO)),
      dbio_lb   = pmax(0, pmin(1, POV  - zcrit * se.POV)),
      dbio_ub   = pmax(0, pmin(1, POV  + zcrit * se.POV))
    )

  idx <- match(area_levels, direct_df$PUMA)

  dgaus_pre[, k] <- direct_df$dgaus_est[idx]
  dbio_pre[, k]  <- direct_df$dbio_est[idx]

  dgaus_qual[, 1, k] <- direct_df$dgaus_lb[idx]
  dgaus_qual[, 2, k] <- direct_df$dgaus_ub[idx]
  dbio_qual[, 1, k]  <- direct_df$dbio_lb[idx]
  dbio_qual[, 2, k]  <- direct_df$dbio_ub[idx]

  cat("Finished", k, "simulation dataset\n")
}

mse_db_area <- numeric(n_area)
mse_dg_area <- numeric(n_area)
cr_db_area  <- numeric(n_area)
cr_dg_area  <- numeric(n_area)
is_db_area  <- numeric(n_area)
is_dg_area  <- numeric(n_area)

for (j in seq_len(n_area)) {
  mse_db_area[j] <- mean((dbio_pre[j, ]  - truth$POV[j])^2,  na.rm = TRUE)
  mse_dg_area[j] <- mean((dgaus_pre[j, ] - truth$INCO[j])^2, na.rm = TRUE)

  cr_db_area[j] <- mean(
    dbio_qual[j, 1, ] <= truth$POV[j] & truth$POV[j] <= dbio_qual[j, 2, ],
    na.rm = TRUE
  )
  cr_dg_area[j] <- mean(
    dgaus_qual[j, 1, ] <= truth$INCO[j] & truth$INCO[j] <= dgaus_qual[j, 2, ],
    na.rm = TRUE
  )

  is_db_area[j] <- mean(
    interval_score(dbio_qual[j, 1, ], dbio_qual[j, 2, ], truth$POV[j]),
    na.rm = TRUE
  )
  is_dg_area[j] <- mean(
    interval_score(dgaus_qual[j, 1, ], dgaus_qual[j, 2, ], truth$INCO[j]),
    na.rm = TRUE
  )
}

mse_db_rep <- colMeans((dbio_pre  - truth$POV)^2,  na.rm = TRUE)
mse_dg_rep <- colMeans((dgaus_pre - truth$INCO)^2, na.rm = TRUE)

cr_db_rep <- sapply(seq_len(n_sim), function(k) {
  mean(dbio_qual[, 1, k] <= truth$POV & truth$POV <= dbio_qual[, 2, k], na.rm = TRUE)
})

cr_dg_rep <- sapply(seq_len(n_sim), function(k) {
  mean(dgaus_qual[, 1, k] <= truth$INCO & truth$INCO <= dgaus_qual[, 2, k], na.rm = TRUE)
})

is_db_rep <- sapply(seq_len(n_sim), function(k) {
  mean(interval_score(dbio_qual[, 1, k], dbio_qual[, 2, k], truth$POV), na.rm = TRUE)
})

is_dg_rep <- sapply(seq_len(n_sim), function(k) {
  mean(interval_score(dgaus_qual[, 1, k], dgaus_qual[, 2, k], truth$INCO), na.rm = TRUE)
})

summary_table_area <- tibble::tibble(
  Response = c("Bernoulli", "Gaussian"),
  Model = c("Direct", "Direct"),
  MSE = c(
    safe_mean(mse_db_area),
    safe_mean(mse_dg_area)
  ),
  Mean_IS = c(
    safe_mean(is_db_area),
    safe_mean(is_dg_area)
  ),
  Mean_CR = c(
    safe_mean(cr_db_area),
    safe_mean(cr_dg_area)
  ),
  Coverage_Gap = c(
    safe_mean(abs(cr_db_area - 0.95)),
    safe_mean(abs(cr_dg_area - 0.95))
  )
)

summary_table_rep <- tibble::tibble(
  Response = c("Bernoulli", "Gaussian"),
  Model = c("Direct", "Direct"),
  MSE = c(
    safe_mean(mse_db_rep),
    safe_mean(mse_dg_rep)
  ),
  Mean_IS = c(
    safe_mean(is_db_rep),
    safe_mean(is_dg_rep)
  ),
  Mean_CR = c(
    safe_mean(cr_db_rep),
    safe_mean(cr_dg_rep)
  ),
  Coverage_Gap = c(
    safe_mean(abs(cr_db_rep - 0.95)),
    safe_mean(abs(cr_dg_rep - 0.95))
  )
)

direct_area_metrics <- tibble::tibble(
  PUMA = area_levels,
  Bernoulli_MSE = mse_db_area,
  Bernoulli_IS = is_db_area,
  Bernoulli_CR = cr_db_area,
  Bernoulli_Coverage_Gap = abs(cr_db_area - 0.95),
  Gaussian_MSE = mse_dg_area,
  Gaussian_IS = is_dg_area,
  Gaussian_CR = cr_dg_area,
  Gaussian_Coverage_Gap = abs(cr_dg_area - 0.95)
)

direct_rep_metrics <- tibble::tibble(
  rep = seq_len(n_sim),
  Bernoulli_MSE = mse_db_rep,
  Bernoulli_IS = is_db_rep,
  Bernoulli_CR = cr_db_rep,
  Bernoulli_Coverage_Gap = abs(cr_db_rep - 0.95),
  Gaussian_MSE = mse_dg_rep,
  Gaussian_IS = is_dg_rep,
  Gaussian_CR = cr_dg_rep,
  Gaussian_Coverage_Gap = abs(cr_dg_rep - 0.95)
)

print(summary_table_area)
print(summary_table_rep)
print(direct_area_metrics)
print(direct_rep_metrics)

# save(
#   truth,
#   area_levels,
#   dgaus_pre,
#   dbio_pre,
#   dgaus_qual,
#   dbio_qual,
#   mse_db_area,
#   mse_dg_area,
#   cr_db_area,
#   cr_dg_area,
#   is_db_area,
#   is_dg_area,
#   mse_db_rep,
#   mse_dg_rep,
#   cr_db_rep,
#   cr_dg_rep,
#   is_db_rep,
#   is_dg_rep,
#   cor_set,
#   summary_table_area,
#   summary_table_rep,
#   direct_area_metrics,
#   direct_rep_metrics,
#   file = file.path("data", "empirical_basis_direct_only_results.RData")
# )

# readr::write_csv(summary_table_area, "basis_direct_summary_table_area.csv")
# readr::write_csv(summary_table_rep, "basis_direct_summary_table_rep.csv")
# readr::write_csv(direct_area_metrics, "basis_direct_area_metrics.csv")
# readr::write_csv(direct_rep_metrics, "basis_direct_rep_metrics.csv")