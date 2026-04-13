source("basis version/packages.r")
source("basis version/functions.R")

data_file <- "source data/IL23.csv"

dir.create("data", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

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

pcells <- pums |>
  dplyr::group_by(PUMA, SEX, BACH) |>
  dplyr::summarise(popsize = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(factor(PUMA, levels = area_levels), SEX, BACH)

predX <- model.matrix(~ SEX + BACH - 1, data = pcells)

puma_sf <- tigris::pumas(state = "IL", year = 2023, class = "sf") |>
  sf::st_transform(4326) |>
  dplyr::mutate(PUMA = sprintf("%05d", as.integer(PUMACE20))) |>
  dplyr::filter(PUMA %in% area_levels) |>
  dplyr::arrange(factor(PUMA, levels = area_levels))

if (!all(area_levels %in% puma_sf$PUMA)) {
  stop("Some PUMAs in truth are missing from the shapefile.")
}

nb <- spdep::poly2nb(puma_sf, queen = TRUE)
W_mat <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)

eig <- eigen(W_mat, symmetric = TRUE)
pos_id <- which(eig$values > 0)

if (length(pos_id) == 0) {
  stop("No positive eigenvalues found for basis construction.")
}

basis_mat <- eig$vectors[, pos_id, drop = FALSE]
rownames(basis_mat) <- puma_sf$PUMA

predPsi <- basis_mat[match(pcells$PUMA, rownames(basis_mat)), , drop = FALSE]

n_sim <- 100
nsim <- 1000
nburn <- 1000
nthin <- 1
sample_size <- 1000

n_area <- nrow(truth)

ubios_pre <- array(NA_real_, dim = c(n_area, n_sim))
mbios_pre <- array(NA_real_, dim = c(n_area, n_sim))
ugaus_pre <- array(NA_real_, dim = c(n_area, n_sim))
mgaus_pre <- array(NA_real_, dim = c(n_area, n_sim))
dgaus_pre <- array(NA_real_, dim = c(n_area, n_sim))
dbio_pre <- array(NA_real_, dim = c(n_area, n_sim))

ubios_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))
mbios_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))
ugaus_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))
mgaus_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))

cor_set <- numeric(n_sim)

for (k in seq_len(n_sim)) {
  set.seed(k)

  prob <- sampling::inclusionprobabilities(
    pums$PWGTP * (1 + 5*(pums$POV == 1)),
    sample_size
  )

  Ind <- sampling::UPsystematic(prob)
  samp <- pums[as.logical(Ind), , drop = FALSE]
  samp$P <- prob[as.logical(Ind)]
  samp$W <- 1 / samp$P
  samp$scaledWGT <- samp$W * nrow(samp) / sum(samp$W)

  compare_df <- samp |>
    dplyr::group_by(PUMA) |>
    dplyr::summarise(
      weighted_means = stats::weighted.mean(INCO, W),
      weighted_means_POV = stats::weighted.mean(POV, W),
      .groups = "drop"
    )

  dgaus_pre[, k] <- compare_df$weighted_means[match(area_levels, compare_df$PUMA)]
  dbio_pre[, k] <- compare_df$weighted_means_POV[match(area_levels, compare_df$PUMA)]

  cor_set[k] <- stats::cor(samp$POV, samp$INCO)

  modwgt <- samp$scaledWGT
  modX <- model.matrix(~ SEX + BACH - 1, data = samp)
  modPsi <- basis_mat[match(samp$PUMA, rownames(basis_mat)), , drop = FALSE]
  modY <- samp$INCO
  modZ <- samp$POV

  fit_ug <- unis_gaus_fast(
    X = modX,
    Y = modY,
    S = modPsi,
    sig2b = 1000,
    wgt = modwgt,
    predX = predX,
    predS = predPsi,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    a = 0.5,
    b = 0.5,
    a_eps = 0.1,
    b_eps = 0.1
  )

  fit_ub <- unis_bios_fast(
    X = modX,
    Y = modZ,
    S = modPsi,
    sig2b = 1000,
    wgt = modwgt,
    n = NULL,
    predX = predX,
    predS = predPsi,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    a = 0.1,
    b = 0.1
  )

  fit_mt <- MTSM_br_tau2_gaus_fast(
    X_1 = modX,
    X_2 = modX,
    Z_1 = modY,
    Z_2 = modZ,
    S = modPsi,
    sig2b = 1000,
    wgt = modwgt,
    n = NULL,
    predX = predX,
    predS = predPsi,
    n_preds = NULL,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    sig2t = 10,
    sig2e = 10,
    tau_2_init = 0,
    a_eps = 0.1,
    b_eps = 0.1,
    aeta = 0.1,
    beta = 0.1,
    alambda = 2,
    blambda = 1
  )

  res_ug <- gaus_post(
    fit_ug$Preds,
    fit_ug$sig2.chain,
    truth$INCO,
    pcells$PUMA,
    pcells$popsize
  )

  res_mg <- gaus_post(
    fit_mt$preds_gaus.chain,
    fit_mt$sig2.chain,
    truth$INCO,
    pcells$PUMA,
    pcells$popsize
  )

  res_ub <- bios_post(
    fit_ub$Preds,
    truth$POV,
    pcells$PUMA,
    pcells$popsize
  )

  res_mb <- bios_post(
    fit_mt$preds_bios.chain,
    truth$POV,
    pcells$PUMA,
    pcells$popsize
  )

  ugaus_pre[, k] <- res_ug$est
  mgaus_pre[, k] <- res_mg$est
  ubios_pre[, k] <- res_ub$est
  mbios_pre[, k] <- res_mb$est

  ugaus_qual[, 1, k] <- res_ug$lb
  ugaus_qual[, 2, k] <- res_ug$ub
  mgaus_qual[, 1, k] <- res_mg$lb
  mgaus_qual[, 2, k] <- res_mg$ub
  ubios_qual[, 1, k] <- res_ub$lb
  ubios_qual[, 2, k] <- res_ub$ub
  mbios_qual[, 1, k] <- res_mb$lb
  mbios_qual[, 2, k] <- res_mb$ub

  cat("Finished", k, "simulation dataset\n")
}

save(
  truth,
  area_levels,
  basis_mat,
  puma_sf,
  pcells,
  ubios_pre, mbios_pre, ugaus_pre, mgaus_pre, dbio_pre, dgaus_pre,
  ubios_qual, mbios_qual, ugaus_qual, mgaus_qual,
  cor_set,
  file = file.path("data", "empirical_basis_tau_in_gauss_results 2.RData")
)