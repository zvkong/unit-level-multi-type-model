source("basis version/packages.r")
source("basis version/functions.r")
data_file <- "source data/IL23.csv"

dir.create("data", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

pums_raw <- readr::read_csv(data_file, show_col_types = FALSE)

keep_vars <- intersect(
  c("PUMA", "PWGTP", "SEX", "POVPIP", "PINCP", "SCHL", "AGEP"),
  names(pums_raw)
)

pums <- pums_raw |>
  dplyr::select(all_of(keep_vars)) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    PUMA = sprintf("%05d", as.integer(PUMA)),
    PWGTP = as.numeric(PWGTP),
    POVPIP = as.numeric(POVPIP),
    PINCP = as.numeric(PINCP),
    SCHL = as.numeric(SCHL),
    SEX = factor(SEX)
  )

if ("AGEP" %in% names(pums)) {
  pums <- pums |>
    dplyr::mutate(AGEP = as.numeric(AGEP)) |>
    dplyr::filter(AGEP >= 18)
}

pums <- pums |>
  dplyr::mutate(
    LOGINCO = log(PINCP),
    POV = dplyr::if_else(POVPIP < 100, 1, 0),
    BACH = factor(dplyr::if_else(SCHL >= 21, 1, 0))
  ) |>
  dplyr::filter(is.finite(LOGINCO)) |>
  dplyr::mutate(
    INCO = (LOGINCO - min(LOGINCO)) / (max(LOGINCO) - min(LOGINCO))
  ) |>
  dplyr::arrange(PUMA, SEX, BACH)

area_levels <- sort(unique(pums$PUMA))

if (file.exists("IL_poststrat_cells.csv")) {
  pcells <- readr::read_csv("IL_poststrat_cells.csv", show_col_types = FALSE) |>
    dplyr::mutate(
      PUMA = sprintf("%05d", as.integer(PUMA)),
      SEX = factor(SEX, levels = levels(pums$SEX)),
      BACH = factor(BACH, levels = levels(pums$BACH)),
      popsize = as.numeric(popsize)
    ) |>
    dplyr::filter(PUMA %in% area_levels) |>
    dplyr::arrange(factor(PUMA, levels = area_levels), SEX, BACH)
} else {
  pcells <- pums |>
    dplyr::group_by(PUMA, SEX, BACH) |>
    dplyr::summarise(popsize = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(factor(PUMA, levels = area_levels), SEX, BACH)
}

direct_area <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(
    ht_inco = stats::weighted.mean(INCO, PWGTP),
    ht_pov = stats::weighted.mean(POV, PWGTP),
    .groups = "drop"
  ) |>
  dplyr::arrange(PUMA)

predX <- model.matrix(~ SEX + BACH - 1, data = pcells)

puma_sf <- tigris::pumas(state = "IL", year = 2021, class = "sf") |>
  sf::st_transform(4326) |>
  dplyr::mutate(PUMA = sprintf("%05d", as.integer(PUMACE10))) |>
  dplyr::filter(PUMA %in% area_levels) |>
  dplyr::arrange(factor(PUMA, levels = area_levels))

if (!all(area_levels %in% puma_sf$PUMA)) {
  stop("Some PUMAs in analysis data are missing from the shapefile.")
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

modX <- model.matrix(~ SEX + BACH - 1, data = pums)
modPsi <- basis_mat[match(pums$PUMA, rownames(basis_mat)), , drop = FALSE]
modY <- pums$INCO
modZ <- pums$POV
modwgt <- pums$PWGTP * nrow(pums) / sum(pums$PWGTP)

nburn <- 1000
nsim <- 1000
nthin <- 1

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

save(
  pums,
  pcells,
  puma_sf,
  area_levels,
  basis_mat,
  direct_area,
  fit_ug,
  fit_ub,
  fit_mt,
  file = file.path("data", "region_basis_tau_in_gauss_results.RData")
)