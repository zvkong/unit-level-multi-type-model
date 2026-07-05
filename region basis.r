## Inferred generator for data/region_basis_tau_in_gauss_results.RData
##
## This script is reconstructed from the fingerprints observed in the saved
## RData workspace.  It is intentionally not identical to the current
## older region basis.r file.
##
## Main inferred differences:
## - pums uses the older compact column layout:
##   PUMA, PWGTP, SEX, BACH, INCO, POV, scaledWGT
## - INCO uses min-max scaling, not z-score scaling.
## - POV uses POVPIP <= 100, matching the saved POV count.
## - pcells are built from Census ACS B15001 cells, not by grouping pums.
## - the output uses a workspace-style save of analysis objects, matching the
##   many extra objects found in the stored RData.
## - basis_mat is unique only up to eigenvector column signs.  This script does
##   not use the existing RData to force signs; exact signs may differ across
##   R/LAPACK builds while representing the same basis space.
##
## run_models controls whether the script runs MCMC and saves the region-basis
## RData used by region basis results.r.

required_pkgs <- c(
  "readr", "dplyr", "tidyr", "tibble", "stringr",
  "tidycensus", "tigris", "sf", "spdep", "jsonlite"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(tidycensus)
  library(tigris)
  library(sf)
  library(spdep)
  library(jsonlite)
})

source("functions.r")

data_file <- "source data/IL23.csv"
out_dir <- "data"
out_file <- file.path(out_dir, "region_basis_tau_in_gauss_results.RData")
fig_dir <- "figs"
run_models <- TRUE

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

fetch_acs_raw_cells <- function(vars_needed, year = 2023, state = "17") {
  var_api <- paste0(vars_needed, "E")
  key <- Sys.getenv("CENSUS_API_KEY")
  key_query <- if (nzchar(key)) paste0("&key=", key) else ""
  api_url <- paste0(
    "https://api.census.gov/data/", year, "/acs/acs1",
    "?get=NAME,", paste(var_api, collapse = ","),
    "&for=public%20use%20microdata%20area:*",
    "&in=state:", state,
    key_query
  )

  acs_json <- jsonlite::fromJSON(api_url)
  header <- as.character(acs_json[1, ])
  body <- as.data.frame(acs_json[-1, , drop = FALSE], stringsAsFactors = FALSE)
  names(body) <- header

  body |>
    dplyr::mutate(GEOID = paste0(state, .data[["public use microdata area"]])) |>
    dplyr::select(GEOID, NAME, dplyr::all_of(var_api)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(var_api),
      names_to = "variable",
      values_to = "estimate"
    ) |>
    dplyr::mutate(
      variable = sub("E$", "", variable),
      estimate = as.numeric(estimate)
    ) |>
    dplyr::select(GEOID, NAME, variable, estimate) |>
    dplyr::arrange(GEOID, factor(variable, levels = vars_needed))
}

pums_raw <- readr::read_csv(data_file, show_col_types = FALSE)

pums <- pums_raw |>
  dplyr::select(PUMA, PWGTP, AGEP, SEX, POVPIP, PINCP, SCHL) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    PUMA = sprintf("%05d", as.integer(PUMA)),
    PWGTP = as.numeric(PWGTP),
    AGEP = as.numeric(AGEP),
    POVPIP = as.numeric(POVPIP),
    SEX = factor(SEX),
    LOGINCO = log(as.numeric(PINCP)),
    POV = dplyr::if_else(POVPIP <= 100, 1, 0),
    BACH = factor(dplyr::if_else(SCHL >= 21, 1, 0))
  ) |>
  dplyr::filter(AGEP >= 18, is.finite(LOGINCO)) |>
  dplyr::mutate(
    INCO = (LOGINCO - min(LOGINCO)) / (max(LOGINCO) - min(LOGINCO)),
    scaledWGT = PWGTP * dplyr::n() / sum(PWGTP)
  ) |>
  dplyr::select(PUMA, PWGTP, SEX, BACH, INCO, POV, scaledWGT) |>
  dplyr::arrange(PUMA, SEX, BACH)

puma_sf <- tigris::pumas(state = "IL", cb = TRUE, year = 2020, class = "sf") |>
  sf::st_transform(4326)

if ("PUMACE20" %in% names(puma_sf)) {
  puma_sf <- puma_sf |>
    dplyr::mutate(PUMA = sprintf("%05d", as.integer(PUMACE20)))
} else if ("PUMACE10" %in% names(puma_sf)) {
  puma_sf <- puma_sf |>
    dplyr::mutate(PUMA = sprintf("%05d", as.integer(PUMACE10)))
} else {
  stop("PUMA code column not found in shapefile.")
}

puma_sf <- puma_sf |>
  dplyr::filter(PUMA %in% unique(pums$PUMA)) |>
  dplyr::arrange(PUMA)

area_levels <- puma_sf$PUMA

nb <- spdep::poly2nb(puma_sf, queen = TRUE)
A <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
A <- (A + t(A)) / 2

eig <- eigen(A, symmetric = TRUE)
pos_idx <- which(eig$values > 1e-8)

basis_mat <- eig$vectors[, pos_idx, drop = FALSE]
rownames(basis_mat) <- area_levels
colnames(basis_mat) <- paste0("bf", seq_len(ncol(basis_mat)))

truth <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(
    INCO = mean(INCO),
    POV = mean(POV),
    .groups = "drop"
  ) |>
  dplyr::arrange(factor(PUMA, levels = area_levels))

area_levels <- truth$PUMA
basis_mat <- basis_mat[area_levels, , drop = FALSE]

v <- tidycensus::load_variables(2023, "acs1", cache = TRUE) |>
  dplyr::filter(stringr::str_starts(name, "B15001_"))

sex_total_vars <- v |>
  dplyr::filter(label %in% c("Estimate!!Total:!!Male:", "Estimate!!Total:!!Female:")) |>
  dplyr::transmute(
    var = name,
    SEX = dplyr::if_else(stringr::str_detect(label, "Male"), "Male", "Female")
  )

bach_vars <- v |>
  dplyr::filter(
    stringr::str_detect(label, "!!Male:|!!Female:"),
    stringr::str_detect(label, "Bachelor's degree|Graduate or professional degree")
  ) |>
  dplyr::transmute(
    var = name,
    SEX = dplyr::case_when(
      stringr::str_detect(label, "!!Male") ~ "Male",
      stringr::str_detect(label, "!!Female") ~ "Female",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(SEX))

vars_needed <- c(sex_total_vars$var, bach_vars$var)

raw_cells <- fetch_acs_raw_cells(vars_needed, year = 2023, state = "17")

sex_totals <- raw_cells |>
  dplyr::inner_join(sex_total_vars, by = c("variable" = "var")) |>
  dplyr::group_by(GEOID, NAME, SEX) |>
  dplyr::summarise(pop_total = sum(estimate), .groups = "drop")

bach <- raw_cells |>
  dplyr::inner_join(bach_vars, by = c("variable" = "var")) |>
  dplyr::group_by(GEOID, NAME, SEX) |>
  dplyr::summarise(pop_bach = sum(estimate), .groups = "drop")

cells <- sex_totals |>
  dplyr::left_join(bach, by = c("GEOID", "NAME", "SEX")) |>
  dplyr::mutate(pop_bach = dplyr::coalesce(pop_bach, 0)) |>
  dplyr::transmute(
    PUMA = sprintf("%05d", as.integer(substr(GEOID, nchar(GEOID) - 4, nchar(GEOID)))),
    puma_name = NAME,
    SEX = factor(dplyr::if_else(SEX == "Male", "1", "2"), levels = levels(pums$SEX)),
    pop_bach = pop_bach,
    pop_bach0 = pmax(pop_total - pop_bach, 0)
  ) |>
  tidyr::pivot_longer(cols = c(pop_bach, pop_bach0), names_to = "BACH", values_to = "pop") |>
  dplyr::mutate(
    BACH = factor(dplyr::if_else(BACH == "pop_bach", "1", "0"), levels = levels(pums$BACH))
  ) |>
  dplyr::filter(PUMA %in% area_levels) |>
  dplyr::select(PUMA, SEX, BACH, pop, puma_name) |>
  dplyr::arrange(factor(PUMA, levels = area_levels), SEX, BACH)

pcells <- pums |>
  dplyr::distinct(PUMA, SEX, BACH) |>
  dplyr::left_join(cells, by = c("PUMA", "SEX", "BACH")) |>
  dplyr::mutate(popsize = pop) |>
  dplyr::select(PUMA, SEX, BACH, popsize) |>
  dplyr::arrange(factor(PUMA, levels = area_levels), SEX, BACH)

stopifnot(!anyNA(pcells$popsize))

predX <- model.matrix(~ SEX + BACH - 1, data = pcells)

modX <- model.matrix(~ SEX + BACH - 1, data = pums)
modY <- pums$INCO
modZ <- pums$POV

nburn <- 1000
nsim <- 1000
nthin <- 1

if (!run_models) {
  message("Pre-model objects have been generated.")
  message("Set run_models <- TRUE to run MCMC and save ", out_file)
} else {
  unis_wage <- unis_gaus_fast(
    X = modX,
    Y = modY,
    S = basis_mat[match(pums$PUMA, rownames(basis_mat)), , drop = FALSE],
    sig2b = 1000,
    wgt = pums$scaledWGT,
    predX = predX,
    predS = basis_mat[match(pcells$PUMA, rownames(basis_mat)), , drop = FALSE],
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    a = 0.5,
    b = 0.5,
    a_eps = 0.1,
    b_eps = 0.1
  )

  unis_pov <- unis_bios_fast(
    X = modX,
    Y = modZ,
    S = basis_mat[match(pums$PUMA, rownames(basis_mat)), , drop = FALSE],
    sig2b = 1000,
    wgt = pums$scaledWGT,
    n = NULL,
    predX = predX,
    predS = basis_mat[match(pcells$PUMA, rownames(basis_mat)), , drop = FALSE],
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    a = 0.1,
    b = 0.1
  )

  mult_br <- MTSM_br_tau2_gaus_fast(
    X_1 = modX,
    X_2 = modX,
    Z_1 = modY,
    Z_2 = modZ,
    S = basis_mat[match(pums$PUMA, rownames(basis_mat)), , drop = FALSE],
    sig2b = 1000,
    wgt = pums$scaledWGT,
    n = NULL,
    predX = predX,
    predS = basis_mat[match(pcells$PUMA, rownames(basis_mat)), , drop = FALSE],
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

  true_mean_dummy <- stats::setNames(rep(NA_real_, length(area_levels)), area_levels)

  res_ug <- gaus_post(unis_wage$Preds, unis_wage$sig2.chain, true_mean_dummy, pcells$PUMA, pcells$popsize)
  res_mg <- gaus_post(mult_br$preds_gaus.chain, mult_br$sig2.chain, true_mean_dummy, pcells$PUMA, pcells$popsize)
  res_ub <- bios_post(unis_pov$Preds, true_mean_dummy, pcells$PUMA, pcells$popsize)
  res_mb <- bios_post(mult_br$preds_bios.chain, true_mean_dummy, pcells$PUMA, pcells$popsize)

  ht_by_puma <- pums |>
    dplyr::group_by(PUMA) |>
    dplyr::summarise(
      ht_inco = stats::weighted.mean(INCO, PWGTP),
      ht_pov = stats::weighted.mean(POV, PWGTP),
      .groups = "drop"
    )

  gaus_mapdat <- puma_sf |>
    dplyr::left_join(
      tibble::tibble(PUMA = area_levels, univariate = res_ug$est, multitype = res_mg$est),
      by = "PUMA"
    ) |>
    dplyr::left_join(ht_by_puma |> dplyr::select(PUMA, ht_inco), by = "PUMA")

  bios_mapdat <- puma_sf |>
    dplyr::left_join(
      tibble::tibble(PUMA = area_levels, univariate = res_ub$est, multitype = res_mb$est),
      by = "PUMA"
    ) |>
    dplyr::left_join(ht_by_puma |> dplyr::select(PUMA, ht_pov), by = "PUMA")

  sigma2_gaus_mapdat <- puma_sf |>
    dplyr::left_join(
      tibble::tibble(PUMA = area_levels, univariate = as.numeric(res_ug$sigma2), multitype = as.numeric(res_mg$sigma2)),
      by = "PUMA"
    )

  sigma2_bios_mapdat <- puma_sf |>
    dplyr::left_join(
      tibble::tibble(PUMA = area_levels, univariate = as.numeric(res_ub$sigma2), multitype = as.numeric(res_mb$sigma2)),
      by = "PUMA"
    )

  fit_ug <- unis_wage
  fit_ub <- unis_pov
  fit_mt <- mult_br

  helper_objects <- c(
    "required_pkgs", "missing_pkgs", "out_dir", "out_file",
    "fetch_acs_raw_cells", "run_models", "helper_objects", "objects_to_save"
  )
  objects_to_save <- setdiff(ls(all.names = TRUE), helper_objects)
  save(list = objects_to_save, file = out_file)
}
