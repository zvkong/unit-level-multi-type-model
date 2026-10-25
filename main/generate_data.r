library(dplyr)
library(tidyr)
library(sf)
library(tigris)
library(jsonlite)

## Settings
data_dir <- "./main/data"
census_key <- Sys.getenv("CENSUS_API_KEY")
if (!nzchar(census_key)) stop("Set CENSUS_API_KEY before generating the ACS population data.")
options(timeout = 600, tigris_use_cache = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
download_dir <- tempfile("illinois_acs_2023_")
dir.create(download_dir)

## Original 2023 one-year Illinois person records
pums_url <- "https://www2.census.gov/programs-surveys/acs/data/pums/2023/1-Year/csv_pil.zip"
pums_zip <- file.path(download_dir, "csv_pil.zip")
download.file(pums_url, pums_zip, mode = "wb")
unzip(pums_zip, files = "psam_p17.csv", exdir = download_dir)
# Copy the original CSV without changing its columns, values or row order.
stopifnot(file.copy(file.path(download_dir, "psam_p17.csv"),
  file.path(data_dir, "IL23.csv"), overwrite = TRUE))

## Empirical study: 2023 TIGER/Line PUMA boundaries
puma_2023 <- tigris::pumas(state = "IL", cb = FALSE, year = 2023, class = "sf") |>
  sf::st_transform(4326)
stopifnot(nrow(puma_2023) == 88L, "PUMACE20" %in% names(puma_2023))
saveRDS(puma_2023, file.path(data_dir, "il_puma_2023_sf.rds"))

## Region study: 2020 cartographic PUMA boundaries
puma_2020 <- tigris::pumas(state = "IL", cb = TRUE, year = 2020, class = "sf") |>
  sf::st_transform(4326)
stopifnot(nrow(puma_2020) == 88L, "PUMACE20" %in% names(puma_2020))
saveRDS(puma_2020, file.path(data_dir, "il_puma_2020_cb_sf.rds"))

## B15001 mappings: sex totals and bachelor's/graduate degrees across adult age groups
maps <- list(
  sex_total_vars = tibble::tibble(
    var = c("B15001_002", "B15001_043"),
    SEX = c("Male", "Female")
  ),
  bach_vars = tibble::tibble(
    var = paste0("B15001_", sprintf("%03d", c(
      9, 10, 17, 18, 25, 26, 33, 34, 41, 42,
      50, 51, 58, 59, 66, 67, 74, 75, 82, 83
    ))),
    SEX = rep(c("Male", "Female"), each = 10)
  )
)
saveRDS(maps, file.path(data_dir, "acs2023_b15001_var_maps.rds"))

## Original 2023 one-year ACS population estimates for Illinois PUMAs
vars_needed <- c(maps$sex_total_vars$var, maps$bach_vars$var)
var_api <- paste0(vars_needed, "E")
api_url <- paste0(
  "https://api.census.gov/data/2023/acs/acs1?get=NAME,",
  paste(var_api, collapse = ","),
  "&for=public%20use%20microdata%20area:*&in=state:17&key=",
  utils::URLencode(census_key, reserved = TRUE)
)
response <- tryCatch(
  suppressWarnings(jsonlite::fromJSON(api_url)),
  # API error: report the failure without displaying the API key.
  error = function(e) stop("Census request failed. Check CENSUS_API_KEY and the connection.",
    call. = FALSE)
)
body <- as.data.frame(response[-1L, , drop = FALSE], stringsAsFactors = FALSE)
names(body) <- as.character(response[1L, ])
raw_cells <- body |>
  dplyr::mutate(GEOID = paste0("17", .data[["public use microdata area"]])) |>
  dplyr::select(GEOID, NAME, dplyr::all_of(var_api)) |>
  tidyr::pivot_longer(dplyr::all_of(var_api), names_to = "variable", values_to = "estimate") |>
  dplyr::mutate(variable = sub("E$", "", variable), estimate = as.numeric(estimate)) |>
  dplyr::select(GEOID, NAME, variable, estimate) |>
  dplyr::arrange(GEOID, factor(variable, levels = vars_needed))
stopifnot(nrow(raw_cells) == 1936L, length(unique(raw_cells$GEOID)) == 88L,
  all(is.finite(raw_cells$estimate) & raw_cells$estimate >= 0),
  sum(raw_cells$estimate[raw_cells$variable %in% maps$sex_total_vars$var]) == 9848062)
saveRDS(raw_cells, file.path(data_dir, "acs2023_il_puma_b15001_raw_cells.rds"))

message("Generated all five input files in ", data_dir)
