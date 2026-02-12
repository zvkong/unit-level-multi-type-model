# Setup
source("packages.r")    # libraries
source("functions.r")   # model + poststrat functions

# Choose the variables you need (examples below). Replace/extend as needed.
pums_vars <- c(
  "PINCP",
  "PWGTP",     # person weight
  "PUMA",      # PUMA code (2010 PUMA for 2021 ACS)
  "AGEP",
  "SEX",
  "RAC1P",
  "SCHL",
  "POVPIP"     # poverty ratio
)

pums21 <- tidycensus::get_pums(
  variables = pums_vars,
  state = "IL",
  year = 2021,
  survey = "acs1",
  recode = TRUE
)
# Missouri PUMA shapefile (sf)
# For ACS 2021, geography corresponds to 2010 PUMAs (PUMACE10).
puma_sf <- tigris::pumas(
  state = "IL",
  year = 2021,
  class = "sf"
) %>%
  st_transform(4326) %>%
  mutate(
    PUMA = as.integer(PUMACE10)
  )

puma_sf <- dplyr::mutate(puma_sf, PUMA = as.character(.data$PUMACE10))
pums_samp <- pums21 %>% 
  select(PUMA, PWGTP, SEX, RAC1P, POVPIP, PINCP, SCHL, AGEP) %>%
  filter(AGEP >= 18) %>%
  tidyr::drop_na() |>
  dplyr::mutate(
  PWGTP   = as.numeric(PWGTP),
  SEX     = factor(SEX),
  LOGINCO = log(as.numeric(PINCP)),
  BACH    = dplyr::if_else(SCHL >= 21, 1L, 0L),
  POV     = dplyr::if_else(POVPIP <= 100, 1L, 0L)
  ) |>
  dplyr::filter(is.finite(LOGINCO)) |>
  dplyr::mutate(
  INCO = (LOGINCO - min(LOGINCO)) / (max(LOGINCO) - min(LOGINCO)),
  scaledWGT = PWGTP * dplyr::n() / sum(PWGTP)
  ) %>%
  select(PUMA, SEX, BACH, INCO, POV, PWGTP, scaledWGT) %>%
  arrange(PUMA, SEX, BACH)

nrow(pums_samp)
puma_sf <- tigris::pumas(
  state = "IL",
  year = 2021,
  class = "sf"
) %>%
  st_transform(4326) %>%
  mutate(
    PUMA = as.integer(PUMACE10)  # align key with PUMS "PUMA"
  )

puma_sf <- dplyr::mutate(puma_sf, PUMA = as.character(.data$PUMACE10))

## load the variables
v <- load_variables(2021, "acs1", cache = TRUE) %>%
  filter(str_starts(name, "B15001_"))
## total male and female
SEX_total_vars <- v %>%
  filter(label %in% 
           c("Estimate!!Total:!!Male:", "Estimate!!Total:!!Female:")) %>%
  transmute(var = name,
            SEX = if_else(str_detect(label, "Male"), "Male", "Female"))

bach_vars <- v %>%
  filter(str_detect(label, "!!Male:|!!Female:"),
         str_detect(label, 
                    "Bachelor's degree|Graduate or professional degree")) %>%
  transmute(
    var = name, 
    SEX = case_when(
      str_detect(label, "!!Male") ~ "Male",
      str_detect(label, "!!Female") ~ "Female",
      TRUE ~ NA_character_ ## TRUE is if others
    )
  ) %>% 
    filter(!is.na(SEX))

vars_needed <- c(SEX_total_vars$var, bach_vars$var)

raw <- get_acs(
    geography = "puma",
    state = "IL",
    variables = vars_needed,
    year = 2021,
    survey = "acs1",
    cache = TRUE
) %>% select(GEOID, NAME, variable, estimate)

SEX_totals <- raw %>%
  inner_join(SEX_total_vars, by = c("variable" = "var")) %>%
  group_by(GEOID, NAME, SEX) %>%
    summarise(pop_total = sum(estimate), .groups = "drop")
bach <- raw %>%
  inner_join(bach_vars, by = c("variable" = "var")) %>%
  group_by(GEOID, NAME, SEX) %>%
  summarise(pop_bach = sum(estimate), .groups = "drop")

cells <- SEX_totals %>%
  left_join(bach, by = c("GEOID","NAME","SEX")) %>%
  mutate(pop_bach = coalesce(pop_bach, 0)) %>%
  transmute(
    PUMA = substr(GEOID, 3, nchar(GEOID)),
    puma_name = NAME,
    SEX = SEX,
    pop_bach = pop_bach,
    pop_bach0 = pmax(pop_total - pop_bach, 0)
  ) %>%
  pivot_longer(
    cols = c(pop_bach, pop_bach0),
    names_to = "BACH",
    values_to = "pop"
  ) %>%
  mutate(BACH = if_else(BACH == "pop_bach", 1L, 0L),
        SEX = factor(if_else(SEX == "Male", "1", "2"), levels = c("1", "2"))) %>%
  select(PUMA, SEX, BACH, pop, puma_name) %>%
  arrange(PUMA, SEX, desc(BACH))

pcells <- pums_samp %>% 
    distinct(PUMA, SEX, BACH) %>%
    left_join(cells, by = c("PUMA", "SEX", "BACH")) %>%
    mutate(popsize = pop) %>%
    select(PUMA, SEX, BACH, popsize) %>%
    arrange(PUMA, SEX, BACH)

predX <- model.matrix(~ SEX + BACH - 1, data = pcells)
predS <- model.matrix(~ as.factor(pcells$PUMA) - 1)
pums <- pums_samp

# unit-level design
modX <- model.matrix(~ SEX + BACH - 1, data = pums)
modS <- model.matrix(~ as.factor(pums$PUMA) - 1)
modY <- pums$INCO
modZ <- pums$POV
modW <- pums$scaledWGT


# Fit models
nsim  <- 3000
nburn <- 1000
nthin <- 1

unis_wage <- unis_gaus(
  X = modX, Y = modY, S = modS,
  sig2b = 1000, wgt = modW, n = NULL,
  predX = predX, predS = predS,
  nburn = nburn, nsim = nsim, nthin = nthin,
  a = 0.5, b = 0.5, a_eps = 0.1, b_eps = 0.1
)

unis_pov <- unis_bios(
  X = modX, Y = modZ, S = modS,
  sig2b = 1000, wgt = modW, n = NULL,
  predX = predX, predS = predS,
  nburn = nburn, nsim = nsim, nthin = nthin,
  a = 0.1, b = 0.1
)

# multi-type (shared RE)
mult_br <- MTSM_br(
  X_1 = modX, X_2 = modX, Z_1 = modY, Z_2 = modZ, S = modS,
  sig2b = 1000, wgt = modW, n = NULL,
  predX = predX, predS = predS,
  nburn = nburn, nsim = nsim, nthin = nthin,
    sig2t = 1,
    sig2e = 5,
    tau_1_init = -1,
    tau_2_init = 1,
    a_eps = 0.1,
    b_eps = 0.1,
    aeta  = 0.1,
    beta  = 0.1,
    alambda = 1,
    blambda = 1
)

# Poststratify to PUMA
C <- length(unique(pcells$PUMA))
true_mean_dummy <- rep(NA_real_, C)

# Gaussian
res_ug <- gaus_post(
  preds = unis_wage$Preds,
  sig2chain = unis_wage$sig2.chain,
  true_mean = true_mean_dummy,
  region = pcells$PUMA,
  popsize = pcells$popsize
)
res_mg <- gaus_post(
  preds = mult_br$preds_gaus.chain,
  sig2chain = mult_br$sig2.chain,
  true_mean = true_mean_dummy,
  region = pcells$PUMA,
  popsize = pcells$popsize
)

# Binomial
res_ub <- bios_post(
  preds = unis_pov$Preds,
  true_mean = true_mean_dummy,
  region = pcells$PUMA,
  popsize = pcells$popsize
)
res_mb <- bios_post(
  preds = mult_br$preds_bios.chain,
  true_mean = true_mean_dummy,
  region = pcells$PUMA,
  popsize = pcells$popsize
)

# Horvitz–Thompson direct
ht_by_puma <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(
    ht_inco = stats::weighted.mean(INCO, PWGTP),
    ht_pov  = stats::weighted.mean(POV,  PWGTP),
    .groups = "drop"
  )

# Join for plotting
gaus_mapdat <- puma_sf |>
  dplyr::left_join(
    dplyr::tibble(PUMA = sort(unique(pcells$PUMA)),
                  ugaus = res_ug$est,
                  mgaus = res_mg$est),
    by = "PUMA"
  ) |>
  dplyr::left_join(ht_by_puma |> dplyr::select(PUMA, ht_inco), by = "PUMA")

bios_mapdat <- puma_sf |>
  dplyr::left_join(
    dplyr::tibble(PUMA = sort(unique(pcells$PUMA)),
                  upov = res_ub$est,
                  mpov = res_mb$est),
    by = "PUMA"
  ) |>
  dplyr::left_join(ht_by_puma |> dplyr::select(PUMA, ht_pov), by = "PUMA")

sigma2_gaus_mapdat <- puma_sf |>
  dplyr::left_join(
    dplyr::tibble(PUMA = sort(unique(pcells$PUMA)),
                  ugaus = as.numeric(res_ug$sigma2),
                  mgaus = as.numeric(res_mg$sigma2)),
    by = "PUMA"
  )
sigma2_gaus_ratio <- puma_sf |>
  left_join(
    tibble(PUMA = sort(unique(pcells$PUMA)),
                ratio = as.numeric(res_mg$sigma2/res_ug$sigma2)),
    by = "PUMA"
  )

sigma2_bios_mapdat <- puma_sf |>
  dplyr::left_join(
    dplyr::tibble(PUMA = sort(unique(pcells$PUMA)),
                  upov = as.numeric(res_ub$sigma2),
                  mpov = as.numeric(res_mb$sigma2)),
    by = "PUMA"
  )

sigma2_bios_ratio <- puma_sf |>
  left_join(
    tibble(PUMA = sort(unique(pcells$PUMA)),
           ratio = as.numeric(res_mb$sigma2/res_ub$sigma2)),
    by = "PUMA"
  )


ratio_pal <- RColorBrewer::brewer.pal(9, "Greens")
pal_seq <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
choose_pal <- rev(RColorBrewer::brewer.pal(9, "RdBu"))

robust_limits <- function(x, probs = c(0.02, 0.98)) {
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE))
}

# shared limits for three-way INCO comparison
lims_g_est <- robust_limits(c(gaus_mapdat$ugaus, gaus_mapdat$mgaus, gaus_mapdat$ht_inco))
mid_g_est  <- stats::median(c(gaus_mapdat$ugaus, gaus_mapdat$mgaus, gaus_mapdat$ht_inco), na.rm = TRUE)

# POV: keep [0,1] or use robust_limits(...) if needed
lims_b_est <- c(0, 1)

# variances (two-way)
lims_g_var <- robust_limits(
                            c(sigma2_gaus_mapdat$ugaus,
                              sigma2_gaus_mapdat$mgaus))
lims_b_var <- robust_limits(
                            c(sigma2_bios_mapdat$upov,
                              sigma2_bios_mapdat$mpov))

if (!("ht_inco" %in% names(gaus_mapdat)) || !("ht_pov" %in% names(bios_mapdat))) {
  ht_by_puma <- pums |>
    dplyr::group_by(PUMA) |>
    dplyr::summarise(
      ht_inco = stats::weighted.mean(INCO, PWGTP),
      ht_pov  = stats::weighted.mean(POV,  PWGTP),
      .groups = "drop"
    )
  gaus_mapdat <- gaus_mapdat |>
    dplyr::left_join(ht_by_puma |> dplyr::select(PUMA, ht_inco), by = "PUMA")
  bios_mapdat <- bios_mapdat |>
    dplyr::left_join(ht_by_puma |> dplyr::select(PUMA, ht_pov),  by = "PUMA")
}


ct1 <- 0.99
ct2 <- 0.99
## INCO (transformed income) — build long data + cut
inco_est_long <- gaus_mapdat |>
  dplyr::select(PUMA,
                `Horvitz–Thompson` = ht_inco,
                `Multi-type` = mgaus,
                `Univariate (Gaussian)` = ugaus,
                geometry) |>
  tidyr::pivot_longer(
    cols = c(`Horvitz–Thompson`, `Multi-type`, `Univariate (Gaussian)`),
    names_to = "source", values_to = "estimate"
  ) |>
  dplyr::mutate(source = factor(source,
    levels = c("Horvitz–Thompson", "Multi-type", "Univariate (Gaussian)")
  ))

inco_cut <- stats::quantile(inco_est_long$estimate, ct1, na.rm = TRUE)

inco_est_long <- inco_est_long |>
  dplyr::mutate(
    estimate_plot = ifelse(estimate > inco_cut, NA_real_, estimate)
  )

plot_gaus <- ggplot2::ggplot(inco_est_long) +
  ggplot2::geom_sf(ggplot2::aes(fill = estimate_plot), colour = NA) +
  ggplot2::facet_wrap(~source, nrow = 1) +
  ggplot2::scale_fill_gradientn(
    colours  = choose_pal,
    name     = "Transformed income",
    na.value = "grey90",
    limits   = c(min(inco_est_long$estimate_plot, na.rm = TRUE), inco_cut)
  ) +
  ggplot2::labs(
    title    = "Transformed Income (INCO) by PUMA",
    subtitle = "INCO = min–max scaled log income"
  ) +
  ggplot2::theme_minimal()

## POV — build long data + cut
pov_est_long <- bios_mapdat |>
  dplyr::select(PUMA,
                `Horvitz–Thompson` = ht_pov,
                `Multi-type` = mpov,
                `Univariate (Bernoulli)` = upov,
                geometry) |>
  tidyr::pivot_longer(
    cols = c(`Horvitz–Thompson`, `Multi-type`, `Univariate (Bernoulli)`),
    names_to = "source", values_to = "estimate"
  ) |>
  dplyr::mutate(source = factor(source,
    levels = c("Horvitz–Thompson", "Multi-type", "Univariate (Bernoulli)")
  ))

pov_cut <- stats::quantile(pov_est_long$estimate, ct1, na.rm = TRUE)

pov_est_long <- pov_est_long |>
  dplyr::mutate(
    estimate_plot = ifelse(estimate > pov_cut, NA_real_, estimate)
  )

plot_bios <- ggplot2::ggplot(pov_est_long) +
  ggplot2::geom_sf(ggplot2::aes(fill = estimate_plot), colour = NA) +
  ggplot2::facet_wrap(~source, nrow = 1) +
  ggplot2::scale_fill_gradientn(
    colours  = choose_pal,
    name     = "Poverty rate",
    na.value = "grey90",
    limits   = c(min(pov_est_long$estimate_plot, na.rm = TRUE), pov_cut)
  ) +
  ggplot2::labs(title = "POV by PUMA") +
  ggplot2::theme_minimal()


## INCO variance — build long data + cut
inco_var_long <- sigma2_gaus_mapdat |>
  dplyr::select(PUMA,
                `Multi-type` = mgaus,
                `Univariate (Gaussian)` = ugaus,
                geometry) |>
  tidyr::pivot_longer(
    cols = c(`Multi-type`, `Univariate (Gaussian)`),
    names_to = "source", values_to = "sigma2"
  ) |>
  dplyr::mutate(source = factor(source, 
  levels = c("Multi-type", "Univariate (Gaussian)")))

inco_var_ratio  <- sigma2_gaus_ratio |>
    select(PUMA,
            `Variance Ratio` = ratio,
            geometry)
inco_var_cut <- stats::quantile(inco_var_long$sigma2, ct2, na.rm = TRUE)

inco_var_long <- inco_var_long |>
  dplyr::mutate(
    sigma2_plot = ifelse(sigma2 > inco_var_cut, NA_real_, sigma2)
  )

plot_sigma_gaus <- ggplot2::ggplot(inco_var_long) +
  ggplot2::geom_sf(ggplot2::aes(fill = sigma2_plot), colour = NA) +
  ggplot2::facet_wrap(~source, nrow = 1) +
  ggplot2::scale_fill_gradientn(
    colours  = choose_pal,
    name     = expression(sigma^2),
    na.value = "grey90",
    limits   = c(min(inco_var_long$sigma2_plot, na.rm = TRUE), inco_var_cut)
  ) +
  ggplot2::labs(title = "Posterior Variance (Transformed Gaussian) by PUMA") +
  ggplot2::theme_minimal()

plot_sigma_gaus_ratio <- ggplot2::ggplot(inco_var_ratio) +
  ggplot2::geom_sf(ggplot2::aes(fill = `Variance Ratio`), colour = NA) +
  ggplot2::scale_fill_gradientn(
    colours  = ratio_pal,
    name     = expression(sigma^2[Multi-type] / sigma^2[Univariate]),
    na.value = "grey90",
    limits   = c(0, 2)
  ) +
  ggplot2::labs(title = "Variance Ratio (Multi-type / Univariate)") +
  ggplot2::theme_minimal()
## POV variance — build long data + cut
pov_var_long <- sigma2_bios_mapdat |>
  dplyr::select(PUMA,
                `Multi-type` = mpov,
                `Univariate (Bernoulli)` = upov,
                geometry) |>
  tidyr::pivot_longer(
    cols = c(`Multi-type`, `Univariate (Bernoulli)`),
    names_to = "source", values_to = "sigma2"
  ) |>
  dplyr::mutate(source = factor(source,
  levels = c("Multi-type", "Univariate (Bernoulli)")))

pov_var_ratio  <- sigma2_bios_ratio |>
    select(PUMA,
            `Variance Ratio` = ratio,
            geometry)

pov_var_cut <- stats::quantile(pov_var_long$sigma2, ct2, na.rm = TRUE)

pov_var_long <- pov_var_long |>
  dplyr::mutate(
    sigma2_plot = ifelse(sigma2 > pov_var_cut, NA_real_, sigma2)
  )

plot_sigma_bios <- ggplot2::ggplot(pov_var_long) +
  ggplot2::geom_sf(ggplot2::aes(fill = sigma2_plot), colour = NA) +
  ggplot2::facet_wrap(~source, nrow = 1) +
  ggplot2::scale_fill_gradientn(
    colours  = choose_pal,
    name     = expression(sigma^2),
    na.value = "grey90",
    limits   = c(min(pov_var_long$sigma2_plot, na.rm = TRUE), pov_var_cut)
  ) +
  ggplot2::labs(title = "Posterior Variance (Bernoulli) by PUMA") +
  ggplot2::theme_minimal()

plot_sigma_bios_ratio <- ggplot2::ggplot(sigma2_bios_ratio) +
  ggplot2::geom_sf(ggplot2::aes(fill = `ratio`), colour = NA) +
  ggplot2::scale_fill_gradientn(
    colours  = ratio_pal,
    name     = expression(sigma^2[Multi-type] / sigma^2[Univariate]),
    na.value = "gray",
    limits   = c(0, 1)
  ) +
  ggplot2::labs(title = "Variance Ratio (Multi-type / Univariate)") +
  ggplot2::theme_minimal()
  
  
# Render
plot_gaus; plot_bios
plot_sigma_gaus; plot_sigma_bios

plot(sigma2_bios_mapdat$mpov, sigma2_bios_mapdat$upov)
abline(0,1)
plot(sigma2_gaus_mapdat$mgau, sigma2_gaus_mapdat$ugau)
abline(0,1)

plot(mult_br$tau_1.chain, type = "l")
plot(mult_br$tau_2.chain, type = "l")

boxplot(res_mg$sigma2/res_ug$sigma2)
abline(h=1)
boxplot(res_mb$sigma2/res_ub$sigma2)
abline(h=1)
median(res_mg$sigma2/res_ug$sigma2)
median(res_mb$sigma2/res_ub$sigma2)

## spatial plot for variance
p_gaus_sigma <- plot_sigma_gaus|plot_sigma_gaus_ratio
p_bios_sigma <- plot_sigma_bios|plot_sigma_bios_ratio 

par(mfrow=c(1,2))
plot(sigma2_gaus_mapdat$mgaus, sigma2_gaus_mapdat$ugaus,
    size = 1, pch = 16,
    ylab = expression(sigma^2[Univariate]),
    xlab = expression(sigma^2[Multi-type]),
    main = "Gaussian Variance Comparison")
abline(0,1, col = "red")
plot(sigma2_bios_mapdat$mpov, sigma2_bios_mapdat$upov,
    size = 1, pch = 16,
    ylab = expression(sigma^2[Univariate]),
    xlab = expression(sigma^2[Multi-type]),
    main = "Bernoulli Variance Comparison")
abline(0,1, col = "red")
boxplot(sigma2_gaus_mapdat$mgaus/sigma2_gaus_mapdat$ugaus, 
    main = "Gaussian Variance Ratio(Multi-type/Univariate)")
abline(h=1, col = "red")
boxplot(sigma2_bios_mapdat$mpov/sigma2_bios_mapdat$upov, 
    main = "Bernoulli Variance Ratio(Multi-type/Univariate)")
abline(h=1, col = "red")

## plot the ratio map for the variance
## save the ggplots

# for (p in c("plot_gaus", "plot_bios")) {
#   ggsave(file.path("figs", paste0(p, ".png")),
#          plot = get(p), width = 7, height = 4, dpi = 300)
# }

# save.image("region_bc.rdata")
