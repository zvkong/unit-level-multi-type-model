source("packages.r")
source("functions.r")

pums_vars <- c(
  "PINCP", "SERIALNO", "SPORDER", "PWGTP", "PUMA", 
  "AGEP", "SEX", "RAC1P", "HISP", "SCHL", "HINCP", "POVPIP"
)

pums21 <- tidycensus::get_pums(
  variables = pums_vars,
  state = "IL",
  year = 2021,
  survey = "acs1",
  recode = TRUE
) |>
  dplyr::mutate(PUMA = sprintf("%05d", as.integer(PUMA)))

puma_sf <- tigris::pumas(
  state = "IL",
  year = 2021,
  class = "sf"
) |>
  sf::st_transform(4326) |>
  dplyr::mutate(PUMA = sprintf("%05d", as.integer(.data$PUMACE10)))

nb <- spdep::poly2nb(puma_sf)
W_mat <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)

eig <- eigen(W_mat)

K_basis <- sum(
  eig$values > 0
)
basis_mat <- eig$vectors[, 1:K_basis, drop = FALSE]
rownames(basis_mat) <- puma_sf$PUMA

valid_pumas <- rownames(basis_mat)

pums_samp <- pums21 |>
  dplyr::select(PUMA, PWGTP, SEX, RAC1P, POVPIP, PINCP, SCHL, AGEP) |>
  dplyr::filter(AGEP >= 18) |>
  tidyr::drop_na() |>
  dplyr::mutate(PINCP = as.numeric(PINCP)) |>
  dplyr::filter(PINCP > 0) |>
  dplyr::filter(PUMA %in% valid_pumas) |>
  dplyr::mutate(
    PWGTP = as.numeric(PWGTP),
    SEX = factor(SEX),
    LOGINCO = log(PINCP),
    BACH = dplyr::if_else(as.integer(SCHL) >= 21, 1L, 0L),
    POV = dplyr::if_else(POVPIP <= 100, 1L, 0L)
  ) |>
  dplyr::filter(is.finite(LOGINCO)) |>
  dplyr::mutate(
    INCO = (LOGINCO - min(LOGINCO)) / (max(LOGINCO) - min(LOGINCO)),
    scaledWGT = PWGTP * dplyr::n() / sum(PWGTP)
  ) |>
  dplyr::select(PUMA, SEX, BACH, INCO, POV, PWGTP, scaledWGT) |>
  dplyr::arrange(PUMA, SEX, BACH)

v <- tidycensus::load_variables(2021, "acs1", cache = TRUE) |>
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

raw <- tidycensus::get_acs(
  geography = "puma",
  state = "IL",
  variables = vars_needed,
  year = 2021,
  survey = "acs1",
  cache = TRUE
) |>
  dplyr::select(GEOID, NAME, variable, estimate)

sex_totals <- raw |>
  dplyr::inner_join(sex_total_vars, by = c("variable" = "var")) |>
  dplyr::group_by(GEOID, NAME, SEX) |>
  dplyr::summarise(pop_total = sum(estimate), .groups = "drop")

bach <- raw |>
  dplyr::inner_join(bach_vars, by = c("variable" = "var")) |>
  dplyr::group_by(GEOID, NAME, SEX) |>
  dplyr::summarise(pop_bach = sum(estimate), .groups = "drop")

cells <- sex_totals |>
  dplyr::left_join(bach, by = c("GEOID", "NAME", "SEX")) |>
  dplyr::mutate(pop_bach = dplyr::coalesce(pop_bach, 0)) |>
  dplyr::transmute(
    PUMA = sprintf("%05d", as.integer(substr(GEOID, 3, nchar(GEOID)))),
    puma_name = NAME,
    SEX = SEX,
    pop_bach = pop_bach,
    pop_bach0 = pmax(pop_total - pop_bach, 0)
  ) |>
  tidyr::pivot_longer(
    cols = c(pop_bach, pop_bach0),
    names_to = "BACH",
    values_to = "pop"
  ) |>
  dplyr::mutate(
    BACH = dplyr::if_else(BACH == "pop_bach", 1L, 0L),
    SEX = factor(dplyr::if_else(SEX == "Male", "1", "2"), levels = c("1", "2"))
  ) |>
  dplyr::select(PUMA, SEX, BACH, pop, puma_name) |>
  dplyr::arrange(PUMA, SEX, dplyr::desc(BACH))

pcells <- pums_samp |>
  dplyr::distinct(PUMA, SEX, BACH) |>
  dplyr::left_join(cells, by = c("PUMA", "SEX", "BACH")) |>
  dplyr::mutate(popsize = pop) |>
  dplyr::select(PUMA, SEX, BACH, popsize) |>
  dplyr::arrange(PUMA, SEX, BACH)

predX <- model.matrix(~ SEX + BACH - 1, data = pcells)
predS <- basis_mat[as.character(pcells$PUMA), , drop = FALSE]

pums <- pums_samp
modX <- model.matrix(~ SEX + BACH - 1, data = pums)
modS <- basis_mat[as.character(pums$PUMA), , drop = FALSE]
modY <- pums$INCO
modZ <- pums$POV
modW <- pums$scaledWGT

nsim <- 2000
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

mult_br <- MTSM_br(
  X_1 = modX, X_2 = modX, Z_1 = modY, Z_2 = modZ, S = modS,
  sig2b = 1000, wgt = modW, n = NULL,
  predX = predX, predS = predS, n_preds = NULL,
  nburn = nburn, nsim = nsim, nthin = nthin,
  sig2t = 1, sig2e = 5, tau_1_init = -1,
  a_eps = 0.1, b_eps = 0.1, aeta = 0.1, beta = 0.1,
  alambda = 1, blambda = 1
)

puma_levels <- sort(unique(pcells$PUMA))
true_mean_dummy <- stats::setNames(rep(NA_real_, length(puma_levels)), puma_levels)

res_ug <- gaus_post(
  preds = unis_wage$Preds, sig2chain = unis_wage$sig2.chain,
  true_mean = true_mean_dummy, region = pcells$PUMA, popsize = pcells$popsize
)

res_mg <- gaus_post(
  preds = mult_br$preds_gaus.chain, sig2chain = mult_br$sig2.chain,
  true_mean = true_mean_dummy, region = pcells$PUMA, popsize = pcells$popsize
)

res_ub <- bios_post(
  preds = unis_pov$Preds, true_mean = true_mean_dummy,
  region = pcells$PUMA, popsize = pcells$popsize
)

res_mb <- bios_post(
  preds = mult_br$preds_bios.chain, true_mean = true_mean_dummy,
  region = pcells$PUMA, popsize = pcells$popsize
)

ht_by_puma <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(
    ht_inco = stats::weighted.mean(INCO, PWGTP),
    ht_pov = stats::weighted.mean(POV, PWGTP),
    .groups = "drop"
  )

gaus_mapdat <- puma_sf |>
  dplyr::left_join(dplyr::tibble(PUMA = puma_levels, ugaus = res_ug$est, mgaus = res_mg$est), by = "PUMA") |>
  dplyr::left_join(ht_by_puma |> dplyr::select(PUMA, ht_inco), by = "PUMA")

bios_mapdat <- puma_sf |>
  dplyr::left_join(dplyr::tibble(PUMA = puma_levels, upov = res_ub$est, mpov = res_mb$est), by = "PUMA") |>
  dplyr::left_join(ht_by_puma |> dplyr::select(PUMA, ht_pov), by = "PUMA")

sigma2_gaus_mapdat <- puma_sf |>
  dplyr::left_join(dplyr::tibble(PUMA = puma_levels, ugaus = as.numeric(res_ug$sigma2), mgaus = as.numeric(res_mg$sigma2)), by = "PUMA")

sigma2_gaus_ratio <- puma_sf |>
  dplyr::left_join(dplyr::tibble(PUMA = puma_levels, ratio = as.numeric(res_mg$sigma2 / res_ug$sigma2)), by = "PUMA")

sigma2_bios_mapdat <- puma_sf |>
  dplyr::left_join(dplyr::tibble(PUMA = puma_levels, upov = as.numeric(res_ub$sigma2), mpov = as.numeric(res_mb$sigma2)), by = "PUMA")

sigma2_bios_ratio <- puma_sf |>
  dplyr::left_join(dplyr::tibble(PUMA = puma_levels, ratio = as.numeric(res_mb$sigma2 / res_ub$sigma2)), by = "PUMA")

save.image("region_basis_results.RData")
load("region_basis_results.RData")
ratio_pal <- RColorBrewer::brewer.pal(9, "Greens")
fill_pal <- rev(RColorBrewer::brewer.pal(9, "RdBu"))

ct_est <- 0.99
ct_var <- 0.99

inco_est_long <- gaus_mapdat |>
  dplyr::select(PUMA, `Horvitz–Thompson` = ht_inco, `Multi-type` = mgaus, `Univariate (Gaussian)` = ugaus, geometry) |>
  tidyr::pivot_longer(cols = c(`Horvitz–Thompson`, `Multi-type`, `Univariate (Gaussian)`), names_to = "source", values_to = "estimate") |>
  dplyr::mutate(source = factor(source, levels = c("Horvitz–Thompson", "Multi-type", "Univariate (Gaussian)")))

inco_cut <- stats::quantile(inco_est_long$estimate, ct_est, na.rm = TRUE)
inco_est_long <- inco_est_long |> dplyr::mutate(estimate_plot = ifelse(estimate > inco_cut, NA_real_, estimate))

plot_gaus <- ggplot2::ggplot(inco_est_long) +
  ggplot2::geom_sf(ggplot2::aes(fill = estimate_plot), colour = NA) +
  ggplot2::facet_wrap(~source, nrow = 1) +
  ggplot2::scale_fill_gradientn(colours = fill_pal, name = "Transformed income", na.value = "grey90", limits = c(min(inco_est_long$estimate_plot, na.rm = TRUE), inco_cut)) +
  ggplot2::labs(title = "Transformed Income by PUMA", subtitle = "INCO = min–max scaled log income") +
  ggplot2::theme_minimal()

pov_est_long <- bios_mapdat |>
  dplyr::select(PUMA, `Horvitz–Thompson` = ht_pov, `Multi-type` = mpov, `Univariate (Bernoulli)` = upov, geometry) |>
  tidyr::pivot_longer(cols = c(`Horvitz–Thompson`, `Multi-type`, `Univariate (Bernoulli)`), names_to = "source", values_to = "estimate") |>
  dplyr::mutate(source = factor(source, levels = c("Horvitz–Thompson", "Multi-type", "Univariate (Bernoulli)")))

pov_cut <- stats::quantile(pov_est_long$estimate, ct_est, na.rm = TRUE)
pov_est_long <- pov_est_long |> dplyr::mutate(estimate_plot = ifelse(estimate > pov_cut, NA_real_, estimate))

plot_bios <- ggplot2::ggplot(pov_est_long) +
  ggplot2::geom_sf(ggplot2::aes(fill = estimate_plot), colour = NA) +
  ggplot2::facet_wrap(~source, nrow = 1) +
  ggplot2::scale_fill_gradientn(colours = fill_pal, name = "Poverty rate", na.value = "grey90", limits = c(min(pov_est_long$estimate_plot, na.rm = TRUE), pov_cut)) +
  ggplot2::labs(title = "Poverty Rate by PUMA") +
  ggplot2::theme_minimal()

inco_var_long <- sigma2_gaus_mapdat |>
  dplyr::select(PUMA, `Multi-type` = mgaus, `Univariate (Gaussian)` = ugaus, geometry) |>
  tidyr::pivot_longer(cols = c(`Multi-type`, `Univariate (Gaussian)`), names_to = "source", values_to = "sigma2") |>
  dplyr::mutate(source = factor(source, levels = c("Multi-type", "Univariate (Gaussian)")))

inco_var_cut <- stats::quantile(inco_var_long$sigma2, ct_var, na.rm = TRUE)
inco_var_long <- inco_var_long |> dplyr::mutate(sigma2_plot = ifelse(sigma2 > inco_var_cut, NA_real_, sigma2))

plot_sigma_gaus <- ggplot2::ggplot(inco_var_long) +
  ggplot2::geom_sf(ggplot2::aes(fill = sigma2_plot), colour = NA) +
  ggplot2::facet_wrap(~source, nrow = 1) +
  ggplot2::scale_fill_gradientn(colours = fill_pal, name = expression(sigma^2), na.value = "grey90", limits = c(min(inco_var_long$sigma2_plot, na.rm = TRUE), inco_var_cut)) +
  ggplot2::labs(title = "Posterior Variance by PUMA") +
  ggplot2::theme_minimal()

plot_sigma_gaus_ratio <- ggplot2::ggplot(sigma2_gaus_ratio) +
  ggplot2::geom_sf(ggplot2::aes(fill = ratio), colour = NA) +
  ggplot2::scale_fill_gradientn(colours = ratio_pal, name = expression(sigma^2[Multi-type] / sigma^2[Univariate]), na.value = "grey90", limits = c(0, 2)) +
  ggplot2::labs(title = "Variance Ratio") +
  ggplot2::theme_minimal()

pov_var_long <- sigma2_bios_mapdat |>
  dplyr::select(PUMA, `Multi-type` = mpov, `Univariate (Bernoulli)` = upov, geometry) |>
  tidyr::pivot_longer(cols = c(`Multi-type`, `Univariate (Bernoulli)`), names_to = "source", values_to = "sigma2") |>
  dplyr::mutate(source = factor(source, levels = c("Multi-type", "Univariate (Bernoulli)")))

pov_var_cut <- stats::quantile(pov_var_long$sigma2, ct_var, na.rm = TRUE)
pov_var_long <- pov_var_long |> dplyr::mutate(sigma2_plot = ifelse(sigma2 > pov_var_cut, NA_real_, sigma2))

plot_sigma_bios <- ggplot2::ggplot(pov_var_long) +
  ggplot2::geom_sf(ggplot2::aes(fill = sigma2_plot), colour = NA) +
  ggplot2::facet_wrap(~source, nrow = 1) +
  ggplot2::scale_fill_gradientn(colours = fill_pal, name = expression(sigma^2), na.value = "grey90", limits = c(min(pov_var_long$sigma2_plot, na.rm = TRUE), pov_var_cut)) +
  ggplot2::labs(title = "Posterior Variance by PUMA") +
  ggplot2::theme_minimal()

plot_sigma_bios_ratio <- ggplot2::ggplot(sigma2_bios_ratio) +
  ggplot2::geom_sf(ggplot2::aes(fill = ratio), colour = NA) +
  ggplot2::scale_fill_gradientn(colours = ratio_pal, name = expression(sigma^2[Multi-type] / sigma^2[Univariate]), na.value = "grey90", limits = c(0, 1)) +
  ggplot2::labs(title = "Variance Ratio") +
  ggplot2::theme_minimal()

p_gaus_sigma <- plot_sigma_gaus | plot_sigma_gaus_ratio
p_bios_sigma <- plot_sigma_bios | plot_sigma_bios_ratio 

for (p in c("p_gaus_sigma", "p_bios_sigma")) {
  ggsave(file.path("figs", paste0(p, ".png")),
         plot = get(p), width = 5.5, height = 1, dpi = 500)
}

plot(sigma2_gaus_mapdat$mgaus, sigma2_gaus_mapdat$ugaus, size = 1, pch = 16, 
    ylab = expression(sigma^2[Univariate]), xlab = expression(sigma^2[Multi-type]),
    main = "Gaussian Variance Comparison")
abline(0, 1, col = 'red')

p_sig_gaus <- ggplot(sigma2_gaus_mapdat, aes(x = mgaus, y = ugaus)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  labs(x = expression(sigma^2[Multi-type]), y = expression(sigma^2[Univariate]),
       title = "Gaussian Variance Comparison") +
  theme_minimal()

p_sig_bios <- ggplot(sigma2_bios_mapdat, aes(x = mpov, y = upov)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  labs(x = expression(sigma^2[Multi-type]), y = expression(sigma^2[Univariate]),
       title = "Bernoulli Variance Comparison") +
  theme_minimal()

plot(sigma2_bios_mapdat$mpov, sigma2_bios_mapdat$upov, size = 1, pch = 16, 
    ylab = expression(sigma^2[Univariate]), xlab = expression(sigma^2[Multi-type]),
    main = "Bernoulli Variance Comparison")
abline(0, 1, col = 'red')

p_b_box <- boxplot(sigma2_gaus_mapdat$mgaus/sigma2_gaus_mapdat$ugaus, main = "Gaussian Variance Ratio(Multi-type/Univariate)")
abline(h = 1, col = 'red')

boxplot(sigma2_bios_mapdat$mpov/sigma2_bios_mapdat$upov, main = "Bernoulli Variance Ratio(Multi-type/Univariate)")
abline(h = 1, col = 'red')

par(mfrow = c(1, 2))

load("region_bc.rdata")

p_gaussian <- (plot_gaus | p_sig_gaus) + patchwork::plot_layout(widths = c(3, 1))
p_binomial <- (plot_bios | p_sig_bios) + patchwork::plot_layout(widths = c(3, 1))

for (p in c("p_gaussian", "p_binomial")) {
  ggsave(file.path("figs", paste0(p, ".png")),
         plot = get(p), width = 10.5, height = 3, dpi = 300)
}

ggplot2::ggsave(
  filename = "p_gaussian_combined.png", 
  plot = p_gaussian,                    
  width = 16,                           
  height = 5,                           
  dpi = 300                             
)

ggplot2::ggsave(
  filename = "p_binomial_combined.png", 
  plot = p_binomial,                    
  width = 16,                           
  height = 5,                           
  dpi = 300                             
)