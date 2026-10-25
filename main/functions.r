## Model samplers

# rmvn_prec: Draw a multivariate normal vector from a precision matrix and information vector.
rmvn_prec <- function(b, Q) {
  Q <- (Q + t(Q)) / 2
  L <- chol(Q)
  y <- forwardsolve(t(L), as.vector(b))
  z <- stats::rnorm(length(y))
  as.vector(backsolve(L, y + z))
}

# uni_gaus: Sample the survey-weighted Gaussian spatial model.
uni_gaus <- function(X, Y, S,
                     sig2b = 1000, wgt = NULL, n = NULL,
                     predX, predS,
                     nburn = 1000, nsim = 1000, nthin = 1,
                     a = 0.1, b = 0.1,
                     a_eps = 0.1, b_eps = 0.1, on_draw = NULL) {
  N <- length(Y)
  p <- ncol(X)
  r <- ncol(S)
  npred <- nrow(predX)

  if (is.null(wgt)) wgt <- rep(1, N)
  if (is.null(n)) n <- rep(1, N)

  X_w_X <- crossprod(X, wgt * X)
  S_w_S <- crossprod(S, wgt * S)
  X_w_S <- crossprod(X, wgt * S)
  S_w_X <- t(X_w_S)
  tX_wY <- crossprod(X, wgt * Y)
  tS_wY <- crossprod(S, wgt * Y)
  w_sum <- sum(wgt)
  Ip <- diag(p)
  Ir <- diag(r)

  Beta <- rep(1, p)
  Sigma2_u <- 1
  U <- rep(0, r)
  Mu <- rep(0, N)
  sig2 <- 1

  n_keep <- nsim / nthin
  Beta.chain <- array(0, dim = c(p, n_keep))
  Sigma2_u.chain <- numeric(n_keep)
  U.chain <- array(0, dim = c(r, n_keep))
  sig2.chain <- numeric(n_keep)
  preds.chain <- array(0, dim = c(npred, n_keep))

  pb <- utils::txtProgressBar(min = 0, max = nsim + nburn, style = 3)

  for (index in seq_len(nsim + nburn)) {
    sig2 <- 1 / stats::rgamma(
      1,
      shape = a_eps + w_sum / 2,
      rate = b_eps + 0.5 * sum(wgt * (Y - Mu)^2)
    )

    Q_Beta <- X_w_X / sig2 + Ip / sig2b
    b_Beta <- (tX_wY - X_w_S %*% U) / sig2
    Beta <- rmvn_prec(b_Beta, Q_Beta)

    Sigma2_u <- 1 / stats::rgamma(
      1,
      shape = a + r / 2,
      rate = b + 0.5 * sum(U^2)
    )

    Q_U <- S_w_S / sig2 + Ir / Sigma2_u
    b_U <- (tS_wY - S_w_X %*% Beta) / sig2
    U <- rmvn_prec(b_U, Q_U)
    U <- U - mean(U)

    Mu <- as.vector(X %*% Beta + S %*% U)
    preds <- as.vector(predX %*% Beta + predS %*% U)

    utils::setTxtProgressBar(pb, index)

    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) / nthin
      Beta.chain[, pos] <- Beta
      U.chain[, pos] <- U
      Sigma2_u.chain[pos] <- Sigma2_u
      sig2.chain[pos] <- sig2
      preds.chain[, pos] <- preds
      if (!is.null(on_draw)) on_draw(pos, gaussian = preds, sig2 = sig2)
    }
  }

  close(pb)

  list(
    Beta.chain = Beta.chain,
    U.chain = U.chain,
    Sigma2_u.chain = Sigma2_u.chain,
    sig2.chain = sig2.chain,
    Preds = preds.chain
  )
}

# uni_bio: Sample the weighted binomial spatial model with Polya-Gamma updates.
uni_bio <- function(X, Y, S,
                    sig2b = 1000, wgt = NULL, n = NULL,
                    predX, predS,
                    nburn = 1000, nsim = 5000, nthin = 1,
                    a = 0.1, b = 0.1, on_draw = NULL) {
  N <- length(Y)
  p <- ncol(X)
  r <- ncol(S)
  npred <- nrow(predX)

  if (is.null(wgt)) wgt <- rep(1, N)
  if (is.null(n)) n <- rep(1, N)

  k <- wgt * (Y - n / 2)
  tX_k <- crossprod(X, k)
  tS_k <- crossprod(S, k)
  Ip <- diag(p)
  Ir <- diag(r)

  Mu <- rep(1, N)
  Beta <- rep(1, p)
  Sigma2_u <- 1
  U <- rep(1, r)

  n_keep <- nsim / nthin
  Beta.chain <- array(0, dim = c(p, n_keep))
  Sigma2_u.chain <- numeric(n_keep)
  U.chain <- array(0, dim = c(r, n_keep))
  preds.chain <- array(0, dim = c(npred, n_keep))
  logit_bios.chain <- array(0, dim = c(npred, n_keep))

  pb <- utils::txtProgressBar(min = 0, max = nsim + nburn, style = 3)

  for (index in seq_len(nsim + nburn)) {
    omega <- BayesLogit::rpg.gamma(N, wgt * n, Mu)

    X_w_X <- crossprod(X, omega * X)
    S_w_S <- crossprod(S, omega * S)
    S_w_X <- crossprod(S, omega * X)

    Q_Beta <- X_w_X + Ip / sig2b
    b_Beta <- tX_k - t(S_w_X) %*% U
    Beta <- rmvn_prec(b_Beta, Q_Beta)

    Sigma2_u <- 1 / stats::rgamma(
      1,
      shape = a + r / 2,
      rate = b + 0.5 * sum(U^2)
    )

    Q_U <- S_w_S + Ir / Sigma2_u
    b_U <- tS_k - S_w_X %*% Beta
    U <- rmvn_prec(b_U, Q_U)
    U <- U - mean(U)

    Mu <- as.vector(X %*% Beta + S %*% U)
    logit_bios <- as.vector(predX %*% Beta + predS %*% U)
    preds <- stats::plogis(logit_bios)

    utils::setTxtProgressBar(pb, index)

    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) / nthin
      Beta.chain[, pos] <- Beta
      U.chain[, pos] <- U
      Sigma2_u.chain[pos] <- Sigma2_u
      preds.chain[, pos] <- preds
      logit_bios.chain[, pos] <- logit_bios
      if (!is.null(on_draw)) on_draw(pos, probability = preds)
    }
  }

  close(pb)

  list(
    Beta.chain = Beta.chain,
    U.chain = U.chain,
    Sigma2_u.chain = Sigma2_u.chain,
    Preds = preds.chain,
    logit_bios.chain = logit_bios.chain
  )
}

# multi: Sample the joint Gaussian-binomial model with shared area effects.
multi <- function(X_1, X_2, Z_1, Z_2, S,
                  sig2b = 1000, wgt = NULL, n = NULL,
                  predX, predS, n_preds = NULL,
                  nburn = 1000, nsim = 5000, nthin = 1,
                  sig2t = 2000, sig2e = 10, tau_2_init = 1,
                  a_eps = 0.1, b_eps = 0.1,
                  aeta = 0.1, beta = 0.1,
                  alambda = 0.1, blambda = 0.1, on_draw = NULL) {
  N <- nrow(X_1)
  r <- ncol(S)
  npred <- nrow(predX)

  if (is.null(wgt)) wgt <- rep(1, N)
  if (is.null(n)) n <- rep(1, N)

  w_sum <- sum(wgt)
  X1_w_X1 <- crossprod(X_1, wgt * X_1)
  S_w_S <- crossprod(S, wgt * S)
  X1_w_S <- crossprod(X_1, wgt * S)
  S_w_X1 <- t(X1_w_S)
  tX1_wZ1 <- crossprod(X_1, wgt * Z_1)
  tS_wZ1 <- crossprod(S, wgt * Z_1)

  k <- wgt * (Z_2 - n / 2)
  tX2_k <- crossprod(X_2, k)
  tS_k <- crossprod(S, k)

  p_1 <- ncol(X_1)
  p_2 <- ncol(X_2)
  Ip1 <- diag(p_1)
  Ip2 <- diag(p_2)
  Ir <- diag(r)

  tau_2 <- tau_2_init
  Beta_1 <- rep(1, p_1)
  Beta_2 <- rep(1, p_2)
  lambda <- rep(0, r)
  eta <- rep(1, r)
  Mu_1 <- rep(1, N)
  Mu_2 <- rep(1, N)
  sig2 <- 1
  sig2l <- 1

  n_keep <- nsim / nthin
  tau_2.chain <- numeric(n_keep)
  Beta_1.chain <- array(0, dim = c(p_1, n_keep))
  Beta_2.chain <- array(0, dim = c(p_2, n_keep))
  Sigma2_lambda.chain <- numeric(n_keep)
  Sigma2_eta.chain <- numeric(n_keep)
  sig2.chain <- numeric(n_keep)
  lambda.chain <- array(0, dim = c(r, n_keep))
  eta.chain <- array(0, dim = c(r, n_keep))
  preds_gaus.chain <- array(0, dim = c(npred, n_keep))
  preds_bios.chain <- array(0, dim = c(npred, n_keep))
  logit_bios.chain <- array(0, dim = c(npred, n_keep))

  pb <- utils::txtProgressBar(min = 0, max = nsim + nburn, style = 3)

  for (index in seq_len(nsim + nburn)) {
    sig2 <- 1 / stats::rgamma(
      1,
      shape = a_eps + w_sum / 2,
      rate = b_eps + 0.5 * sum(wgt * (Z_1 - Mu_1)^2)
    )

    omega <- BayesLogit::rpg.gamma(N, wgt * n, Mu_2)
    OM <- crossprod(S, omega * S)
    X2_w_X2 <- crossprod(X_2, omega * X_2)
    S_w_X2 <- crossprod(S, omega * X_2)

    sig2e <- 1 / stats::rgamma(
      1,
      shape = aeta + r / 2,
      rate = beta + 0.5 * sum(eta^2)
    )

    Q_eta <- tau_2^2 * (S_w_S / sig2) + Ir / sig2e + OM
    b_eta <- tau_2 * (tS_wZ1 - S_w_X1 %*% Beta_1) / sig2 +
      (tS_k - S_w_X2 %*% Beta_2 - OM %*% lambda)
    eta <- rmvn_prec(b_eta, Q_eta)
    eta <- eta - mean(eta)

    Q_Beta_1 <- (X1_w_X1 / sig2) + Ip1 / sig2b
    b_Beta_1 <- (tX1_wZ1 - tau_2 * X1_w_S %*% eta) / sig2
    Beta_1 <- rmvn_prec(b_Beta_1, Q_Beta_1)

    sig2l <- 1 / stats::rgamma(
      1,
      shape = alambda + r / 2,
      rate = blambda + 0.5 * sum(lambda^2)
    )
    Q_lambda <- OM + Ir / sig2l
    b_lambda <- tS_k - S_w_X2 %*% Beta_2 - OM %*% eta
    lambda <- rmvn_prec(b_lambda, Q_lambda)
    lambda <- lambda - mean(lambda)

    Q_Beta_2 <- X2_w_X2 + Ip2 / sig2b
    b_Beta_2 <- tX2_k - t(S_w_X2) %*% eta - t(S_w_X2) %*% lambda
    Beta_2 <- rmvn_prec(b_Beta_2, Q_Beta_2)

    M_sig2 <- S_w_S / sig2
    var_tau_2 <- 1 / (sum(eta * (M_sig2 %*% eta)) + 1 / sig2t)
    mean_tau_2 <- var_tau_2 * sum(eta * ((tS_wZ1 - S_w_X1 %*% Beta_1) / sig2))
    tau_2 <- stats::rnorm(1, mean_tau_2, sqrt(var_tau_2))

    Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_2 * S %*% eta)
    Mu_2 <- as.vector(X_2 %*% Beta_2 + S %*% eta + S %*% lambda)

    preds_gaus <- as.vector(predX %*% Beta_1 + tau_2 * predS %*% eta)
    logit_bios <- as.vector(predX %*% Beta_2 + predS %*% eta + predS %*% lambda)
    preds_bios <- stats::plogis(logit_bios)

    utils::setTxtProgressBar(pb, index)

    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) / nthin
      tau_2.chain[pos] <- tau_2
      Beta_1.chain[, pos] <- Beta_1
      Beta_2.chain[, pos] <- Beta_2
      lambda.chain[, pos] <- lambda
      eta.chain[, pos] <- eta
      sig2.chain[pos] <- sig2
      Sigma2_lambda.chain[pos] <- sig2l
      Sigma2_eta.chain[pos] <- sig2e
      preds_gaus.chain[, pos] <- preds_gaus
      preds_bios.chain[, pos] <- preds_bios
      logit_bios.chain[, pos] <- logit_bios
      if (!is.null(on_draw))
        on_draw(pos, gaussian = preds_gaus, sig2 = sig2, probability = preds_bios)
    }
  }

  close(pb)

  list(
    Beta_1.chain = Beta_1.chain,
    Beta_2.chain = Beta_2.chain,
    lambda.chain = lambda.chain,
    eta.chain = eta.chain,
    Sigma2_lambda.chain = Sigma2_lambda.chain,
    Sigma2_eta.chain = Sigma2_eta.chain,
    sig2.chain = sig2.chain,
    preds_gaus.chain = preds_gaus.chain,
    preds_bios.chain = preds_bios.chain,
    logit_bios.chain = logit_bios.chain,
    tau_2.chain = tau_2.chain
  )
}

## Transformation and posterior summaries

# boxcox_transform: Transform positive incomes to the Box-Cox scale.
boxcox_transform <- function(y, lambda) {
  if (any(!is.finite(y) | y <= 0)) stop("Box-Cox income must be strictly positive")
  if (abs(lambda) < 1e-8) log(y) else (y^lambda - 1) / lambda
}

# sample_boxcox_lambda: Estimate sample lambda using a survey-weighted working-regression profile.
sample_boxcox_lambda <- function(y, x, weights, grid = seq(-2, 2, by = 0.01)) {
  y <- as.numeric(y)
  weights <- as.numeric(weights)
  if (length(y) < 2L || length(weights) != length(y) || nrow(x) != length(y) ||
      any(!is.finite(y) | y <= 0) || any(!is.finite(weights) | weights <= 0))
    stop("Lambda estimation requires positive sample incomes and inverse-inclusion weights")
  if (length(grid) < 3L || any(!is.finite(grid)) || is.unsorted(grid, strictly = TRUE))
    stop("Lambda grid must be finite and strictly increasing")
  design <- model.matrix(~ ., as.data.frame(x))
  if (nrow(design) != length(y) || any(!is.finite(design)))
    stop("Incomplete sample covariates for lambda estimation")
  n <- length(y)
  w <- weights / max(weights)
  w <- w * (n / sum(w))
  root_w <- sqrt(w)
  decomposition <- qr(design * root_w)
  if (decomposition$rank != ncol(design) || n <= ncol(design))
    stop("The sample lambda working regression is not identifiable")
  log_y <- log(y)
  jacobian_sum <- sum(w * log_y)
  center <- jacobian_sum / n
  # Profile RSS: Evaluate the centered weighted residual sum for one lambda.
  scaled_rss <- vapply(grid, function(lambda) {
    z <- if (lambda == 0) log_y - center else expm1(lambda * (log_y - center)) / lambda
    if (any(!is.finite(z))) return(Inf)
    sum(qr.resid(decomposition, z * root_w)^2)
  }, numeric(1))
  log_rss <- log(scaled_rss) + 2 * grid * center
  profile <- -n / 2 * (log_rss - log(n)) + (grid - 1) * jacobian_sum
  valid <- is.finite(profile) & is.finite(scaled_rss) & scaled_rss > 0
  if (!any(valid)) stop("No finite sample Box-Cox profile values")
  profile[!valid] <- -Inf
  index <- which.max(profile)
  if (index %in% c(1L, length(grid)))
    stop("Sample lambda optimum is at the search boundary; review the grid explicitly")
  list(lambda = unname(grid[index]), index = index, grid = grid, profile = profile,
       sample_size = n, design_columns = colnames(design),
       method = "Survey-weighted Gaussian working-regression profile with weighted Jacobian",
       weight_source = "1 / actual inclusion probability, normalized to sample size",
       source = "current sample only")
}

# check_prediction_lambda: Check the original-income predictive second-moment restriction.
check_prediction_lambda <- function(lambda) {
  if (length(lambda) != 1L || !is.finite(lambda)) stop("Invalid prediction lambda")
  if (lambda < 0 && lambda >= -2)
    stop("Estimated negative lambda gives an infinite second income moment under the existing Gaussian inverse-domain prediction rule; no fit was run")
  invisible(lambda)
}

# poststratify_cells: Aggregate cell proportion draws using population-share weights.
poststratify_cells <- function(cell_draws, area_weights, areas) {
  result <- as.matrix(area_weights %*% cell_draws)
  rownames(result) <- areas
  colnames(result) <- seq_len(ncol(result))
  result
}

# summarize_draws: Calculate posterior means, 95% intervals and variances.
summarize_draws <- function(post) {
  stopifnot(is.matrix(post), ncol(post) >= 2L, all(is.finite(post)),
            !is.null(rownames(post)))
  data.frame(PUMA = rownames(post), estimate = rowMeans(post),
             lower = apply(post, 1, quantile, 0.025),
             upper = apply(post, 1, quantile, 0.975), variance = apply(post, 1, var),
             row.names = NULL)
}

# score_intervals: Calculate squared error, interval score and coverage.
score_intervals <- function(estimate, lower, upper, truth) {
  if (any(!is.finite(c(estimate, lower, upper, truth))) || any(lower > upper))
    stop("Invalid estimate, truth or interval; no records are silently excluded")
  data.frame(squared_error = (estimate - truth)^2,
             interval_score = upper - lower + 40 * pmax(lower - truth, 0) +
               40 * pmax(truth - upper, 0),
             covered = as.numeric(lower <= truth & truth <= upper))
}

# emp_poverty_direct_ht: Estimate poverty from the HT total divided by known population size.
emp_poverty_direct_ht <- function(samp, areas, population_sizes) {
  N <- population_sizes[match(as.character(areas), names(population_sizes))]
  if (any(!is.finite(N) | N <= 0)) stop("Missing known PUMA population sizes")
  design <- survey::svydesign(ids = ~1, weights = ~W, data = samp)
  total <- survey::svyby(~POV, ~PUMA, design, survey::svytotal,
                         vartype = "se", keep.names = FALSE)
  ix <- match(as.character(areas), as.character(total$PUMA))
  estimate <- as.numeric(total$POV[ix] / N)
  se <- as.numeric(survey::SE(total)[ix] / N)
  lower <- pmax(0, estimate - qnorm(0.975) * se)
  upper <- pmin(1, estimate + qnorm(0.975) * se)
  if (any(!is.finite(c(estimate, lower, upper, se))) || any(lower > upper))
    stop("Incomplete HT Direct estimates or intervals")
  data.frame(PUMA = areas, response = "poverty", model = "Direct",
    estimate = estimate, lower = lower, upper = upper, variance = se^2)
}

# predict_cell_poverty: Draw finite-population poverty proportions for each cell.
predict_cell_poverty <- function(probability, popsize) {
  draws <- matrix(0, nrow(probability), ncol(probability))
  for (t in seq_len(ncol(probability))) {
    draws[, t] <- rbinom(length(popsize), popsize, probability[, t]) / pmax(popsize, 1)
  }
  draws
}

# predict_poverty_region: Draw cell counts and aggregate regional poverty proportions.
predict_poverty_region <- function(probability, popsize, area_weights) {
  N_int <- as.integer(round(popsize))
  counts <- matrix(rbinom(length(N_int) * ncol(probability),
    size = rep(N_int, times = ncol(probability)), prob = as.vector(probability)),
    nrow = length(N_int), ncol = ncol(probability))
  draws <- as.matrix(area_weights %*% counts)
  rownames(draws) <- rownames(area_weights)
  draws
}

# md_weighted_quantile: Invert the weighted sample CDF with midpoint handling at weight boundaries.
md_weighted_quantile <- function(x, w, p) {
  take <- which(w > 0)
  take <- take[order(x[take])]
  x <- x[take]
  w <- w[take] / max(w[take])
  if (p == 0) return(x[1L])
  if (p == 1) return(tail(x, 1L))
  total <- sum(w)
  target <- p * total
  cumulative <- cumsum(w)
  tolerance <- 8 * .Machine$double.eps * total
  exact <- which(abs(cumulative - target) <= tolerance)
  if (length(exact)) {
    k <- exact[1L]
    return(if (k < length(x)) (x[k] + x[k + 1L]) / 2 else x[k])
  }
  x[which(cumulative >= target)[1L]]
}

# md_direct: Estimate weighted medians and approximate Woodruff intervals.
md_direct <- function(samp, areas) {
  # Area median: Calculate the current sample weighted median for one PUMA.
  estimate <- vapply(areas, function(area) {
    use <- as.character(samp$PUMA) == area
    md_weighted_quantile(samp$PINCP_NUM[use], samp$W[use], 0.5)
  }, numeric(1))
  # Woodruff: form an interval for the median CDF, then invert the sample CDF.
  samp$median_cdf <- as.numeric(samp$PINCP_NUM <= estimate[match(as.character(samp$PUMA), areas)])
  design <- survey::svydesign(ids = ~1, weights = ~W, data = samp)
  value <- survey::svyby(~median_cdf, ~PUMA, design, survey::svymean,
    vartype = "se", keep.names = FALSE)
  ix <- match(areas, as.character(value$PUMA))
  probability <- as.numeric(value$median_cdf[ix])
  se_probability <- as.numeric(survey::SE(value)[ix])
  critical <- stats::qnorm(0.975)
  raw_lower <- probability - critical * se_probability
  raw_upper <- probability + critical * se_probability
  lower_p <- pmax(0, raw_lower)
  upper_p <- pmin(1, raw_upper)
  # Area interval: Invert the two CDF confidence limits for one PUMA.
  limits <- vapply(seq_along(areas), function(j) {
    use <- as.character(samp$PUMA) == areas[j]
    # CDF inversion: Calculate the weighted quantile at one probability.
    vapply(c(lower_p[j], upper_p[j]), function(p)
      md_weighted_quantile(samp$PINCP_NUM[use], samp$W[use], p), numeric(1))
  }, numeric(2))
  result <- data.frame(PUMA = areas, response = "income", model = "Direct",
    estimate = as.numeric(estimate), lower = limits[1L, ], upper = limits[2L, ],
    variance = ((limits[2L, ] - limits[1L, ]) / (2 * critical))^2, row.names = NULL)
  if (any(!is.finite(as.matrix(result[c("estimate", "lower", "upper", "variance")]))) ||
      any(result$lower > result$upper)) stop("Undefined Direct median intervals")
  attr(result, "diagnostics") <- data.frame(PUMA = areas,
    cdf_at_median = probability, cdf_se = se_probability,
    clipped_lower_probability = raw_lower < 0, clipped_upper_probability = raw_upper > 1,
    zero_width = result$lower == result$upper)
  result
}

# pure_r_empirical_medians: Predict all-unit income medians and optional joint poverty draws.
pure_r_empirical_medians <- function(mu, sig2, lambda, popsize, cell_area,
                                     chunk_units = 100000L,
                                     probability = NULL) {
  cells <- nrow(mu)
  draws <- ncol(mu)
  stopifnot(length(popsize) == cells, length(cell_area) == cells,
    length(sig2) == draws, all(is.finite(sig2) & sig2 >= 0))
  population <- as.numeric(popsize)
  areas <- max(cell_area)
  area_population <- numeric(areas)
  for (j in seq_len(cells)) area_population[cell_area[j]] <-
    area_population[cell_area[j]] + population[j]
  log_case <- abs(lambda) < 1e-8
  median_draws <- matrix(NA_real_, areas, draws,
    dimnames = list(as.character(seq_len(areas)), as.character(seq_len(draws))))
  poverty <- if (is.null(probability)) NULL else matrix(0, cells, draws)
  if (!is.null(colnames(mu))) colnames(median_draws) <- colnames(mu)
  audit <- matrix(0, cells, 6L, dimnames = list(rownames(mu),
    c("proposals", "invalid", "affected_units", "generated_units",
      "fixed_observed_units", "total_population_units")))
  audit[, "total_population_units"] <- population * draws
  for (t in seq_len(draws)) {
    sd <- sqrt(sig2[t])
    area_values <- lapply(area_population, numeric)
    filled <- numeric(areas)
    for (j in seq_len(cells)) {
      remaining <- population[j]
      mean <- mu[j, t]
      if (sig2[t] == 0 && !log_case && 1 + lambda * mean <= 0)
        stop("Degenerate predictive distribution outside inverse domain")
      d <- cell_area[j]
      while (remaining > 0) {
        size <- min(chunk_units, remaining)
        z <- stats::rnorm(size, mean, sd)
        bad <- rep(FALSE, size)
        if (!log_case) bad <- 1 + lambda * z <= 0
        audit[j, "proposals"] <- audit[j, "proposals"] + size
        audit[j, "affected_units"] <- audit[j, "affected_units"] + sum(bad)
        retries <- 0L
        while (any(bad)) {
          retries <- retries + 1L
          if (retries > 10000L) stop("Excessive rejection rate; inspect the fitted iteration")
          rejected <- sum(bad)
          audit[j, "invalid"] <- audit[j, "invalid"] + rejected
          audit[j, "proposals"] <- audit[j, "proposals"] + rejected
          z[bad] <- stats::rnorm(rejected, mean, sd)
          bad <- 1 + lambda * z <= 0
        }
        if (any(!is.finite(z))) stop("Non-finite transformed prediction")
        if (log_case) y <- exp(z) else y <- exp(log1p(lambda * z) / lambda)
        if (any(!is.finite(y) | y <= 0)) stop("Invalid original-dollar prediction")
        area_values[[d]][filled[d] + seq_len(size)] <- y
        filled[d] <- filled[d] + size
        audit[j, "generated_units"] <- audit[j, "generated_units"] + size
        remaining <- remaining - size
      }
    }
    if (any(filled != area_population)) stop("Predicted population counts disagree")
    if (!is.null(probability)) poverty[, t] <-
      rbinom(cells, popsize, probability[, t]) / pmax(popsize, 1)
    # median() averages the two central original-dollar values when N is even.
    for (d in seq_len(areas)) median_draws[d, t] <- stats::median(area_values[[d]])
  }
  list(median = median_draws, poverty = poverty, audit = audit)
}

# pure_r_region_medians: Predict population medians using conditional order statistics.
pure_r_region_medians <- function(mu, sig2, lambda, popsize, cell_area) {
  # inverse: Transform one predicted income back to original dollars.
  inverse <- function(z) {
    if (abs(lambda) < 1e-8) y <- exp(z) else y <- exp(log1p(lambda * z) / lambda)
    if (!is.finite(y) || y <= 0) stop("Nonfinite inverse Box-Cox prediction")
    y
  }
  # normal_mass: Calculate normal interval probabilities using stable tail calculations.
  normal_mass <- function(a, b) {
    if (a >= b) return(0)
    if (a >= 0) {
      pnorm(a, 0, 1, lower.tail = FALSE, log.p = FALSE) -
        pnorm(b, 0, 1, lower.tail = FALSE, log.p = FALSE)
    } else {
      pnorm(b, 0, 1, lower.tail = TRUE, log.p = FALSE) -
        pnorm(a, 0, 1, lower.tail = TRUE, log.p = FALSE)
    }
  }
  # conditional_normal: Draw a normal value conditional on interval bounds.
  conditional_normal <- function(mu, sd, lo, hi) {
    a <- (lo - mu) / sd
    b <- (hi - mu) / sd
    lower <- a < 0
    pa <- pnorm(a, 0, 1, lower.tail = lower, log.p = FALSE)
    pb <- pnorm(b, 0, 1, lower.tail = lower, log.p = FALSE)
    u <- runif(1L, 0, 1)
    p <- pa + u * (pb - pa)
    if (p <= 0 || p >= 1 || pa == pb)
      stop("Numerically unresolved conditional normal interval")
    mu + sd * qnorm(p, 0, 1, lower.tail = lower, log.p = FALSE)
  }

  cells <- nrow(mu)
  draws <- ncol(mu)
  stopifnot(length(popsize) == cells, length(cell_area) == cells,
    length(sig2) == draws,
    all(is.finite(popsize) & popsize >= 0 & popsize == floor(popsize)),
    all(popsize <= .Machine$integer.max), all(is.finite(sig2) & sig2 >= 0))
  pop <- as.integer(popsize)
  area <- as.integer(cell_area) - 1L
  areas <- max(area) + 1L
  indices <- area_counts <- ranks <- vector("list", areas)
  for (d in seq_len(areas)) indices[[d]] <- integer(0)
  for (j in seq_len(cells)) indices[[area[j] + 1L]] <- c(indices[[area[j] + 1L]], j)
  lo <- if (abs(lambda) < 1e-8 || lambda < 0) -Inf else -1 / lambda
  hi <- if (abs(lambda) < 1e-8 || lambda > 0) Inf else -1 / lambda
  for (d in seq_len(areas)) {
    total <- 0
    area_counts[[d]] <- integer(length(indices[[d]]))
    for (i in seq_along(indices[[d]])) {
      j <- indices[[d]][i]
      total <- total + pop[j]
      area_counts[[d]][i] <- pop[j]
    }
    if (total < 1 || total > .Machine$integer.max) stop("Invalid area population")
    ranks[[d]] <- integer(0)
    if (total %% 2 == 0) ranks[[d]] <- as.integer(total / 2)
    ranks[[d]] <- c(ranks[[d]], as.integer(floor(total / 2) + 1))
  }
  medians <- matrix(0, areas, draws)
  audit <- matrix(0, draws, 3L,
    dimnames = list(NULL, c("partition_nodes", "leaf_generated_units", "max_depth")))
  for (t in seq_len(draws)) {
    for (d in seq_len(areas)) {
      sd <- sqrt(sig2[t])
      means <- numeric(length(indices[[d]]))
      nodes <- generated <- 0
      max_depth <- 0L
      for (i in seq_along(indices[[d]])) {
        j <- indices[[d]][i]
        if (!is.finite(mu[j, t])) stop("Nonfinite fitted mean")
        if (sig2[t] == 0 && (mu[j, t] <= lo || mu[j, t] >= hi))
          stop("Degenerate Gaussian outside inverse domain")
        means[i] <- mu[j, t]
      }

      # select: Select central order statistics with shared recursive partitions.
      select <- function(n, lower, upper, requested, depth = 0L) {
        nodes <<- nodes + 1
        max_depth <<- max(max_depth, depth)
        if (depth > 120L) stop("Order-statistic recursion did not resolve")
        latent <- 0L
        for (value in n) latent <- latent + value
        for (k in requested) if (k < 1L || k > latent) stop("Invalid conditional rank")
        if (sd == 0) {
          mass_values <- numeric(0)
          mass_counts <- integer(0)
          for (j in seq_along(n)) if (n[j] != 0L) {
            mass_values <- c(mass_values, inverse(means[j]))
            mass_counts <- c(mass_counts, n[j])
          }
          ordered <- order(mass_values, mass_counts)
          mass_values <- mass_values[ordered]
          mass_counts <- mass_counts[ordered]
          answer <- numeric(0)
          for (k in requested) {
            accumulated <- 0L
            for (i in seq_along(mass_values)) {
              accumulated <- accumulated + mass_counts[i]
              if (accumulated >= k) {
                answer <- c(answer, mass_values[i])
                break
              }
            }
          }
          return(answer)
        }
        if (latent <= 32L) {
          values <- numeric(0)
          for (j in seq_along(n)) if (n[j] > 0L) for (i in seq_len(n[j]))
            values <- c(values, inverse(conditional_normal(means[j], sd, lower, upper)))
          generated <<- generated + latent
          values <- sort(values)
          return(values[requested])
        }
        if (is.finite(lower) && is.finite(upper)) {
          pivot <- lower + (upper - lower) / 2
        } else {
          center <- 0
          for (j in seq_along(n)) center <- center + n[j] * means[j] / latent
          pivot <- center
          if (is.finite(lower) && pivot <= lower) pivot <- lower + sd
          if (is.finite(upper) && pivot >= upper) pivot <- upper - sd
        }
        if (!(pivot > lower && pivot < upper) || !is.finite(pivot))
          stop("Numerically unresolved partition boundaries")
        left <- right <- integer(length(n))
        for (j in seq_along(n)) {
          if (n[j] == 0L) next
          a <- (lower - means[j]) / sd
          b <- (upper - means[j]) / sd
          mass <- normal_mass(a, b)
          if (!(mass > 0)) stop("Conditional interval has zero numerical mass")
          probability <- normal_mass(a, (pivot - means[j]) / sd) / mass
          if (!is.finite(probability) || probability < -1e-12 || probability > 1 + 1e-12)
            stop("Invalid binomial partition probability")
          probability <- max(0, min(1, probability))
          left[j] <- as.integer(rbinom(1L, n[j], probability))
          right[j] <- n[j] - left[j]
        }
        n_left <- 0L
        for (value in left) n_left <- n_left + value
        left_ranks <- right_ranks <- integer(0)
        for (k in requested) {
          if (k <= n_left) left_ranks <- c(left_ranks, k) else
            right_ranks <- c(right_ranks, k - n_left)
        }
        answer <- numeric(0)
        if (length(left_ranks))
          answer <- select(left, lower, pivot, left_ranks, depth + 1L)
        if (length(right_ranks)) {
          other <- select(right, pivot, upper, right_ranks, depth + 1L)
          answer <- c(answer, other)
        }
        answer
      }
      values <- select(area_counts[[d]], lo, hi, ranks[[d]])
      medians[d, t] <- if (length(values) == 1L) values[1L] else
        values[1L] / 2 + values[2L] / 2
      audit[t, 1L] <- audit[t, 1L] + nodes
      audit[t, 2L] <- audit[t, 2L] + generated
      audit[t, 3L] <- max(audit[t, 3L], max_depth)
    }
  }
  list(medians = medians, audit = audit)
}
