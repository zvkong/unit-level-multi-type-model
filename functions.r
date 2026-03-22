large_abs <- function(v, q) {
  if (!is.numeric(q) || length(q) != 1L || q <= 0 || q > 1) {
    stop("q must be a single number in (0, 1].")
  }
  n_keep <- max(1L, ceiling(q * length(v)))
  order(abs(v), decreasing = TRUE)[seq_len(n_keep)]
}

rmvn_prec <- function(b, Q) {
  Q <- (Q + t(Q)) / 2 
  L <- chol(Q)
  y <- forwardsolve(t(L), as.vector(b))
  z <- stats::rnorm(length(y))
  as.vector(backsolve(L, y + z))
}

interval_score <- function(lower, upper, truth, alpha = 0.05) {
  width        <- upper - lower
  penalty_low  <- (2 / alpha) * pmax(0, lower - truth)
  penalty_high <- (2 / alpha) * pmax(0, truth - upper)
  width + penalty_low + penalty_high
}

gaus_post <- function(preds, sig2chain, true_mean, region, popsize) {
  region  <- as.character(region)
  regions <- unique(region)
  C       <- length(regions)
  nsim    <- ncol(preds)

  idx_by_reg <- split(seq_along(region), region)
  N_by_reg   <- sapply(idx_by_reg[regions], function(id) sum(popsize[id]))

  row_idx <- integer(length(region))
  for(i in seq_along(regions)) { row_idx[idx_by_reg[[regions[i]]]] <- i }
  
  W_post <- Matrix::sparseMatrix(
    i = row_idx, j = seq_along(region), x = popsize / N_by_reg[row_idx], dims = c(C, length(region))
  )

  post_mean <- as.matrix(W_post %*% preds)
  rownames(post_mean) <- regions

  mu_sd <- sqrt(outer(1 / N_by_reg, sig2chain))
  
  post <- matrix(stats::rnorm(C * nsim, mean = as.vector(post_mean), sd = as.vector(mu_sd)), C, nsim)
  rownames(post) <- regions

  if (!is.null(names(true_mean))) true_mean <- true_mean[regions]

  est    <- rowMeans(post)
  sigma2 <- apply(post, 1, stats::var)
  lb     <- apply(post, 1, stats::quantile, probs = 0.025)
  ub     <- apply(post, 1, stats::quantile, probs = 0.975)
  cr     <- mean(lb <= true_mean & true_mean <= ub)

  est_param    <- rowMeans(post_mean)
  var_param    <- apply(post_mean, 1, stats::var)
  lb_param     <- apply(post_mean, 1, stats::quantile, probs = 0.025)
  ub_param     <- apply(post_mean, 1, stats::quantile, probs = 0.975)
  cr_param     <- mean(lb_param <= true_mean & true_mean <= ub_param)

  list(est = est, lb = lb, ub = ub, cr = cr, sigma2 = sigma2, post = post,
       est_nonoise = est_param, lb_nonoise = lb_param, ub_nonoise = ub_param,
       cr_nonoise = cr_param, sigma2_nonoise = var_param, post_nonoise = post_mean)
}

bios_post <- function(preds, true_mean, region, popsize) {
  region  <- as.character(region)
  regions <- unique(region)
  C       <- length(regions)
  R_sim   <- ncol(preds)

  idx_by_reg <- split(seq_along(region), region)
  N_by_reg   <- sapply(idx_by_reg[regions], function(id) sum(popsize[id]))

  row_idx <- integer(length(region))
  for(i in seq_along(regions)) { row_idx[idx_by_reg[[regions[i]]]] <- i }
  
  W_post <- Matrix::sparseMatrix(
    i = row_idx, j = seq_along(region), x = popsize / N_by_reg[row_idx], dims = c(C, length(region))
  )
  
  post_mean <- as.matrix(W_post %*% preds)
  rownames(post_mean) <- regions

  N_int <- as.integer(round(popsize))
  s_matrix <- matrix(
    stats::rbinom(n = length(N_int) * R_sim, size = rep(N_int, times = R_sim), prob = as.vector(preds)), 
    nrow = length(N_int), ncol = R_sim
  )

  W_sum <- Matrix::sparseMatrix(
    i = row_idx, j = seq_along(region), x = 1 / N_by_reg[row_idx], dims = c(C, length(region))
  )
  
  post <- as.matrix(W_sum %*% s_matrix)
  rownames(post) <- regions

  if (!is.null(names(true_mean))) true_mean <- true_mean[regions]

  est    <- rowMeans(post)
  sigma2 <- apply(post, 1, stats::var)
  lb     <- apply(post, 1, stats::quantile, probs = 0.025)
  ub     <- apply(post, 1, stats::quantile, probs = 0.975)
  cr     <- mean(lb <= true_mean & true_mean <= ub)
  
  est_param <- rowMeans(post_mean)
  var_param <- apply(post_mean, 1, stats::var)
  lb_param  <- apply(post_mean, 1, stats::quantile, probs = 0.025)
  ub_param  <- apply(post_mean, 1, stats::quantile, probs = 0.975)
  cr_param  <- mean(lb_param <= true_mean & true_mean <= ub_param)

  list(est = est, lb = lb, ub = ub, cr = cr, sigma2 = sigma2, post = post,
       est_nonoise = est_param, lb_nonoise = lb_param, ub_nonoise = ub_param,
       cr_nonoise = cr_param, sigma2_nonoise = var_param, post_nonoise = post_mean)
}

unis_bios <- function(X, Y, S, sig2b = 1000, wgt = NULL, n = NULL,
                      predX, predS, nburn = 1000, nsim = 5000, nthin = 1,
                      a = 0.1, b = 0.1) {
  N <- length(Y); p <- ncol(X); r <- ncol(S); npred <- nrow(predX)
  if (is.null(wgt)) wgt <- rep(1, N)
  if (is.null(n))   n   <- rep(1, N)

  k <- wgt * (Y - n / 2)
  tX_k <- crossprod(X, k)
  tS_k <- crossprod(S, k)
  Ip <- diag(p); Ir <- diag(r)

  Mu <- rep(1, N); Beta <- rep(1, p); Sigma2_u <- 1; U <- rep(1, r)
  n_keep <- nsim / nthin
  Beta.chain <- array(0, dim = c(p, n_keep))
  Sigma2_u.chain <- numeric(n_keep)
  U.chain <- array(0, dim = c(r, n_keep))
  Mu.chain <- array(0, dim = c(N, n_keep))
  preds.chain <- array(0, dim = c(npred, n_keep))
  logit_bios.chain <- array(0, dim = c(npred, n_keep))

  pb <- txtProgressBar(min = 0, max = nsim + nburn, style = 3)

  for (index in seq_len(nsim + nburn)) {
    omega <- BayesLogit::rpg.gamma(N, wgt, Mu)
    
    X_w_X <- crossprod(X, omega * X)
    S_w_S <- crossprod(S, omega * S)
    S_w_X <- crossprod(S, omega * X)

    Q_Beta  <- X_w_X + Ip / sig2b
    b_Beta  <- tX_k - t(S_w_X) %*% U
    Beta    <- rmvn_prec(b_Beta, Q_Beta)

    b_prime   <- b + 0.5 * sum(U^2)
    Sigma2_u  <- 1 / stats::rgamma(1, shape = a + r / 2, rate = b_prime)

    Q_U  <- S_w_S + Ir / Sigma2_u
    b_U  <- tS_k - S_w_X %*% Beta
    U    <- rmvn_prec(b_U, Q_U)

    Mu <- as.vector(X %*% Beta + S %*% U)
    logit_bios <- predX %*% Beta + predS %*% U
    preds      <- stats::plogis(logit_bios)

    utils::setTxtProgressBar(pb, index)
    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) / nthin
      Beta.chain[, pos] <- Beta; U.chain[, pos] <- U; Sigma2_u.chain[pos] <- Sigma2_u
      Mu.chain[, pos] <- Mu; preds.chain[, pos] <- preds; logit_bios.chain[, pos] <- logit_bios
    }
  }
  list(Beta.chain=Beta.chain, U.chain=U.chain, Sigma2_u.chain=Sigma2_u.chain,
       Mu.chain=Mu.chain, Preds=preds.chain, logit_bios.chain=logit_bios.chain)
}

unis_gaus <- function(X, Y, S, sig2b = 1000, wgt = NULL, n = NULL,
                      predX, predS, nburn = 1000, nsim = 5000, nthin = 1,
                      a = 0.1, b = 0.1, a_eps = 0.1, b_eps = 0.1) {
  N <- length(Y); p <- ncol(X); r <- ncol(S); npred <- nrow(predX)
  if (is.null(wgt)) wgt <- rep(1, N)

  X_w_X <- crossprod(X, wgt * X)
  S_w_S <- crossprod(S, wgt * S)
  X_w_S <- crossprod(X, wgt * S)
  S_w_X <- t(X_w_S)
  tX_wY <- crossprod(X, wgt * Y)
  tS_wY <- crossprod(S, wgt * Y)
  w_sum <- sum(wgt)
  Ip <- diag(p); Ir <- diag(r)

  Beta <- rep(1, p); Sigma2_u <- 1; U <- rep(0, r); Mu <- rep(0, N); sig2 <- 1
  n_keep <- nsim / nthin
  Beta.chain <- array(0, dim = c(p, n_keep)); U.chain <- array(0, dim = c(r, n_keep))
  Sigma2_u.chain <- numeric(n_keep); sig2.chain <- numeric(n_keep)
  Mu.chain <- array(0, dim = c(N, n_keep)); preds.chain <- array(0, dim = c(npred, n_keep))

  pb <- txtProgressBar(min = 0, max = nsim + nburn, style = 3)

  for (index in seq_len(nsim + nburn)) {
    sig2 <- 1 / stats::rgamma(1, shape = a_eps + w_sum / 2, rate = b_eps + 0.5 * sum(wgt * (Y - Mu)^2))

    Q_Beta  <- X_w_X / sig2 + Ip / sig2b
    b_Beta  <- (tX_wY - X_w_S %*% U) / sig2
    Beta    <- rmvn_prec(b_Beta, Q_Beta)

    b_prime  <- b + 0.5 * sum(U^2)
    Sigma2_u <- 1 / stats::rgamma(1, shape = a + r / 2, rate = b_prime)

    Q_U  <- S_w_S / sig2 + Ir / Sigma2_u
    b_U  <- (tS_wY - S_w_X %*% Beta) / sig2
    U    <- rmvn_prec(b_U, Q_U)

    Mu    <- as.vector(X %*% Beta + S %*% U)
    preds <- predX %*% Beta + predS %*% U

    utils::setTxtProgressBar(pb, index)
    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) / nthin
      Beta.chain[, pos] <- Beta; U.chain[, pos] <- U; Sigma2_u.chain[pos] <- Sigma2_u
      Mu.chain[, pos] <- Mu; sig2.chain[pos] <- sig2; preds.chain[, pos] <- preds
    }
  }
  list(Beta.chain=Beta.chain, U.chain=U.chain, Sigma2_u.chain=Sigma2_u.chain,
       Mu.chain=Mu.chain, sig2.chain=sig2.chain, Preds=preds.chain)
}

# MTSM_br <- function(X_1, X_2, Z_1, Z_2, S, sig2b = 1000, wgt = NULL, n_binom = NULL,
#                     predX, predS, n_preds, nburn = 1000, nsim = 5000, nthin = 1,
#                     sig2t = 10, sig2e = 10, tau_1_init = 1,
#                     a_eps = 0.1, b_eps = 0.1, aeta = 0.1, beta = 0.1,
#                     alambda = 0.1, blambda = 0.1) {

#   N <- nrow(X_1); r <- ncol(S); npred <- nrow(predX)
#   if (is.null(wgt)) wgt <- rep(1, N)
#   if (is.null(n_binom)) n_binom <- rep(1, N)

#   w_sum   <- sum(wgt)
#   X1_w_X1 <- crossprod(X_1, wgt * X_1)
#   S_w_S   <- crossprod(S, wgt * S)
#   X1_w_S  <- crossprod(X_1, wgt * S)
#   S_w_X1  <- t(X1_w_S)
#   tX1_wZ1 <- crossprod(X_1, wgt * Z_1)
#   tS_wZ1  <- crossprod(S, wgt * Z_1)
  
#   k       <- wgt * (Z_2 - n_binom / 2)
#   tX2_k   <- crossprod(X_2, k)
#   tS_k    <- crossprod(S, k)

#   p_1 <- ncol(X_1); p_2 <- ncol(X_2)
#   Ip1 <- diag(p_1); Ip2 <- diag(p_2); Ir <- diag(r)

#   tau_1 <- tau_1_init
#   Beta_1 <- rep(1, p_1); Beta_2 <- rep(1, p_2)
#   lambda <- rep(0, r); eta <- rep(1, r)
#   Mu_1 <- rep(1, N); Mu_2 <- rep(1, N); sig2 <- 1
  
#   n_keep <- nsim / nthin
#   tau_1.chain <- Sigma2_lambda.chain <- Sigma2_eta.chain <- sig2.chain <- numeric(n_keep)
#   Beta_1.chain <- array(0, dim = c(p_1, n_keep)); Beta_2.chain <- array(0, dim = c(p_2, n_keep))
#   lambda.chain <- array(0, dim = c(r, n_keep)); eta.chain <- array(0, dim = c(r, n_keep))
#   Mu_1.chain <- Mu_2.chain <- array(0, dim = c(N, n_keep))
#   preds_gaus.chain <- preds_bios.chain <- logit_bios.chain <- array(0, dim = c(npred, n_keep))

#   pb <- txtProgressBar(min = 0, max = nsim + nburn, style = 3)

#   for (index in seq_len(nsim + nburn)) {
#     sig2 <- 1 / stats::rgamma(1, shape = a_eps + w_sum / 2, rate = b_eps + 0.5 * sum(wgt * (Z_1 - Mu_1)^2))

#     omega <- BayesLogit::rpg.gamma(N, wgt, Mu_2)
#     OM <- crossprod(S, omega * S)
#     X2_w_X2 <- crossprod(X_2, omega * X_2)
#     S_w_X2 <- crossprod(S, omega * X_2)

#     sig2e <- 1 / stats::rgamma(1, shape = aeta + r / 2, rate = beta + 0.5 * sum(eta^2))
    
#     Q_eta <- tau_1^2 * (S_w_S / sig2) + Ir / sig2e + OM
#     b_eta <- tau_1 * (tS_wZ1 - S_w_X1 %*% Beta_1) / sig2 + (tS_k - S_w_X2 %*% Beta_2 - OM %*% lambda)
#     eta <- rmvn_prec(b_eta, Q_eta)
#     eta <- eta - mean(eta)

#     Q_Beta_1  <- (X1_w_X1 / sig2) + Ip1 / sig2b
#     b_Beta_1  <- (tX1_wZ1 - tau_1 * X1_w_S %*% eta) / sig2
#     Beta_1    <- rmvn_prec(b_Beta_1, Q_Beta_1)

#     sig2l <- 1 / stats::rgamma(1, shape = alambda + r / 2, rate = blambda + 0.5 * sum(lambda^2))
#     Q_lambda <- OM + Ir / sig2l
#     b_lambda <- tS_k - S_w_X2 %*% Beta_2 - OM %*% eta
#     lambda <- rmvn_prec(b_lambda, Q_lambda)
#     lambda <- lambda - mean(lambda)

#     Q_Beta_2  <- X2_w_X2 + Ip2 / sig2b
#     b_Beta_2  <- tX2_k - t(S_w_X2) %*% eta - t(S_w_X2) %*% lambda
#     Beta_2    <- rmvn_prec(b_Beta_2, Q_Beta_2)

#     M_sig2 <- S_w_S / sig2
#     var_tau_1 <- 1 / (sum(eta * (M_sig2 %*% eta)) + 1 / sig2t) 
#     mean_tau_1 <- var_tau_1 * sum(eta * ((tS_wZ1 - S_w_X1 %*% Beta_1) / sig2))
#     tau_1 <- stats::rnorm(1, mean_tau_1, sqrt(var_tau_1))

#     Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta)
#     Mu_2 <- as.vector(X_2 %*% Beta_2 + S %*% eta + S %*% lambda)
    
#     preds_gaus <- predX %*% Beta_1 + tau_1 * predS %*% eta
#     logit_bios <- predX %*% Beta_2 + predS %*% eta + predS %*% lambda
#     preds_bios <- stats::plogis(logit_bios)

#     utils::setTxtProgressBar(pb, index)
#     if (index > nburn && (index - nburn) %% nthin == 0) {
#       pos <- (index - nburn) / nthin
#       tau_1.chain[pos] <- tau_1; sig2.chain[pos] <- sig2
#       Sigma2_lambda.chain[pos] <- sig2l; Sigma2_eta.chain[pos] <- sig2e
#       Beta_1.chain[, pos] <- Beta_1; Beta_2.chain[, pos] <- Beta_2
#       lambda.chain[, pos] <- lambda; eta.chain[, pos] <- eta
#       Mu_1.chain[, pos] <- Mu_1; Mu_2.chain[, pos] <- Mu_2
#       preds_gaus.chain[, pos] <- preds_gaus; preds_bios.chain[, pos] <- preds_bios; logit_bios.chain[, pos] <- logit_bios
#     }
#   }

#   list(Beta_1.chain=Beta_1.chain, Beta_2.chain=Beta_2.chain, lambda.chain=lambda.chain, eta.chain=eta.chain,
#        Sigma2_lambda.chain=Sigma2_lambda.chain, Sigma2_eta.chain=Sigma2_eta.chain, sig2.chain=sig2.chain,
#        Mu_1.chain=Mu_1.chain, Mu_2.chain=Mu_2.chain, preds_gaus.chain=preds_gaus.chain, preds_bios.chain=preds_bios.chain,
#        logit_bios.chain=logit_bios.chain, tau_1.chain=tau_1.chain)
# }
## Multitype spatial model: Binomial with random effect only in binomial block
MTSM_br <- function(X_1, X_2, Z_1, Z_2, S, area_idx = NULL,
                    sig2b = 1000, wgt = NULL, n = NULL,
                    predX, predS, n_preds = NULL,
                    nburn = 1000, nsim = 5000, nthin = 1,
                    sig2t = 10, tau_1 = 1, tau_2_init = 1,
                    a_eps = 0.1, b_eps = 0.1,
                    aeta = 0.1, beta = 0.1,
                    alambda = 0.1, blambda = 0.1) {

  N <- nrow(X_1)
  if (nrow(X_2) != N) stop("X_1 and X_2 must have the same number of rows.")
  if (length(Z_1) != N || length(Z_2) != N) stop("Z_1 and Z_2 must have length N.")
  if (nrow(S) != N) stop("S must have N rows.")
  if (is.null(wgt)) wgt <- rep(1, N)
  if (is.null(n)) n <- rep(1, N)

  npred <- nrow(predX)
  if (nrow(predS) != npred) stop("predS must have the same number of rows as predX.")
  if (ncol(predS) != ncol(S)) stop("predS must have the same number of columns as S.")

  p_1 <- ncol(X_1)
  p_2 <- ncol(X_2)
  r <- ncol(S)

  Ip1 <- diag(p_1)
  Ip2 <- diag(p_2)
  Ir <- diag(r)

  w_sum <- sum(wgt)

  X1_w_X1 <- crossprod(X_1, wgt * X_1)
  S_w_S   <- crossprod(S, wgt * S)
  X1_w_S  <- crossprod(X_1, wgt * S)
  S_w_X1  <- t(X1_w_S)
  tX1_wZ1 <- crossprod(X_1, wgt * Z_1)
  tS_wZ1  <- crossprod(S, wgt * Z_1)

  kappa <- wgt * (Z_2 - n / 2)
  tX2_k <- crossprod(X_2, kappa)
  tS_k  <- crossprod(S, kappa)

  use_opt <- !is.null(area_idx)
  if (use_opt) {
    area_idx <- as.character(area_idx)
    area_fac <- factor(area_idx, levels = unique(area_idx))
    first_idx <- match(levels(area_fac), area_idx)
    S_unique <- S[first_idx, , drop = FALSE]
  }

  n_keep <- floor(nsim / nthin)

  tau_2 <- tau_2_init
  Beta_1 <- rep(0, p_1)
  Beta_2 <- rep(0, p_2)
  eta <- rep(0, r)
  lambda <- rep(0, r)
  Mu_1 <- rep(0, N)
  Mu_2 <- rep(0, N)
  sig2 <- 1
  sig2e <- 1
  sig2l <- 1

  tau_1.chain <- rep(tau_1, n_keep)
  tau_2.chain <- numeric(n_keep)
  Sigma2_lambda.chain <- numeric(n_keep)
  Sigma2_eta.chain <- numeric(n_keep)
  sig2.chain <- numeric(n_keep)

  Beta_1.chain <- array(0, dim = c(p_1, n_keep))
  Beta_2.chain <- array(0, dim = c(p_2, n_keep))
  eta.chain <- array(0, dim = c(r, n_keep))
  lambda.chain <- array(0, dim = c(r, n_keep))
  Mu_1.chain <- array(0, dim = c(N, n_keep))
  Mu_2.chain <- array(0, dim = c(N, n_keep))
  preds_gaus.chain <- array(0, dim = c(npred, n_keep))
  preds_bios.chain <- array(0, dim = c(npred, n_keep))
  logit_bios.chain <- array(0, dim = c(npred, n_keep))

  pb <- txtProgressBar(min = 0, max = nsim + nburn, style = 3)

  for (index in seq_len(nsim + nburn)) {
    sig2 <- 1 / stats::rgamma(
      1,
      shape = a_eps + w_sum / 2,
      rate = b_eps + 0.5 * sum(wgt * (Z_1 - Mu_1)^2)
    )

    omega <- BayesLogit::rpg.gamma(N, wgt * n, Mu_2)
    X2_w_X2 <- crossprod(X_2, omega * X_2)

    if (use_opt) {
      omega_area <- as.vector(rowsum(omega, area_fac, reorder = FALSE))
      OM <- crossprod(S_unique, omega_area * S_unique)

      omega_X2_area <- rowsum(omega * X_2, area_fac, reorder = FALSE)
      S_w_X2 <- crossprod(S_unique, omega_X2_area)
    } else {
      OM <- crossprod(S, omega * S)
      S_w_X2 <- crossprod(S, omega * X_2)
    }

    sig2e <- 1 / stats::rgamma(
      1,
      shape = aeta + r / 2,
      rate = beta + 0.5 * sum(eta^2)
    )

    Q_eta <- (tau_1^2) * (S_w_S / sig2) + (tau_2^2) * OM + Ir / sig2e
    b_eta <- tau_1 * (tS_wZ1 - S_w_X1 %*% Beta_1) / sig2 +
      tau_2 * (tS_k - S_w_X2 %*% Beta_2 - OM %*% lambda)
    eta <- rmvn_prec(b_eta, Q_eta)

    Q_Beta_1 <- X1_w_X1 / sig2 + Ip1 / sig2b
    b_Beta_1 <- (tX1_wZ1 - tau_1 * X1_w_S %*% eta) / sig2
    Beta_1 <- rmvn_prec(b_Beta_1, Q_Beta_1)

    sig2l <- 1 / stats::rgamma(
      1,
      shape = alambda + r / 2,
      rate = blambda + 0.5 * sum(lambda^2)
    )

    Q_lambda <- OM + Ir / sig2l
    b_lambda <- tS_k - S_w_X2 %*% Beta_2 - tau_2 * OM %*% eta
    lambda <- rmvn_prec(b_lambda, Q_lambda)

    Q_Beta_2 <- X2_w_X2 + Ip2 / sig2b
    b_Beta_2 <- tX2_k - t(S_w_X2) %*% (tau_2 * eta + lambda)
    Beta_2 <- rmvn_prec(b_Beta_2, Q_Beta_2)

    var_tau_2 <- 1 / (sum(eta * (OM %*% eta)) + 1 / sig2t)
    mean_tau_2 <- var_tau_2 * sum(eta * (tS_k - S_w_X2 %*% Beta_2 - OM %*% lambda))
    tau_2 <- stats::rnorm(1, mean_tau_2, sqrt(var_tau_2))

    Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta)
    Mu_2 <- as.vector(X_2 %*% Beta_2 + tau_2 * S %*% eta + S %*% lambda)

    preds_gaus <- as.vector(predX %*% Beta_1 + tau_1 * predS %*% eta)
    logit_bios <- as.vector(predX %*% Beta_2 + tau_2 * predS %*% eta + predS %*% lambda)
    preds_bios <- stats::plogis(logit_bios)

    utils::setTxtProgressBar(pb, index)

    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) / nthin
      tau_2.chain[pos] <- tau_2
      sig2.chain[pos] <- sig2
      Sigma2_eta.chain[pos] <- sig2e
      Sigma2_lambda.chain[pos] <- sig2l
      Beta_1.chain[, pos] <- Beta_1
      Beta_2.chain[, pos] <- Beta_2
      eta.chain[, pos] <- eta
      lambda.chain[, pos] <- lambda
      Mu_1.chain[, pos] <- Mu_1
      Mu_2.chain[, pos] <- Mu_2
      preds_gaus.chain[, pos] <- preds_gaus
      preds_bios.chain[, pos] <- preds_bios
      logit_bios.chain[, pos] <- logit_bios
    }
  }

  close(pb)

  list(
    Beta_1.chain = Beta_1.chain,
    Beta_2.chain = Beta_2.chain,
    eta.chain = eta.chain,
    lambda.chain = lambda.chain,
    Sigma2_eta.chain = Sigma2_eta.chain,
    Sigma2_lambda.chain = Sigma2_lambda.chain,
    sig2.chain = sig2.chain,
    Mu_1.chain = Mu_1.chain,
    Mu_2.chain = Mu_2.chain,
    preds_gaus.chain = preds_gaus.chain,
    preds_bios.chain = preds_bios.chain,
    logit_bios.chain = logit_bios.chain,
    tau_1.chain = tau_1.chain,
    tau_2.chain = tau_2.chain
  )
}

## functions for basis function version
build_adj_basis_abs <- function(area_sf, area_id, q = 0.25,
                                queen = TRUE, style = "B", zero.policy = TRUE) {
  area_sf <- area_sf |>
    dplyr::arrange(.data[[area_id]])

  nb <- spdep::poly2nb(area_sf, queen = queen)
  W <- spdep::nb2mat(nb, style = style, zero.policy = zero.policy)
  W <- (W + t(W)) / 2

  eig <- eigen(W, symmetric = TRUE)
  keep <- large_abs(eig$values, q)

  basis_mat <- eig$vectors[, keep, drop = FALSE]

  rownames(basis_mat) <- as.character(area_sf[[area_id]])
  colnames(basis_mat) <- paste0("bf", seq_len(ncol(basis_mat)))

  list(
    basis = basis_mat,
    W = W,
    eig_values = eig$values,
    kept = keep
  )
}

# MTSM_basis <- function(X_1, X_2, Z_1, Z_2, S, area_idx = NULL, sig2b = 1000,
#                        wgt = NULL, n_binom = NULL, predX, predS, n_preds = NULL,
#                        nburn = 1000, nsim = 5000, nthin = 1,
#                        sig2t = 10, tau_1_init = 0,
#                        a_eps = 0.1, b_eps = 0.1,
#                        aeta = 0.1, beta = 0.1,
#                        alambda = 0.1, blambda = 0.1) {

#   N <- nrow(X_1)
#   r <- ncol(S)
#   npred <- nrow(predX)

#   if (nrow(X_2) != N) stop("X_1 and X_2 must have the same number of rows.")
#   if (length(Z_1) != N || length(Z_2) != N) stop("Z_1 and Z_2 must have length N.")
#   if (nrow(S) != N) stop("S must have N rows.")
#   if (nrow(predS) != npred) stop("predS must have the same number of rows as predX.")
#   if (ncol(predS) != r) stop("predS must have the same number of columns as S.")

#   if (is.null(wgt)) wgt <- rep(1, N)
#   if (is.null(n_binom)) n_binom <- rep(1, N)

#   p_1 <- ncol(X_1)
#   p_2 <- ncol(X_2)

#   Ip1 <- diag(p_1)
#   Ip2 <- diag(p_2)
#   Ir <- diag(r)

#   w_sum <- sum(wgt)

#   X1_w_X1 <- crossprod(X_1, wgt * X_1)
#   S_w_S   <- crossprod(S, wgt * S)
#   X1_w_S  <- crossprod(X_1, wgt * S)
#   S_w_X1  <- t(X1_w_S)
#   tX1_wZ1 <- crossprod(X_1, wgt * Z_1)
#   tS_wZ1  <- crossprod(S, wgt * Z_1)

#   kappa <- wgt * (Z_2 - n_binom / 2)
#   tX2_k <- crossprod(X_2, kappa)
#   tS_k  <- crossprod(S, kappa)

#   use_opt <- !is.null(area_idx)
#   if (use_opt) {
#     area_idx <- as.character(area_idx)
#     area_fac <- factor(area_idx, levels = unique(area_idx))
#     first_idx <- match(levels(area_fac), area_idx)
#     S_unique <- S[first_idx, , drop = FALSE]
#   }

#   n_keep <- floor(nsim / nthin)

#   tau_1 <- tau_1_init
#   Beta_1 <- rep(0, p_1)
#   Beta_2 <- rep(0, p_2)
#   eta <- rep(0, r)
#   lambda <- rep(0, r)
#   Mu_1 <- rep(0, N)
#   Mu_2 <- rep(0, N)
#   sig2 <- 1
#   sig2e <- 1
#   sig2l <- 1

#   tau_1.chain <- numeric(n_keep)
#   Sigma2_lambda.chain <- numeric(n_keep)
#   Sigma2_eta.chain <- numeric(n_keep)
#   sig2.chain <- numeric(n_keep)

#   Beta_1.chain <- array(0, dim = c(p_1, n_keep))
#   Beta_2.chain <- array(0, dim = c(p_2, n_keep))
#   eta.chain <- array(0, dim = c(r, n_keep))
#   lambda.chain <- array(0, dim = c(r, n_keep))
#   Mu_1.chain <- array(0, dim = c(N, n_keep))
#   Mu_2.chain <- array(0, dim = c(N, n_keep))
#   preds_gaus.chain <- array(0, dim = c(npred, n_keep))
#   preds_bios.chain <- array(0, dim = c(npred, n_keep))
#   logit_bios.chain <- array(0, dim = c(npred, n_keep))

#   pb <- txtProgressBar(min = 0, max = nsim + nburn, style = 3)

#   for (index in seq_len(nsim + nburn)) {
#     sig2 <- 1 / stats::rgamma(
#       1,
#       shape = a_eps + w_sum / 2,
#       rate = b_eps + 0.5 * sum(wgt * (Z_1 - Mu_1)^2)
#     )

#     omega <- BayesLogit::rpg.gamma(N, wgt * n_binom, Mu_2)

#     X2_w_X2 <- crossprod(X_2, omega * X_2)

#     if (use_opt) {
#       omega_area <- as.vector(rowsum(omega, area_fac, reorder = FALSE))
#       OM <- crossprod(S_unique, omega_area * S_unique)

#       omega_X2_area <- rowsum(omega * X_2, area_fac, reorder = FALSE)
#       S_w_X2 <- crossprod(S_unique, omega_X2_area)
#     } else {
#       OM <- crossprod(S, omega * S)
#       S_w_X2 <- crossprod(S, omega * X_2)
#     }

#     sig2e <- 1 / stats::rgamma(
#       1,
#       shape = aeta + r / 2,
#       rate = beta + 0.5 * sum(eta^2)
#     )

#     Q_eta <- tau_1^2 * (S_w_S / sig2) + OM + Ir / sig2e
#     b_eta <- tau_1 * (tS_wZ1 - S_w_X1 %*% Beta_1) / sig2 +
#       (tS_k - S_w_X2 %*% Beta_2 - OM %*% lambda)
#     eta <- rmvn_prec(b_eta, Q_eta)

#     Q_Beta_1 <- X1_w_X1 / sig2 + Ip1 / sig2b
#     b_Beta_1 <- (tX1_wZ1 - tau_1 * X1_w_S %*% eta) / sig2
#     Beta_1 <- rmvn_prec(b_Beta_1, Q_Beta_1)

#     sig2l <- 1 / stats::rgamma(
#       1,
#       shape = alambda + r / 2,
#       rate = blambda + 0.5 * sum(lambda^2)
#     )

#     Q_lambda <- OM + Ir / sig2l
#     b_lambda <- tS_k - S_w_X2 %*% Beta_2 - OM %*% eta
#     lambda <- rmvn_prec(b_lambda, Q_lambda)

#     Q_Beta_2 <- X2_w_X2 + Ip2 / sig2b
#     b_Beta_2 <- tX2_k - t(S_w_X2) %*% eta - t(S_w_X2) %*% lambda
#     Beta_2 <- rmvn_prec(b_Beta_2, Q_Beta_2)

#     M_sig2 <- S_w_S / sig2
#     var_tau_1 <- 1 / (sum(eta * (M_sig2 %*% eta)) + 1 / sig2t)
#     mean_tau_1 <- var_tau_1 * sum(eta * ((tS_wZ1 - S_w_X1 %*% Beta_1) / sig2))
#     tau_1 <- stats::rnorm(1, mean_tau_1, sqrt(var_tau_1))

#     Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta)
#     Mu_2 <- as.vector(X_2 %*% Beta_2 + S %*% eta + S %*% lambda)

#     preds_gaus <- as.vector(predX %*% Beta_1 + tau_1 * predS %*% eta)
#     logit_bios <- as.vector(predX %*% Beta_2 + predS %*% eta + predS %*% lambda)
#     preds_bios <- stats::plogis(logit_bios)

#     utils::setTxtProgressBar(pb, index)

#     if (index > nburn && (index - nburn) %% nthin == 0) {
#       pos <- (index - nburn) / nthin
#       tau_1.chain[pos] <- tau_1
#       sig2.chain[pos] <- sig2
#       Sigma2_eta.chain[pos] <- sig2e
#       Sigma2_lambda.chain[pos] <- sig2l
#       Beta_1.chain[, pos] <- Beta_1
#       Beta_2.chain[, pos] <- Beta_2
#       eta.chain[, pos] <- eta
#       lambda.chain[, pos] <- lambda
#       Mu_1.chain[, pos] <- Mu_1
#       Mu_2.chain[, pos] <- Mu_2
#       preds_gaus.chain[, pos] <- preds_gaus
#       preds_bios.chain[, pos] <- preds_bios
#       logit_bios.chain[, pos] <- logit_bios
#     }
#   }

#   list(
#     Beta_1.chain = Beta_1.chain,
#     Beta_2.chain = Beta_2.chain,
#     eta.chain = eta.chain,
#     lambda.chain = lambda.chain,
#     Sigma2_eta.chain = Sigma2_eta.chain,
#     Sigma2_lambda.chain = Sigma2_lambda.chain,
#     sig2.chain = sig2.chain,
#     Mu_1.chain = Mu_1.chain,
#     Mu_2.chain = Mu_2.chain,
#     preds_gaus.chain = preds_gaus.chain,
#     preds_bios.chain = preds_bios.chain,
#     logit_bios.chain = logit_bios.chain,
#     tau_1.chain = tau_1.chain
#   )
# }

MTSM_basis <- function(X_1, X_2, Z_1, Z_2, S, area_idx = NULL,
                       sig2b = 1000, wgt = NULL, n_binom = NULL,
                       predX, predS, n_preds = NULL,
                       nburn = 1000, nsim = 5000, nthin = 1,
                       sig2t = 10, tau_1 = 1, tau_2_init = 1,
                       a_eps = 0.1, b_eps = 0.1,
                       aeta = 0.1, beta = 0.1,
                       alambda = 0.1, blambda = 0.1) {

  N <- nrow(X_1)
  if (nrow(X_2) != N) stop("X_1 and X_2 must have the same number of rows.")
  if (length(Z_1) != N || length(Z_2) != N) stop("Z_1 and Z_2 must have length N.")
  if (nrow(S) != N) stop("S must have N rows.")
  if (is.null(wgt)) wgt <- rep(1, N)
  if (is.null(n_binom)) n_binom <- rep(1, N)

  npred <- nrow(predX)
  if (nrow(predS) != npred) stop("predS must have the same number of rows as predX.")
  if (ncol(predS) != ncol(S)) stop("predS must have the same number of columns as S.")

  p_1 <- ncol(X_1)
  p_2 <- ncol(X_2)
  r <- ncol(S)

  Ip1 <- diag(p_1)
  Ip2 <- diag(p_2)
  Ir <- diag(r)

  w_sum <- sum(wgt)

  X1_w_X1 <- crossprod(X_1, wgt * X_1)
  S_w_S   <- crossprod(S, wgt * S)
  X1_w_S  <- crossprod(X_1, wgt * S)
  S_w_X1  <- t(X1_w_S)
  tX1_wZ1 <- crossprod(X_1, wgt * Z_1)
  tS_wZ1  <- crossprod(S, wgt * Z_1)

  kappa <- wgt * (Z_2 - n_binom / 2)
  tX2_k <- crossprod(X_2, kappa)
  tS_k  <- crossprod(S, kappa)

  use_opt <- !is.null(area_idx)
  if (use_opt) {
    area_idx <- as.character(area_idx)
    area_fac <- factor(area_idx, levels = unique(area_idx))
    first_idx <- match(levels(area_fac), area_idx)
    S_unique <- S[first_idx, , drop = FALSE]
  }

  n_keep <- floor(nsim / nthin)

  tau_2 <- tau_2_init
  Beta_1 <- rep(0, p_1)
  Beta_2 <- rep(0, p_2)
  eta <- rep(0, r)
  lambda <- rep(0, r)
  Mu_1 <- rep(0, N)
  Mu_2 <- rep(0, N)
  sig2 <- 1
  sig2e <- 1
  sig2l <- 1

  tau_1.chain <- rep(tau_1, n_keep)
  tau_2.chain <- numeric(n_keep)
  Sigma2_lambda.chain <- numeric(n_keep)
  Sigma2_eta.chain <- numeric(n_keep)
  sig2.chain <- numeric(n_keep)

  Beta_1.chain <- array(0, dim = c(p_1, n_keep))
  Beta_2.chain <- array(0, dim = c(p_2, n_keep))
  eta.chain <- array(0, dim = c(r, n_keep))
  lambda.chain <- array(0, dim = c(r, n_keep))
  Mu_1.chain <- array(0, dim = c(N, n_keep))
  Mu_2.chain <- array(0, dim = c(N, n_keep))
  preds_gaus.chain <- array(0, dim = c(npred, n_keep))
  preds_bios.chain <- array(0, dim = c(npred, n_keep))
  logit_bios.chain <- array(0, dim = c(npred, n_keep))

  pb <- txtProgressBar(min = 0, max = nsim + nburn, style = 3)

  for (index in seq_len(nsim + nburn)) {
    sig2 <- 1 / stats::rgamma(
      1,
      shape = a_eps + w_sum / 2,
      rate = b_eps + 0.5 * sum(wgt * (Z_1 - Mu_1)^2)
    )

    omega <- BayesLogit::rpg.gamma(N, wgt * n_binom, Mu_2)
    X2_w_X2 <- crossprod(X_2, omega * X_2)

    if (use_opt) {
      omega_area <- as.vector(rowsum(omega, area_fac, reorder = FALSE))
      OM <- crossprod(S_unique, omega_area * S_unique)

      omega_X2_area <- rowsum(omega * X_2, area_fac, reorder = FALSE)
      S_w_X2 <- crossprod(S_unique, omega_X2_area)
    } else {
      OM <- crossprod(S, omega * S)
      S_w_X2 <- crossprod(S, omega * X_2)
    }

    sig2e <- 1 / stats::rgamma(
      1,
      shape = aeta + r / 2,
      rate = beta + 0.5 * sum(eta^2)
    )

    Q_eta <- (tau_1^2) * (S_w_S / sig2) + (tau_2^2) * OM + Ir / sig2e
    b_eta <- tau_1 * (tS_wZ1 - S_w_X1 %*% Beta_1) / sig2 +
      tau_2 * (tS_k - S_w_X2 %*% Beta_2 - OM %*% lambda)
    eta <- rmvn_prec(b_eta, Q_eta)

    Q_Beta_1 <- X1_w_X1 / sig2 + Ip1 / sig2b
    b_Beta_1 <- (tX1_wZ1 - tau_1 * X1_w_S %*% eta) / sig2
    Beta_1 <- rmvn_prec(b_Beta_1, Q_Beta_1)

    sig2l <- 1 / stats::rgamma(
      1,
      shape = alambda + r / 2,
      rate = blambda + 0.5 * sum(lambda^2)
    )

    Q_lambda <- OM + Ir / sig2l
    b_lambda <- tS_k - S_w_X2 %*% Beta_2 - tau_2 * OM %*% eta
    lambda <- rmvn_prec(b_lambda, Q_lambda)

    Q_Beta_2 <- X2_w_X2 + Ip2 / sig2b
    b_Beta_2 <- tX2_k - t(S_w_X2) %*% (tau_2 * eta + lambda)
    Beta_2 <- rmvn_prec(b_Beta_2, Q_Beta_2)

    var_tau_2 <- 1 / (sum(eta * (OM %*% eta)) + 1 / sig2t)
    mean_tau_2 <- var_tau_2 * sum(eta * (tS_k - S_w_X2 %*% Beta_2 - OM %*% lambda))
    tau_2 <- stats::rnorm(1, mean_tau_2, sqrt(var_tau_2))

    Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta)
    Mu_2 <- as.vector(X_2 %*% Beta_2 + tau_2 * S %*% eta + S %*% lambda)

    preds_gaus <- as.vector(predX %*% Beta_1 + tau_1 * predS %*% eta)
    logit_bios <- as.vector(predX %*% Beta_2 + tau_2 * predS %*% eta + predS %*% lambda)
    preds_bios <- stats::plogis(logit_bios)

    utils::setTxtProgressBar(pb, index)

    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) / nthin
      tau_2.chain[pos] <- tau_2
      sig2.chain[pos] <- sig2
      Sigma2_eta.chain[pos] <- sig2e
      Sigma2_lambda.chain[pos] <- sig2l
      Beta_1.chain[, pos] <- Beta_1
      Beta_2.chain[, pos] <- Beta_2
      eta.chain[, pos] <- eta
      lambda.chain[, pos] <- lambda
      Mu_1.chain[, pos] <- Mu_1
      Mu_2.chain[, pos] <- Mu_2
      preds_gaus.chain[, pos] <- preds_gaus
      preds_bios.chain[, pos] <- preds_bios
      logit_bios.chain[, pos] <- logit_bios
    }
  }

  close(pb)

  list(
    Beta_1.chain = Beta_1.chain,
    Beta_2.chain = Beta_2.chain,
    eta.chain = eta.chain,
    lambda.chain = lambda.chain,
    Sigma2_eta.chain = Sigma2_eta.chain,
    Sigma2_lambda.chain = Sigma2_lambda.chain,
    sig2.chain = sig2.chain,
    Mu_1.chain = Mu_1.chain,
    Mu_2.chain = Mu_2.chain,
    preds_gaus.chain = preds_gaus.chain,
    preds_bios.chain = preds_bios.chain,
    logit_bios.chain = logit_bios.chain,
    tau_1.chain = tau_1.chain,
    tau_2.chain = tau_2.chain
  )
}