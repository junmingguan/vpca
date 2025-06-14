# Based on matlab code by Shira Kritchman and Boaz Nadler
# https://www.wisdom.weizmann.ac.il/~nadler/Rank_Estimation/rank_estimation.html

#' @export
KN_s_Wishart <- function(alpha, beta=1) {
  if (beta == 1) {
    s_Wishart <- (-3/2 * log(4 * sqrt(pi) * alpha))^(2/3)
    return(s_Wishart)
  } else {
    stop("beta must be 1.")
  }
}


#' @export
KN_mu_sigma <- function(n, p, beta=1) {
  if (beta == 1) {
    mu_np <- (sqrt(n - 1/2) + sqrt(p - 1/2))^2
    sigma_np <- sqrt(mu_np) * (1/sqrt(n - 1/2) + 1 / sqrt(p - 1/2))^(1/3)
    # mu_np <- (sqrt(n - 1/2) + sqrt(p - 1/2))^2 / n
    # sigma_np <- sqrt(mu_np / n) * (1/sqrt(n - 1/2) + 1 / sqrt(p - 1/2))^(1/3)
  } else {
    stop("beta must be 1.")
  }
  return(list(mu_np = mu_np, sigma_np = sigma_np))
}


#' @export
KN_noiseEst <- function(ell, n, kk) {
  #   ell     -  vector of eigenvalues of the SCM, of length p
  #   n       -  number of samples
  #   kk      -  assumed rank
  max_iter <- 30
  eps_threshold <- 1e-5

  p <- length(ell)
  sigma_0 <- (1/(p - kk)) * sum(ell[(kk + 1):p]) * 1 / (1 - kk / n)

  for (counter in 1:max_iter) {
    current_ells <- ell[1:kk]

    tmp <- current_ells + sigma_0 - (p - kk) / n * sigma_0
    discriminant <- tmp^2 - 4 * current_ells * sigma_0

    if (any(discriminant < 0)) {
      # otherwise get complex valued solutions
      break
    }

    Rho_est <- numeric(kk) # Initialize Rho_est as a numeric vector
    Rho_est <- (tmp + sqrt(discriminant)) / 2

    if (any(current_ells - Rho_est < 0)) {
      break
    }
    Delta_l_rho <- pmax(current_ells - Rho_est, 0)
    sigma_new <- (1/(p - kk)) * (sum(ell[(kk + 1):p]) + sum(Delta_l_rho))

    if (abs(sigma_new - sigma_0) / sigma_0 < eps_threshold) {
      break
    } else {
      sigma_0 <- sigma_new
    }
  }

  sig_hat_kk <- sigma_0
  return(sig_hat_kk)
}

#' @export
KN_rankEst <- function(ell, n, beta, alpha = 0.5, max_kk = NULL, sigma2 = NULL) {
  p <- length(ell)

  if (is.null(max_kk)) {
    max_kk <- min(n, p) - 1
  }
  max_kk <- min(max_kk, min(p, n) - 1)

  s_Wishart <- KN_s_Wishart(alpha, beta)

  sigma_arr <- numeric(max_kk)
  K_val <- 0 #

  if (max_kk < 1) {
    K_val <- 0
  } else {
    for (kk_idx in 1:max_kk) {
      mu_sigma_vals <- KN_mu_sigma(n, p - kk_idx, beta)
      mu_np <- mu_sigma_vals$mu_np
      sigma_np <- mu_sigma_vals$sigma_np
      if (is.null(sigma2)) {
        sig_hat_kk <- KN_noiseEst(ell, n, kk_idx)
      } else {
        sig_hat_kk = sigma2
      }
      sigma_arr[kk_idx] <- sig_hat_kk

      at_least_kk_signals <- n * ell[kk_idx] > sig_hat_kk * (mu_np + s_Wishart * sigma_np)
      if (!at_least_kk_signals) {
        K_val <- kk_idx - 1
        break
      }
      if (kk_idx == max_kk) { # If loop completes without break
        K_val <- kk_idx
      }
    }
  }

  if (K_val > 0) {
    sigma_hat <- sigma_arr[K_val]
  } else {
    sigma_hat <- sum(ell[1:p]) / p
  }

  return(list(K = K_val, sigma_hat = sigma_hat))
}
