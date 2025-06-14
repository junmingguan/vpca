# Find exact tau_unbar using root finding
find_tau_unbar <- function(alpha) {
  phi <- function(z) (log(z + 1)/z) - 0.5
  return(uniroot(
    function(tau) phi(tau) + phi(tau/alpha),
    interval = c(sqrt(alpha), 2.5129),
    extendInt = "yes"
  )$root)
}

#' Iterative EVB analytic solution by Nakajima et al. 2015
#'
#' @export
evb_global_pca <- function(svd_result, L, M, H, sigma2 = NULL) {
  gamma <- svd_result$d
  gamma_sq <- gamma^2
  alpha <- L / M
  tau_unbar = 2.5129 * sqrt(alpha)
  H_upper <- min(floor(L / (1 + alpha)) - 1, H)

  sum_gamma_sq <- sum(gamma_sq)
  sigma_upper <- sum_gamma_sq / (L * M)
  threshold_x <- (1 + tau_unbar) * (1 + alpha/tau_unbar)
  if (H_upper >= length(gamma)) {
    sigma_lower <- 0
  } else {
    sum_gamma_lower <- sum(gamma_sq[(H_upper + 1):length(gamma_sq)])
    sigma_lower_term1 <- gamma_sq[H_upper + 1] / (M * threshold_x)
    sigma_lower_term2 <- sum_gamma_lower / (M * (L - H_upper))
    sigma_lower <- max(sigma_lower_term1, sigma_lower_term2)
  }

  compute_psi1 <- function(x) {
    res <- rep(0, length(x))
    for (i in 1:length(x)) {
      if (x[i] > threshold_x) {
        res[i] = 0.5 * ( x[i] - (1 + alpha) + sqrt( (x[i] - (1 + alpha))^2 - 4 * alpha) )
        res[i] = log(res[i] + 1) + alpha * log(res[i]/alpha + 1) - res[i]
      }
    }
    return(res)
  }

  # Theorem 7 in Nakajima et al. 2015
  objective <- function(sigma_inv_sq) {
    sigma_sq <- 1 / sigma_inv_sq
    x <- gamma_sq / (M * sigma_sq)

    # Compute psi (Eq.41)
    psi <- (x - log(x)) + compute_psi1(x)
    omega <- ( sum(psi[1:H_upper]) + sum( (x - log(x))[(H_upper + 1):length(x)] ) ) / L

    return(omega)
  }

  if (is.null(sigma2)) {
    lower_bound <-1 / sigma_upper
    upper_bound <- 1 / sigma_lower
    optim_result <- optimize(
      f = objective,
      interval = c(lower_bound, upper_bound),
      maximum = FALSE
    )

    sigma2 <- 1/optim_result$minimum
  }

  # (25) in Nakajima et al. 2015
  gamma_evb_threshold <- sqrt(M * sigma2 * threshold_x)

  gamma_evb <- numeric(length(gamma))
  for(h in 1:length(gamma)) {
    if(gamma[h] >= gamma_evb_threshold) {
      term <- ( 1 - (L + M)*sigma2/gamma[h]^2 )
      gamma_evb[h] <- gamma[h]/2 * (term + sqrt( term^2 - 4*L*M*(sigma2/gamma[h]^2)^2 ))
      # gamma_evb[h] <- gamma[h]/2 * (term + sqrt(pmax(term^2 - 4*L*M*(sigma2/gamma[h]^2)^2, 0)))
    } else {
      gamma_evb[h] <- 0
    }
  }
  estimated_rank <- sum(gamma_evb > 0)


  return(list(
    estimated_rank = estimated_rank,
    sigma2 = sigma2,
    gamma_evb = gamma_evb
  ))
}



#' EVB analytic global solution by Nakajima et al. 2015
#'
#' @export
evb_ite_pca <- function(svd_result, L, M, max_iter = 100, tol = 1e-6, local = TRUE) {
  alpha <- L/M

  gamma <- svd_result$d
  total_variance <- sum(gamma^2)
  tau_unbar = 2.5129 * sqrt(alpha)

  sigma2 <- 1e-4 * total_variance/(L*M)

  for(iter in 1:max_iter) {
    if (local) {
      # local EVB threshold; Eq. (28) in Nakajima et al. 2015
      gamma_evb_threshold <- sqrt(sigma2) * (sqrt(M) + sqrt(L))
    } else {
      gamma_evb_threshold <- sqrt(sigma2) * sqrt(M*(1 + tau_unbar)*(1 + alpha/tau_unbar))
    }

    gamma_evb <- numeric(length(gamma))
    for(h in 1:length(gamma)) {
      if(gamma[h] >= gamma_evb_threshold) {
        term <- (1 - (L + M)*sigma2/gamma[h]^2)
        gamma_evb[h] <- gamma[h]/2 * (term + sqrt( term^2 - 4*L*M*(sigma2/gamma[h]^2)^2 ))
        # gamma_evb[h] <- gamma[h]/2 * (term + sqrt(pmax(term^2 - 4*L*M*(sigma2/gamma[h]^2)^2, 0)))
      } else {
        gamma_evb[h] <- 0
      }
    }

    sigma2_new <- (total_variance - sum(gamma * gamma_evb)) / (L*M)

    delta <- abs(sigma2_new - sigma2)/sigma2
    if(delta < tol) break
    sigma2 <- sigma2_new
  }

  return(list(
    estimated_rank = sum(gamma_evb > 0),
    sigma2 = sigma2,
    gamma_evb = gamma_evb,
    iterations = iter
  ))
}
