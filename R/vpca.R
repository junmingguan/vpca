# library(microbenchmark)


#' Update either \eqn{Q(\mathbf{X})}.
#'
#' @param obj vpca.data object.
#'
#' @return vpca.data object.
#' @import Matrix
update_X <- function(obj) {
  return(
    within(
      obj,
      {
        # if (missing(Ealpha)) {
        #   Sigma <- chol2inv(chol(E_tau * E_WtW + diag(q)))
        # }
        # else {
        #   Sigma <- chol2inv(chol(E_tau * E_WtW + diag(E_alpha)))
        # }
        if (alpha_for_X) {
          Sigma <- chol2inv(chol(diag(E_alpha_X) + E_tau * E_WtW))
        } else {
          Sigma <- chol2inv(chol(E_tau * E_WtW + diag(q)))
        }
        # print( 
        #   microbenchmark( E_tau * Sigma %*% t(E_W) %*% (Y - outer(E_mu, rep(1, n))) , E_tau * Sigma %*% ( t(E_W) %*% Y - drop( t(E_W) %*% E_mu ) ) , times=100 ) 
        # ) 
        # E_X <- E_tau * Sigma %*% t(E_W) %*% (Y - outer(E_mu, rep(1, n)))
        if (include_mu) {
          E_X <- E_tau * Sigma %*% ( t(E_W) %*% Y - drop( t(E_W) %*% E_mu ) )
        } else {
          E_X <- E_tau * Sigma %*% ( t(E_W) %*% Y )
        }
        E_XXt <- n * Sigma + E_X %*% t(E_X)
        H_X <- n * q / 2 * log(2 * pi * exp(1)) + n / 2 * determinant(Sigma, logarithm = TRUE)$modulus[1]
        if (alpha_for_X) {
          P_X <- n / 2 * sum(E_log_alpha_X) - n * d / 2 * log(2 * pi) - 0.5 * sum(diag(E_XXt) * E_alpha_X)
        }
        Sigma_X <- Sigma
        rm(Sigma)
      }
    )
  )
}

#' Update \eqn{\mathbf{W}}.
#'
#' @param obj vpca.data object.
#'
#' @return vpca.data object.
#' @import Matrix
update_W <- function(obj) {
  return(
    within(
      obj,
      {
        # if (missing(E_alpha)) {
        #   Sigma <- chol2inv(chol(diag(q) + E_tau * E_XXt))
        # }
        # else {
        #   Sigma <- chol2inv(chol(diag(E_alpha) + E_tau * E_XXt))
        # }
        Sigma <- chol2inv(chol(diag(E_alpha) + E_tau * E_XXt))

        # print( 
        #   microbenchmark( E_tau * Sigma %*% E_X %*% (t(Y) - outer(rep(1, n), E_mu)) , 
        #    E_tau * ( Y %*% t(E_X) %*% Sigma - outer(E_mu, rep(1, n)) %*% t(E_X) %*% Sigma ), 
        #    E_tau * ( Y %*% t(E_X) %*% Sigma - outer(E_mu, drop( Sigma %*% rowSums(E_X) ) ) ),
        #    E_tau * ( Y %*% t(E_X) %*% Sigma ),
        #    times=100 ) 
        # )

        # E_W <- E_tau * Sigma %*% E_X %*% (t(Y) - outer(rep(1, n), E_mu))
        # E_W <- t(E_W)
        if (include_mu) {
          # print(sapply(E_mu, function(x) x * Sigma %*% rowSums(E_X)))
          # E_W <- E_tau * ( Y %*% t(E_X) %*% Sigma - outer(E_mu, rep(1, n)) %*% t(E_X) %*% Sigma )
          E_W <- E_tau * ( Y %*% t(E_X) %*% Sigma - outer(E_mu, drop( Sigma %*% rowSums(E_X) ) ) )
        } else {
          E_W <- E_tau * ( Y %*% t(E_X) %*% Sigma )
        }
        E_WtW <- d * Sigma + t(E_W) %*% E_W
        H_W <- d * q / 2 * log(2 * pi * exp(1)) + d / 2 * determinant(Sigma, logarithm = TRUE)$modulus[1]
        P_W <- d / 2 * sum(E_log_alpha) - d * q / 2 * log(2 * pi) - 0.5 * sum(diag(E_WtW) * E_alpha)
        Sigma_W <- Sigma
        rm(Sigma)
      }
    )
  )
}

#' Update \eqn{Q(\mathbf{\mu})}.
#'
#' @param obj vpca.data object.
#'
#' @return vpca.data object.
#' @import Matrix
update_mu <- function(obj) {
  return(
    within(
      obj,
      {
        E_mu <- ( rowSums(Y) - E_W %*% rowSums(E_X) ) * ( E_tau / (beta + n * E_tau) )
        E_mu_Sum_Sq <- d / (beta + n * E_tau) + sum(E_mu^2)
        logdet_Sigma <- sum(log(d / (beta + n * E_tau)))
        H_mu <- d / 2 * log(2 * pi * exp(1)) + 0.5 * logdet_Sigma
        P_mu <- -d / 2 * log(2 * pi) + d / 2 * log(beta) - beta / 2 * (logdet_Sigma + sum(E_mu^2))
        rm(logdet_Sigma)
      }
    )
  )
}

#' Update \eqn{Q(\mathbf{\alpha_\mathbf{X}})}.
#'
#' @param obj vpca.data object.
#'
#' @return vpca.data object.
#' @import Matrix
update_alpha <- function(obj) {
  return(
    within(
      obj,
      {
        if (alpha_MLE) {
          E_alpha <- d / diag(E_WtW)
          E_log_alpha <- log(E_alpha)
        } else {
          a <- a_alpha + d / 2
          b <- b_alpha + diag(E_WtW) / 2
          E_alpha <- a / b
          E_log_alpha <- digamma(a) - log(b)
          H_alpha <- q * (a + lgamma(a) + (1 - a) * digamma(a)) - sum(log(b))
          P_alpha <- q * a_alpha * log(b_alpha) + (a_alpha - 1) * sum(E_log_alpha) - b_alpha * sum(E_alpha) - q * lgamma(a_alpha)
          rm(a)
          rm(b)
        }

        if (alpha_for_X) {
          if (alpha_MLE) {
            E_alpha_X <- n / diag(E_XXt)
            E_log_alpha_X <- log(E_alpha_X)
          } else {
            a <- a_alpha + n / 2
            b <- b_alpha + diag(E_XXt) / 2
            E_alpha_X <- a / b
            E_log_alpha_X <- digamma(a) - log(b)
            H_alpha_X <- q * (a + lgamma(a) + (1 - a) * digamma(a)) - sum(log(b))
            P_alpha_X <- q * a_alpha * log(b_alpha) + (a_alpha - 1) * sum(E_log_alpha_X) - b_alpha * sum(E_alpha_X) - q * lgamma(a_alpha)
            rm(a)
            rm(b)
          }
        }
      }
    )
  )
}

#' Update \eqn{Q(\tau)}.
#'
#' @param obj vpca.data object.
#'
#' @return vpca.data object.
#' @import Matrix
update_tau <- function(obj) {
  return(
    within(
      obj,
      {
        # print(
        #   microbenchmark( # drop( Reduce(`+`, lapply(1:n, function(j) t(Y[, j]) %*% E_W %*% E_X[, j]))),
        #   sum(t(Y) %*% E_W * t(E_X)),
        #   sum(diag(t(Y) %*% E_W %*% E_X)),
        #    times=100 )
        # )
        # print(sum(t(Y) %*% E_W * t(E_X))-sum(diag(t(Y) %*% E_W %*% E_X)))


        a <- a_tau + n * d / 2
        if (include_mu) {
          E_resid_Sq <- Y_Frob_Sq + n * E_mu_Sum_Sq + sum(diag(E_WtW %*% E_XXt)) +
            2 * t(E_mu) %*% E_W %*% rowSums(E_X) -
            2 * sum(t(Y) %*% E_W * t(E_X)) -
            2 * t(E_mu) %*% rowSums(Y)
        } else {
          E_resid_Sq <- Y_Frob_Sq + sum(diag(E_WtW %*% E_XXt)) -
            2 * sum(t(Y) %*% E_W * t(E_X))
        }
        b <- b_tau + E_resid_Sq / 2
        E_tau <- a / b[1]
        E_log_tau <- digamma(a) - log(b)
        H_tau <- a - log(b) + lgamma(a) + (1 - a) * digamma(a)
        P_tau <- a * log(b) - lgamma(a) + (a - 1) * E_log_tau - b_tau * E_tau
        rm(a)
        rm(b)
      }
    )
  )
}


#' Compute ELBO.
#'
#' @param obj vpca.data object.
#'
#' @return vpca.data object.
#' @import Matrix
update_elbo <- function(obj) {
  # complete likelihood
  return(
    within(
      obj,
      {
        complete_loglik <- d * n / 2 * (E_log_tau - log(2*pi)) +
          - E_resid_Sq * E_tau / 2

        entropy <- H_X + H_W + H_tau

        if (include_mu) {
          entropy <- entropy + H_mu
        }

        prior <- P_tau + P_W
        if (include_mu) {
          prior <- prior + P_mu
        }
        if (!alpha_MLE) {
          prior <- prior + P_W + P_alpha
          entropy <- entropy + H_alpha
          if (alpha_for_X) {
            prior <- prior + P_X + P_alpha_X
            entropy <- entropy + H_alpha_X
          }
        }

        elbo <- drop( complete_loglik + prior + entropy )
        rm(complete_loglik)
        rm(prior)
        rm(entropy)
      }
    )
  )
}

#' Variational PCA algorithm runner.
#'
#' @param init vpca.data object.
#' @param min.iter Minimum number of iterations.
#' @param max.iter Maximum number of iterations.
#' @param epsilon The algorithm terminates if the difference in ELBO is less than \eqn{\epsilon}.
#'
#' @return First moments and ELBOs in all iterations.
#' @import Matrix
#' @export
vpcar <- function(init, min.iter = 1000, max.iter = 3000, epsilon = 0.0001) {
  obj <- init
  count <- 0
  elbo_list <- c()
  elbo_old <- -Inf
  while (TRUE) {

    obj <- update_X(obj)
    if (obj$include_mu) {
      obj <- update_mu(obj)
    }
    obj <- update_W(obj)
    obj <- update_alpha(obj)
    obj <- update_tau(obj)
    obj <- update_elbo(obj)

    count <- count + 1
    if (obj$elbo - elbo_old < 0) {
      cat("\nELBO decreases at iteration no.", count, "\n")
    }
    diff <- abs(obj$elbo - elbo_old)
    elbo_old <- obj$elbo
    elbo_list <- c(elbo_list, elbo_old)
    if (diff <= epsilon && count >= min.iter) {
      break
    }
    if (count == max.iter) {
      cat("\nMaximum number of iteration reached.\n")
      break
    }
  }
  return(
    list(
      E_W = obj$E_W, E_X = obj$E_X, E_mu = obj$E_mu, Sigma_W = obj$Sigma_W, Sigma_X = obj$Sigma_X,
      E_alpha = obj$E_alpha, E_tau = obj$E_tau, elbo_list = elbo_list
    )
  )
}
