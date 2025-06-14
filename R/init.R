#' @export
#' @import Matrix
init_random <- function(Y, q, alpha_MLE = FALSE, include_mu = TRUE, alpha_for_X = FALSE,
                        a_alpha = 0.001, b_alpha = 0.001,
                        a_tau = 0.001, b_tau = 0.001, beta = 0.001) {
  obj <- list(Y = Y, Y_Frob_Sq = sum(Y^2))
  n <- ncol(Y)
  d <- nrow(Y)
  q <- q
  obj$n <- n
  obj$d <- d
  obj$q <- q
  obj$E_X <- Matrix(rnorm(q * n), nrow = q, ncol = n)
  obj$E_XXt <- obj$E_X %*% t(obj$E_X)
  diag(obj$E_XXt) <- diag(obj$E_XXt) + obj$n
  obj$include_mu <- include_mu
  if (include_mu) {
    obj$E_mu <- rowSums(Y) / n
  } else {
    obj$E_mu <- 0
  }
  obj$E_W <-  Matrix(rnorm(d * q), nrow = d, ncol = q)
  obj$E_WtW <- t(obj$E_W) %*% obj$E_W
  diag(obj$E_WtW) <- diag(obj$E_WtW) + obj$d
  obj$alpha_MLE <- alpha_MLE
  obj$alpha_for_X <- alpha_for_X
  obj$a_alpha <- a_alpha
  obj$b_alpha <- b_alpha
  if (alpha_MLE) {
    obj$E_alpha <- d / diag(obj$E_WtW)
    obj$E_log_alpha <- log(obj$E_alpha)
    if (alpha_for_X) {
      obj$E_alpha_X <- d / diag(obj$E_XXt)
      obj$E_log_alpha_X <- log(obj$E_alpha_X)
    }
  } else {
    a <- a_alpha + d / 2
    b <- b_alpha + diag(obj$E_WtW) / 2
    obj$E_alpha <- a / b
    obj$E_log_alpha <- digamma(a) - log(b)
    if (alpha_for_X) {
      a <- a_alpha + d / 2
      b <- b_alpha + diag(obj$E_XXt) / 2
      obj$E_alpha_X <- a / b
      obj$E_log_alpha_X <- digamma(a) - log(b)
    }
  }
  obj$a_tau <- a_tau
  obj$b_tau <- b_tau
  obj$E_tau <- a_tau / b_tau
  obj$beta <- beta
  class(obj) <- c("vpca.data", "list")
  return(obj)
}


#' @export
#' @import Matrix
init_svd <- function(Y, q, alpha_MLE = FALSE, include_mu = TRUE, alpha_for_X = FALSE,
                        a_alpha = 0.001, b_alpha = 0.001,
                        a_tau = 0.001, b_tau = 0.001, beta = 0.001) {
  obj <- list(Y = Y, Y_Frob_Sq = sum(Y^2))
  svd_Y <- svd(Y, nu = q, nv = q)
  n <- ncol(Y)
  d <- nrow(Y)
  q <- q
  obj$n <- n
  obj$d <- d
  obj$q <- q
  obj$E_X <-  Matrix(t(svd_Y$v))
  obj$E_XXt <- obj$E_X %*% t(obj$E_X)
  diag(obj$E_XXt) <- diag(obj$E_XXt) + obj$n
  obj$include_mu <- include_mu
  if (include_mu) {
    obj$E_mu <- rowSums(Y) / n
  } else {
    obj$E_mu <- 0
  }
  obj$E_W <- Matrix(svd_Y$u %*% diag(svd_Y$d[1:q]))
  obj$E_WtW <- t(obj$E_W) %*% obj$E_W
  diag(obj$E_WtW) <- diag(obj$E_WtW) + obj$d
  obj$alpha_MLE <- alpha_MLE
  obj$alpha_for_X <- alpha_for_X
  obj$a_alpha <- a_alpha
  obj$b_alpha <- b_alpha
  if (alpha_MLE) {
    obj$E_alpha <- d / diag(obj$E_WtW)
    obj$E_log_alpha <- log(obj$E_alpha)
    if (alpha_for_X) {
      obj$E_alpha_X <- d / diag(obj$E_XXt)
      obj$E_log_alpha_X <- log(obj$E_alpha_X)
    }
  } else {
    a <- a_alpha + d / 2
    b <- b_alpha + diag(obj$E_WtW) / 2
    obj$E_alpha <- a / b
    obj$E_log_alpha <- digamma(a) - log(b)
    if (alpha_for_X) {
      a <- a_alpha + d / 2
      b <- b_alpha + diag(obj$E_XXt) / 2
      obj$E_alpha_X <- a / b
      obj$E_log_alpha_X <- digamma(a) - log(b)
    }
  }
  obj$a_tau <- a_tau
  obj$b_tau <- b_tau
  obj$E_tau <- a_tau / b_tau
  obj$beta <- beta
  class(obj) <- c("vpca.data", "list")
  return(obj)
}


#' @export
#' @import Matrix
init_from_flashier <- function(f, alpha_MLE = FALSE, include_mu = TRUE, alpha_for_X = FALSE,
                        a_alpha = 0.001, b_alpha = 0.001,
                        a_tau = 0.001, b_tau = 0.001, beta = 0.001) {
  if(inherits(f, "flash")) {
    f <- f[["flash_fit"]]
  }
  
  if (!inherits(f, "flash_fit")) {
    stop("f must be a flash or flash_fit object.")
  }

  obj <- list(Y = f$Y, Y_Frob_Sq = sum(f$Y^2))
  n <- ncol(f$Y)
  d <- nrow(f$Y)
  q <- ncol(f$EF[[1]])
  obj$n <- n
  obj$d <- d
  obj$q <- q
  obj$E_X <- t(Matrix(f$EF[[2]]))
  obj$E_XXt <- obj$E_X %*% t(obj$E_X)
  diag(obj$E_XXt) <- diag(obj$E_XXt) + obj$n
  obj$include_mu <- include_mu
  if (include_mu) {
    obj$E_mu <- rowSums(f$Y) / n
  } else {
    obj$E_mu <- 0
  }
  obj$E_W <-  Matrix(f$EF[[1]])
  obj$E_WtW <- t(obj$E_W) %*% obj$E_W
  diag(obj$E_WtW) <- diag(obj$E_WtW) + obj$d
  obj$alpha_MLE <- alpha_MLE
  obj$alpha_for_X <- alpha_for_X
  obj$a_alpha <- a_alpha
  obj$b_alpha <- b_alpha
  if (alpha_MLE) {
    obj$E_alpha <- d / diag(obj$E_WtW)
    obj$E_log_alpha <- log(obj$E_alpha)
    if (alpha_for_X) {
      obj$E_alpha_X <- d / diag(obj$E_XXt)
      obj$E_log_alpha_X <- log(obj$E_alpha_X)
    }
  } else {
    a <- a_alpha + d / 2
    b <- b_alpha + diag(obj$E_WtW) / 2
    obj$E_alpha <- a / b
    obj$E_log_alpha <- digamma(a) - log(b)
    if (alpha_for_X) {
      a <- a_alpha + d / 2
      b <- b_alpha + diag(obj$E_XXt) / 2
      obj$E_alpha_X <- a / b
      obj$E_log_alpha_X <- digamma(a) - log(b)
    }
  }
  obj$a_tau <- a_tau
  obj$b_tau <- b_tau
  obj$E_tau <- a_tau / b_tau
  obj$beta <- beta
  class(obj) <- c("vpca.data", "list")

  return(obj)
}