#' Update \eqn{Q(X)}.
#'
#' @param t d-by-n data matrix.
#' @param m.W \eqn{\mathbb{E}_Q[\mathbf{W}]}.
#' @param Sigma.W \eqn{\Sigma_\mathbf{W}}.
#' @param m.mu \eqn{\mathbb{E}_Q[\mu]}.
#' @param m.tau \eqn{\mathbb{E}_Q[\tau]}.
#' @section TODO: can inversion be avoided?
#'
#' @return \eqn{\mathbb{E}[\mathbf{x}_n]} and \eqn{\Sigma_\mathbf{x}}.
get.x <- function(t, m.W, Sigma.W, m.mu, m.tau, m.alpha) {
  # t: d by n
  # x: q by n
  # m.W: d by q
  d <- nrow(t)
  n <- ncol(t)
  q <- ncol(m.W)
  m.WtW <- nrow(t) * Sigma.W + t(m.W) %*% m.W
  if (missing(m.alpha)) {
    Sigma <- solve(m.tau * m.WtW + diag(q))
  }
  else {
    Sigma <- solve(m.tau * m.WtW + diag(m.alpha))
  }
  M <- m.tau * Sigma %*% t(m.W) %*% (t - outer(m.mu, rep(1, n))) #(t - m.mu)
  return(list(M, Sigma))
}

#' Update \eqn{Q(\mathbf{\mu})}.
#'
#' @param t d-by-n data matrix.
#' @param m.W \eqn{\mathbb{E}_Q[\mathbf{W}]}.
#' @param m.x \eqn{\mathbb{E}[\mathbf{x}_n]}.
#' @param m.tau \eqn{\mathbb{E}_Q[\tau]}.
#' @param beta \eqn{\beta}.
#'
#' @return \eqn{\mathbb{E}[\mu]} and \eqn{\Sigma_\mu}.
get.mu <- function(t, m.W, m.x, m.tau, beta) {
  d <- nrow(t)
  n <- ncol(t)
  m <- m.tau * rowSums(t - m.W %*% m.x) / (beta + n * m.tau)
  Sigma <- diag(d) / (beta + n * m.tau)
  return(list(m, Sigma))
}

#' Update \eqn{\mathbf{W}}.
#'
#' @param t d-by-n data matrix.
#' @param m.x \eqn{\mathbb{E}[\mathbf{x}_n]}.
#' @param Sigma.x \eqn{\Sigma_\mathbf{x}}.
#' @param m.mu \eqn{\mathbb{E}[\mu]}.
#' @param m.tau \eqn{\mathbb{E}_Q[\tau]}.
#' @param m.alpha \eqn{\mathbb{E}[\mathbf{\alpha}]}.
#' @section TODO: can inversion be avoided?
#'
#' @return \eqn{\mathbb{E}[\mathbf{W}]} and \eqn{\Sigma_{\mathbf{W}}}.
get.W <- function(t, m.x, Sigma.x, m.mu, m.tau, m.alpha) {
  d <- nrow(t)
  n <- ncol(t)
  q <- nrow(m.x)
  m.xxt <- n * Sigma.x + m.x %*% t(m.x)
  if (missing(m.alpha)) {
    Sigma <- solve(diag(q) + m.tau * m.xxt)
  } else {
    Sigma <- solve(diag(m.alpha) + m.tau * m.xxt)
  }
  M <- m.tau * Sigma %*% m.x %*% (t(t) - outer(rep(1, n), m.mu))
  M <- t(M)
  return(list(M, Sigma))
}

#' Update \eqn{Q(\mathbf{\alpha})}.
#'
#' @param m.W \eqn{\mathbb{E}[\mathbf{W}]}.
#' @param Sigma.W \eqn{\Sigma_{\mathbf{W}}}.
#' @param a.alpha \eqn{a_\alpha}.
#' @param b.alpha \eqn{b_\alpha}.
#'
#' @return \eqn{\tilde{a}_\alpha} and \eqn{\tilde{b}_\alpha}.
get.alpha <- function(m.W, Sigma.W, a.alpha, b.alpha) {
   d <- nrow(m.W)
   m.WtW <- d * Sigma.W + t(m.W) %*% m.W
   a <- a.alpha + d / 2
   b <- b.alpha + diag(m.WtW) / 2
   return(list(a, b))
}

#' Get MLE estimate of alpha
#'
#' @param m.W \eqn{\mathbb{E}[\mathbf{W}]}.
#' @param Sigma.W \eqn{\Sigma_{\mathbf{W}}}.
#'
#' @return \eqn{\hat{\alpha}_{\text{MLE}}} of \eqn{P(W \mid \alpha)}
#' @export
get.alpha.mle <- function(m.W, Sigma.W) {
  d <- nrow(m.W)
  m.WtW <- d * Sigma.W + t(m.W) %*% m.W
  return(d / diag(m.WtW))
}

#' Update \eqn{Q(\tau)}.
#'
#' @param t d-by-n data matrix.
#' @param m.x \eqn{\mathbb{E}[\mathbf{x}_n]}.
#' @param Sigma.x \eqn{\Sigma_\mathbf{x}}.
#' @param m.W \eqn{\mathbb{E}[\mathbf{W}]}.
#' @param Sigma.W \eqn{\Sigma_{\mathbf{W}}}.
#' @param m.mu \eqn{\mathbb{E}_Q[\mu]}.
#' @param Sigma.mu \eqn{\Sigma_{\mathbf{\mu}}}.
#' @param a.tau \eqn{a_\tau}.
#' @param b.tau \eqn{b_\tau}.
#'
#' @return \eqn{\tilde{a}_\tau} and \eqn{\tilde{b}_\tau}.
get.tau <- function(t, m.x, Sigma.x, m.W,
                    Sigma.W, m.mu, Sigma.mu, a.tau, b.tau) {
  d <- nrow(t)
  q <- nrow(m.x)
  n <- ncol(t)
  a <- a.tau + n * d / 2
  mu.norm2 <- sum(diag(Sigma.mu)) + sum(m.mu^2)
  m.WtW <- d * Sigma.W + t(m.W) %*% m.W
  m.xxt <- n * Sigma.x + m.x %*% t(m.x)
  b <- b.tau + (sum(t^2) + n * mu.norm2 + sum(diag(m.WtW %*% m.xxt)) +
                  2 * t(m.mu) %*% m.W %*% rowSums(m.x) -
                  2 * sum(diag(t(t) %*% m.W %*% m.x)) -
                  2 * t(m.mu) %*% rowSums(t)
                ) / 2
  return(list(a, b[1]))
}


#' Compute ELBO with given moments,
#'
#' @return ELBO \eqn{\mathcal{L}(\theta)}.
get.elbo <- function(t, m.x, Sigma.x, m.mu, Sigma.mu, m.W, Sigma.W,
                 alpha.param, a.tilde.tau, b.tilde.tau,
                 a.tau, b.tau, beta) {
  # complete likelihood
  n <- ncol(t)
  d <- nrow(t)
  q <- ncol(m.W)
  m.tau <- a.tilde.tau / b.tilde.tau
  E.log.tau <- digamma(a.tilde.tau) - log(b.tilde.tau)

  mu.norm2 <- sum(diag(Sigma.mu)) + sum(m.mu^2)
  m.WtW <- d * Sigma.W + t(m.W) %*% m.W
  m.xxt <- n * Sigma.x + m.x %*% t(m.x)
  complete.lik <- d * n / 2 * (E.log.tau - log(2*pi)) +
    -(sum(t^2) + n * mu.norm2 + sum(diag(m.WtW %*% m.xxt)) +
        2 * t(m.mu) %*% m.W %*% rowSums(m.x) -
        2 * sum(diag(t(t) %*% m.W %*% m.x)) -
        2 * t(m.mu) %*% rowSums(t)) * m.tau / 2 # related to tau update

  entropy <- n * q / 2 * log(2*pi*exp(1)) + n / 2 * log(det(Sigma.x)) + # H(X)
    d / 2 * log(2*pi*exp(1)) + 0.5 * log(det(Sigma.mu)) + # H(mu)
    d * q / 2 * log(2*pi*exp(1)) + d / 2 * log(det(Sigma.W)) + # H(W)
    a.tilde.tau - log(b.tilde.tau) + lgamma(a.tilde.tau) + (1 - a.tilde.tau) * digamma(a.tilde.tau)# H(tau)

  prior.mu <- -d/2 * log(2*pi) + d/2 * log(beta) - beta/2 * (sum(diag(Sigma.mu)) + sum(m.mu^2))
  prior.tau <- a.tilde.tau * log(b.tilde.tau) - lgamma(a.tilde.tau) +
    (a.tilde.tau - 1) * E.log.tau - b.tau * m.tau
  prior <- prior.mu + prior.tau

  # m.alpha <- alpha.param$W$m
  if (!alpha.param$mle) {
    if ("W" %in% names(alpha.param)) {
      m.alpha <- alpha.param$W$m
      a.tilde.alpha <- alpha.param$W$a.tilde
      b.tilde.alpha <- alpha.param$W$b.tilde
      a.alpha <- alpha.param$W$a
      b.alpha <- alpha.param$W$b
      E.log.alpha <- digamma(a.tilde.alpha) - log(b.tilde.alpha)
      prior.W <- d / 2 * sum(E.log.alpha) - d*q/2 * log(2*pi) - 0.5 * sum(diag(m.WtW) * m.alpha)
      prior.alpha <- q * a.alpha * log(b.alpha) + (a.alpha - 1) * sum(E.log.alpha) - b.alpha * sum(m.alpha) - q * lgamma(a.alpha)
      prior <- prior + prior.W + prior.alpha

      H.alpha <- q * (a.tilde.alpha + lgamma(a.tilde.alpha) + (1 - a.tilde.alpha) * digamma(a.tilde.alpha)) - sum(log(b.tilde.alpha)) # H(alpha)
      entropy <- entropy + H.alpha
    }
    if ("x" %in% names(alpha.param)) {
      m.alpha <- alpha.param$x$m
      a.tilde.alpha <- alpha.param$x$a.tilde
      b.tilde.alpha <- alpha.param$x$b.tilde
      a.alpha <- alpha.param$x$a
      b.alpha <- alpha.param$x$b
      E.log.alpha <- digamma(a.tilde.alpha) - log(b.tilde.alpha)
      prior.x <- n / 2 * sum(E.log.alpha) - n*d/2 * log(2*pi) - 0.5 * sum(diag(m.xxt) * m.alpha)
      prior.alpha <- d * a.alpha * log(b.alpha) + (a.alpha - 1) * sum(E.log.alpha) - b.alpha * sum(m.alpha) - d * lgamma(a.alpha)
      prior <- prior + prior.x + prior.alpha
      H.alpha <- d * (a.tilde.alpha + lgamma(a.tilde.alpha) + (1 - a.tilde.alpha) * digamma(a.tilde.alpha)) - sum(log(b.tilde.alpha)) # H(alpha)
      entropy <- entropy + H.alpha
    }
  } else {
    if ("W" %in% names(alpha.param)) {
      m.alpha <- alpha.param$W$m
      prior.W <- d / 2 * sum(log(m.alpha)) - d*q/2 * log(2*pi) - 0.5 * sum(diag(m.WtW) * m.alpha)
      prior <- prior + prior.W
    }
    if ("x" %in% names(alpha.param)) {
      m.alpha <- alpha.param$x$m
      prior.x <- n / 2 * sum(log(m.alpha)) - n*d/2 * log(2*pi) - 0.5 * sum(diag(m.WtW) * m.alpha)
      prior <- prior + prior.x
    }
  }

  return(complete.lik + prior + entropy)
}


#' Variational PCA algorithm runner.
#'
#' @param t d-by-n data matrix.
#' @param fit.alpha.mle Whether to fit \eqn{\alpha} via MLE.
#' @param W.alpha Whether to fit \eqn{\alpha} for \eqn{W}.
#' @param x.alpha Whether to fit \eqn{\alpha} for \eqn{x}.
#' @param a.alpha \eqn{a_\alpha}.
#' @param b.alpha \eqn{b_\alpha}.
#' @param a.tau \eqn{a_\tau}.
#' @param b.tau \eqn{a_\alpha}.
#' @param beta \eqn{\beta}.
#' @param min.iter Minimum number of iterations.
#' @param max.iter Maximum number of iterations.
#' @param epsilon The algorithm terminates if the difference in ELBO is less than \eqn{\epsilon}.
#'
#' @return \eqn{\mathbb{E}[W]} and ELBOs in all iterations.
#' @export
vpca <- function(t, fit.alpha.mle = F, W.alpha = T, x.alpha = F, a.alpha = 0.001, b.alpha = 0.001,
                 a.tau = 0.001, b.tau = 0.001, beta = 0.001,
                 min.iter = 1000, max.iter = 3000, epsilon = 0.0001) {
  if (!W.alpha & !x.alpha) {
    stop("Must fit alpha for at least one of W and alpha.")
  }
  d <- nrow(t)
  n <- ncol(t)
  q <- d - 1

  # initialization
  m.x <- matrix(rnorm(q*n), nrow = q, ncol = n)
  Sigma.x <- diag(q)
  m.mu <- rowSums(t) / n
  Sigma.mu <- diag(d)
  m.W <-  matrix(rnorm(d * q), nrow = d, ncol = q)
  Sigma.W <- diag(q)
  # svd.t <- svd((t - m.mu) %*% t(t - m.mu))
  # m.W <- svd.t$u[,1:q] %*% diag(sqrt(svd.t$d[1:q]))
  if (fit.alpha.mle) {
    alpha.param <- list(mle=T)
    if (W.alpha) {
      W.alpha.param <- list(m=rep(a.alpha/b.alpha, q))
      alpha.param$W <- W.alpha.param
    }
    if (x.alpha) {
      x.alpha.param <- list(m=rep(a.alpha/b.alpha, q))
      alpha.param$x <- x.alpha.param
    }
  } else {
    alpha.param <- list(mle=F)
    W.alpha.param <- list(m=a.alpha / rep(b.alpha, q), a.tilde=a.alpha, b.tilde=rep(b.alpha, q),
                        a=a.alpha, b=b.alpha)
    if (W.alpha) {
      W.alpha.param <- list(m=a.alpha / rep(b.alpha, q), a.tilde=a.alpha, b.tilde=rep(b.alpha, q),
                            a=a.alpha, b=b.alpha)
      alpha.param$W <- W.alpha.param
    }
    if (x.alpha) {
      x.alpha.param <- list(m=a.alpha / rep(b.alpha, q), a.tilde=a.alpha, b.tilde=rep(b.alpha, q),
                            a=a.alpha, b=b.alpha)
      alpha.param$x <- x.alpha.param
    }
  }

  a.tilde.tau <- a.tau
  b.tilde.tau <- b.tau
  m.tau <- a.tilde.tau / b.tilde.tau

  elbo.old <- get.elbo(t, m.x, Sigma.x, m.mu, Sigma.mu, m.W, Sigma.W,
                       alpha.param, a.tilde.tau, b.tilde.tau,
                       a.tau, b.tau, beta)[1]
  elbo <- c(elbo.old)

  count <- 0
  while (T) {
    if (x.alpha) {
      x.res <- get.x(t, m.W, Sigma.W, m.mu, m.tau, alpha.param$x$m)
    }
    else {
      x.res <- get.x(t, m.W, Sigma.W, m.mu, m.tau)
    }
    m.x.new <- x.res[[1]]; Sigma.x.new <- x.res[[2]]
    mu.res <- get.mu(t, m.W, m.x.new, m.tau, beta)
    m.mu.new <- mu.res[[1]]; Sigma.mu.new <- mu.res[[2]]

    if (W.alpha) {
      W.res <- get.W(t, m.x.new, Sigma.x.new, m.mu.new, m.tau, alpha.param$W$m)
    }
    else {
      W.res <- get.W(t, m.x.new, Sigma.x.new, m.mu.new, m.tau)
    }
    m.W.new <- W.res[[1]]; Sigma.W.new <- W.res[[2]]

    if (fit.alpha.mle) {
      if (W.alpha) {
        alpha.res <- get.alpha.mle(m.W.new, Sigma.W.new)
        alpha.param$W$m <- alpha.res
      }
      if (x.alpha) {
        alpha.res <- get.alpha.mle(t(m.x.new), Sigma.x.new)
        alpha.param$x$m <- alpha.res
      }
    } else {
      if (W.alpha) {
        alpha.res <- get.alpha(m.W.new, Sigma.W.new, a.alpha, b.alpha)
        a.tilde.alpha.new <- alpha.res[[1]]; b.tilde.alpha.new <- alpha.res[[2]]
        m.alpha.new <- a.tilde.alpha.new / b.tilde.alpha.new
        alpha.param$W$m <- m.alpha.new
        alpha.param$W$a.tilde <- a.tilde.alpha.new
        alpha.param$W$b.tilde <- b.tilde.alpha.new
      }
      if (x.alpha) {
        alpha.res <- get.alpha(t(m.x.new), Sigma.x.new, a.alpha, b.alpha)
        a.tilde.alpha.new <- alpha.res[[1]]; b.tilde.alpha.new <- alpha.res[[2]]
        m.alpha.new <- a.tilde.alpha.new / b.tilde.alpha.new
        alpha.param$x$m <- m.alpha.new
        alpha.param$x$a.tilde <- a.tilde.alpha.new
        alpha.param$x$b.tilde <- b.tilde.alpha.new
      }
    }


    tau.res <- get.tau(t, m.x.new, Sigma.x.new, m.W.new,
                       Sigma.W.new, m.mu.new, Sigma.mu.new, a.tau, b.tau)
    a.tilde.tau.new <- tau.res[[1]]; b.tilde.tau.new <- tau.res[[2]]
    m.tau.new <- a.tilde.tau.new / b.tilde.tau.new

    m.x <- m.x.new
    Sigma.x <- Sigma.x.new
    m.mu <- m.mu.new
    Sigma.mu <- Sigma.mu.new
    m.W <- m.W.new
    Sigma.W <- Sigma.W.new
    # a.tilde.alpha <- a.tilde.alpha.new
    # b.tilde.alpha <- b.tilde.alpha.new
    # m.alpha <- m.alpha.new
    a.tilde.tau <- a.tilde.tau.new
    b.tilde.tau <- b.tilde.tau.new
    m.tau <- m.tau.new

    elbo.new <- get.elbo(t, m.x, Sigma.x, m.mu, Sigma.mu, m.W, Sigma.W,
                         alpha.param, a.tilde.tau, b.tilde.tau,
                         a.tau, b.tau, beta)[1]
    count <- count + 1
    diff <- abs(elbo.new - elbo.old)
    elbo.old <- elbo.new
    elbo <- c(elbo, elbo.new)
    if (diff <= epsilon & count >= min.iter) {
      break
    }
    if (count == max.iter) {
      print("Maximum number of iteration reached.")
      break
    }
  }
  return(list(m.W, Sigma.W, m.x, Sigma.x, m.tau, m.mu, alpha.param, elbo))
}
