#' Title Update \eqn{Q(X)}.
#'
#' @param t d-by-n data matrix.
#' @param m.W \eqn{\mathbb{E}_Q[\mathbf{W}]}.
#' @param Sigma.W \eqn{\Sigma_\mathbf{W}}.
#' @param m.mu \eqn{\mathbb{E}_Q[\mu]}.
#' @param m.tau \eqn{\mathbb{E}_Q[\tau]}.
#' @section TODO: can inversion be avoided?
#'
#' @return \eqn{\mathbb{E}[\mathbf{x}_n]} and \eqn{\Sigma_\mathbf{x}}.
get.x <- function(t, m.W, Sigma.W, m.mu, m.tau) {
  # t: d by n
  # x: q by n
  # m.W: d by q
  M <- matrix(0, nrow = ncol(m.W), ncol = ncol(t))
  m.WtW <- nrow(t) * Sigma.W + t(m.W) %*% m.W
  Sigma.inv <- m.tau * m.WtW + diag(nrow(m.WtW))
  for (i in 1:ncol(t)) {
    M[, i] <- m.tau * solve(Sigma.inv, t(m.W) %*% (t[, i] - m.mu))
  }
  return(list(M, solve(Sigma.inv)))
}

#' Title Update \eqn{Q(\mathbf{\mu})}.
#'
#' @param t d-by-n data matrix.
#' @param m.W \eqn{\mathbb{E}_Q[\mathbf{W}]}.
#' @param m.x \eqn{\mathbb{E}[\mathbf{x}_n]}.
#' @param m.tau \eqn{\mathbb{E}_Q[\tau]}.
#' @param beta \eqn{\beta}.
#'
#' @return \eqn{\mathbb{E}[\mu]} and \eqn{\Sigma_\mu}.
get.mu <- function(t, m.W, m.x, m.tau, beta) {
  m <- m.tau * rowSums(t - m.W %*% m.x) / (beta + ncol(t) * m.tau)
  Sigma <- diag(nrow(t)) / (beta + ncol(t) * m.tau)
  return(list(m, Sigma))
}

#' Title Update \eqn{\mathbf{W}}.
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
  # M: d by q
  M <- matrix(0, nrow = nrow(t), ncol = nrow(m.x))
  m.xxt <- ncol(t) * Sigma.x + m.x %*% t(m.x)
  Sigma.inv <- diag(m.alpha) + m.tau * m.xxt
  for (k in 1:nrow(t)) {
    for (n in 1:ncol(t)) {
      M[k, ] <- M[k, ] + (t[k, n] - m.mu[k]) * m.x[, n]
    }
    M[k, ] <- m.tau * solve(Sigma.inv, M[k, ])
  }
  return(list(M, solve(Sigma.inv)))
}

#' Title Update \eqn{Q(\mathbf{\alpha})}.
#'
#' @param m.W \eqn{\mathbb{E}[\mathbf{W}]}.
#' @param Sigma.W \eqn{\Sigma_{\mathbf{W}}}.
#' @param a.alpha \eqn{a_\alpha}.
#' @param b.alpha \eqn{b_\alpha}.
#'
#' @return \eqn{\tilde{a}_\alpha} and \eqn{\tilde{b}_\alpha}.
get.alpha <- function(m.W, Sigma.W, a.alpha, b.alpha) {
   a <- a.alpha + nrow(m.W) / 2
   tr.Sigma.W <- sum(diag(Sigma.W))
   b <- rep(0, ncol(m.W))
   for (i in 1:ncol(m.W)) {
     b[i] <- b.alpha + (nrow(m.W) * Sigma.W[i, i] + sum(m.W[, i]^2)) / 2
   }
   # b <- b.alpha + (colSums(m.W^2) + tr.Sigma.W) / 2
   return(list(a, b))
}

#' Title Update \eqn{Q(\tau)}.
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
  a <- a.tau + nrow(t) * ncol(t) / 2
  mu.norm2 <- sum(diag(Sigma.mu)) + sum(m.mu^2)
  m.WtW <- nrow(t) * Sigma.W + t(m.W) %*% m.W
  m.xxt <- ncol(t) * Sigma.x + m.x %*% t(m.x)
  b <- b.tau + (sum(t^2) + ncol(t) * mu.norm2 + sum(diag(m.WtW %*% m.xxt)) +
                  2 * t(m.mu) %*% m.W %*% rowSums(m.x) -
                  2 * sum(diag(t(t) %*% m.W %*% m.x)) -
                  2 * t(m.mu) %*% rowSums(t)
                ) / 2
  return(list(a, b[1]))
}


#' Title Compute ELBO with given moments,
#'
#' @return ELBO \eqn{\mathcal{L}(\theta)}.
get.elbo <- function(t, m.x, Sigma.x, m.mu, Sigma.mu, m.W, Sigma.W,
                 a.tilde.alpha, b.tilde.alpha, a.tilde.tau, b.tilde.tau,
                 a.alpha, b.alpha, a.tau, b.tau, beta) {
  # complete likelihood
  n <- ncol(t)
  d <- nrow(t)
  q <- ncol(m.W)
  m.alpha <- a.tilde.alpha / b.tilde.alpha
  m.tau <- a.tilde.tau / b.tilde.tau
  E.log.alpha <- digamma(a.tilde.alpha) - log(b.tilde.alpha)
  E.log.tau <- digamma(a.tilde.tau) - log(b.tilde.tau)

  complete.lik <- d * n / 2 * (E.log.tau - log(2*pi)) +
    m.tau * (sum(diag(t(t) %*% m.W %*% m.x)) + t(m.mu) %*% rowSums(t) -
               0.5 * sum(t^2) - 0.5 * n * (sum(diag(Sigma.mu)) + sum(m.mu^2)))

  entropy <- n * q / 2 * log(2*pi*exp(1)) + n / 2 * log(det(Sigma.x)) + # H(X)
    d / 2 * log(2*pi*exp(1)) + 0.5 * log(det(Sigma.mu)) + # H(mu)
    d * q / 2 * log(2*pi*exp(1)) + q / 2 * log(det(Sigma.W)) + # H(W)
    q * (a.tilde.alpha + log(gamma(a.tilde.alpha)) + (1 - a.tilde.alpha) * digamma(a.tilde.alpha)) - sum(log(b.tilde.alpha)) + # H(alpha)
    a.tilde.tau - log(b.tilde.tau) + log(gamma(a.tilde.tau)) + (1 - a.tilde.tau) * digamma(a.tilde.tau) # H(tau)

  prior.W.alpha <- sum(E.log.alpha) + d / 2 * sum(E.log.alpha) - d*q/2 * log(2*pi) -
    0.5 * sum(diag(Sigma.W)) * sum(m.alpha) - 0.5 * sum(m.alpha * colSums(m.W^2))
  prior.x <- -n*q/2 * log(2*pi) - 0.5 * n * sum(diag(Sigma.x)) - 0.5 * sum(m.x^2)
  prior.mu <- -d/2 * log(2*pi) + d/2 * log(beta) - beta/2 * (sum(diag(Sigma.mu)) + sum(m.mu^2))
  prior.tau <- a.tilde.tau * log(b.tilde.tau) - log(gamma(a.tilde.tau)) +
    (a.tilde.tau - 1) * E.log.tau - b.tau * m.tau

  prior <- prior.W.alpha + prior.x + prior.mu + prior.tau

  return(complete.lik + prior - entropy)
}


#' Title Variational PCA algorithm runner.
#'
#' @param t d-by-n data matrix.
#' @param a.alpha Hyperparameter \eqn{a_\alpha}.
#' @param b.alpha \eqn{b_\alpha}.
#' @param a.tau \eqn{a_\tau}.
#' @param b.tau \eqn{a_\alpha}.
#' @param beta \eqn{\beta}.
#' @param max.iter Maximum number of iterations.
#' @param epsilon The algorithm terminates if the difference in ELBO is less than \eqn{\epsilon}.
#'
#' @return \eqn{\mathbb{E}[W]} and ELBOs in all iterations.
#' @export
vpca <- function(t, a.alpha = 0.001, b.alpha = 0.001,
                 a.tau = 0.001, b.tau = 0.001, beta = 0.001,
                 max.iter = 1000, epsilon = 0.001) {
  d <- nrow(t)
  n <- ncol(t)
  q <- d - 1

  # initialization
  m.x <- matrix(rnorm(q * n), nrow = q, ncol = n)
  Sigma.x <- diag(q)
  m.mu <- rnorm(d)
  Sigma.mu <- diag(d)
  m.W <-matrix(rnorm(d * q), nrow = d, ncol = q)
  Sigma.W <- diag(q)
  a.tilde.alpha <- a.alpha
  b.tilde.alpha <- rep(b.alpha, q)
  m.alpha <- a.tilde.alpha / b.tilde.alpha
  a.tilde.tau <- a.tau
  b.tilde.tau <- b.tau
  m.tau <- a.tilde.tau / b.tilde.tau

  elbo.old <- get.elbo(t, m.x, Sigma.x, m.mu, Sigma.mu, m.W, Sigma.W,
                       a.tilde.alpha, b.tilde.alpha, a.tilde.tau, b.tilde.tau,
                       a.alpha, b.alpha, a.tau, b.tau, beta)[1]
  elbo <- c(elbo.old)

  count <- 0
  while (T) {
    x.res <- get.x(t, m.W, Sigma.W, m.mu, m.tau)
    m.x.new <- x.res[[1]]; Sigma.x.new <- x.res[[2]]
    mu.res <- get.mu(t, m.W, m.x.new, m.tau, beta)
    m.mu.new <- mu.res[[1]]; Sigma.mu.new <- mu.res[[2]]
    W.res <- get.W(t, m.x.new, Sigma.x.new, m.mu.new, m.tau, m.alpha)
    m.W.new <- W.res[[1]]; Sigma.W.new <- W.res[[2]]
    alpha.res <- get.alpha(m.W.new, Sigma.W.new, a.alpha, b.alpha)
    a.tilde.alpha.new <- alpha.res[[1]]; b.tilde.alpha.new <- alpha.res[[2]]
    m.alpha.new <- a.tilde.alpha.new / b.tilde.alpha.new
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
    m.alpha <- m.alpha.new
    m.tau <- m.tau.new

    elbo.new <- get.elbo(t, m.x, Sigma.x, m.mu, Sigma.mu, m.W, Sigma.W,
                         a.tilde.alpha, b.tilde.alpha, a.tilde.tau, b.tilde.tau,
                         a.alpha, b.alpha, a.tau, b.tau, beta)[1]

    count <- count + 1
    diff <- abs(elbo.new - elbo.old)
    elbo.old <- elbo.new
    elbo <- c(elbo, elbo.new)
    if (diff < epsilon || count == max.iter) {
      if (count == max.iter) {
        warning("Maximum number of iteration reached.")
      }
      break
    }
  }
  return(list(m.W, elbo))
}
