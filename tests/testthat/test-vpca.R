d <- 3
n <- 10
q <- 2
t <- matrix(1, d, n)
beta <- 0.001

test_that("Variational update for X", {
  m.W <- matrix(1, d, q)
  Sigma.W <- diag(q)
  m.mu <- rep(0, d)
  m.tau <- 1
  res <- get.x(t, m.W, Sigma.W, m.mu, m.tau)
  expected.Sigma <- solve(diag(2)+3*diag(2)+t(m.W)%*%m.W)

  expect_equal(res[[2]], expected.Sigma, tolerance = 0.01)
  expect_equal(res[[1]], m.tau*expected.Sigma%*%t(m.W)%*%t, tolerance = 0.01)
})

test_that("Variational update for mu", {
  m.W <- matrix(1, d, q)
  m.x <- matrix(1, q, n)
  Sigma.W <- diag(q)
  m.tau <- 1
  res <- get.mu(t, m.W, m.x, m.tau, beta)
  expected.Sigma <- diag(d) / (beta + n * m.tau)
  expected.m <- m.tau * rowSums(t - m.W %*% m.x) / (beta + n * m.tau)
  expect_equal(res[[2]], expected.Sigma, tolerance = 0.01)
  expect_equal(res[[1]], expected.m, tolerance = 0.01)
})

test_that("Variational update for W", {
  m.W <- matrix(1, d, q)
  m.x <- matrix(1, q, n)
  Sigma.W <- diag(q)
  m.tau <- 1
  res <- get.mu(t, m.x, Sigma.x, m.mu, m.tau, m.alpha)
  expected.Sigma <- diag(d) / (beta + n * m.tau)
  expected.m <- m.tau * rowSums(t - m.W %*% m.x) / (beta + n * m.tau)
  expect_equal(res[[2]], expected.Sigma, tolerance = 0.01)
  expect_equal(res[[1]], expected.m, tolerance = 0.01)
})
