% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vpca.R
\name{vpca}
\alias{vpca}
\title{Title Variational PCA algorithm runner.}
\usage{
vpca(
  t,
  a.alpha = 0.001,
  b.alpha = 0.001,
  a.tau = 0.001,
  b.tau = 0.001,
  beta = 0.001,
  max.iter = 1000,
  epsilon = 0.001
)
}
\arguments{
\item{t}{d-by-n data matrix.}

\item{a.alpha}{Hyperparameter \eqn{a_\alpha}.}

\item{b.alpha}{\eqn{b_\alpha}.}

\item{a.tau}{\eqn{a_\tau}.}

\item{b.tau}{\eqn{a_\alpha}.}

\item{beta}{\eqn{\beta}.}

\item{max.iter}{Maximum number of iterations.}

\item{epsilon}{The algorithm terminates if the difference in ELBO is less than \eqn{\epsilon}.}
}
\value{
\eqn{\mathbb{E}[W]} and ELBOs in all iterations.
}
\description{
Title Variational PCA algorithm runner.
}
